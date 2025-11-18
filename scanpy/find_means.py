"""
find_means.py

Scanpy-friendly utilities to inspect distributions, detect modes (local maxima),
and propose cut thresholds at local minima between modes.

Supports two use-cases:
1) Gene expression per cell (Seurat-like functionality)
2) QC metrics from adata.obs such as total_counts, pct_counts_mt/perc_mito,
   n_genes_by_counts, doublet_score, etc.

Features
- Works with AnnData and sparse matrices for gene expression.
- For a given vector (gene expression or obs metric), plots a histogram of cells
  and a KDE curve, detects modes and valleys, and proposes multiple thresholds.
- Human-in-the-loop: proposes several threshold candidates; a person selects one.
- Flexible threshold application for QC metrics (>, <, between, outside).

Dependencies: scanpy, anndata, numpy, scipy, matplotlib, seaborn, scikit-learn

Example (QC metric)
-------------------
import scanpy as sc
from scanpy.find_means import (
    suggest_thresholds_for_obs,
    interactive_select_threshold,
    apply_threshold_to_obs_metric,
)

adata = sc.read_h5ad("your_data.h5ad")
res = suggest_thresholds_for_obs(
    adata,
    key="total_counts",
    transform="log1p",       # often helpful for skewed metrics
    bins=100,
    min_peak_prominence=0.01,
    min_peak_distance=0.05,
    kde_bandwidth=None,
    show=True,
    save="scanpy/figures/total_counts_thresholds.png",
)
choice = interactive_select_threshold(res["thresholds"])  # human selects
if choice is not None:
    # For total_counts you might want to keep cells with counts above a floor
    apply_threshold_to_obs_metric(
        adata,
        key="total_counts",
        threshold=choice,
        direction=">=",  # mark cells that pass the QC
        key_added="qc_keep_total_counts",
        original_scale=res["original_scale"],  # converts back if transform was used
    )

Example (Gene expression)
-------------------------
from scanpy.find_means import suggest_thresholds_for_gene, apply_threshold_to_obs
res = suggest_thresholds_for_gene(adata, gene="Gad1", bins=100, show=True)
thr = interactive_select_threshold(res["thresholds"])  # choose a cutoff
if thr is not None:
    apply_threshold_to_obs(adata, "Gad1", thr, key_added="Gad1_positive")
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple, Dict, Any, Union

import numpy as np
import scipy.sparse as sp
from scipy.signal import find_peaks, argrelextrema
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import scanpy as sc  # type: ignore
    from anndata import AnnData  # type: ignore
except Exception:  # pragma: no cover - allows import without scanpy during static analysis
    AnnData = Any  # type: ignore


@dataclass
class ThresholdSuggestion:
    name: str
    peaks: List[float]
    valleys: List[float]
    thresholds: List[float]
    histogram_bins: np.ndarray
    histogram_counts: np.ndarray
    kde_x: np.ndarray
    kde_y: np.ndarray
    transform: str = "none"
    original_scale: Optional[Dict[str, float]] = None  # for inverse transforms


# ----------------------
# Internal helper methods
# ----------------------

def _prepare_hist(values: np.ndarray, bins: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        edges = np.linspace(0.0, 1.0, bins)
        counts, _ = np.histogram([], bins=edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        width = np.diff(edges)
        return counts, edges, centers, width
    vmin, vmax = float(np.min(vals)), float(np.max(vals))
    if vmax == vmin:
        edges = np.linspace(vmin, vmax + 1e-8, max(5, bins))
    else:
        edges = np.linspace(vmin, vmax, bins)
    counts, edges = np.histogram(vals, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = np.diff(edges)
    return counts, edges, centers, width


def _fit_kde(values: np.ndarray, bandwidth: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Fit a 1D KDE on positive values and return grid x and density y.

    For very sparse distributions with many zeros, omit zeros from KDE fitting.
    """
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 0.0])

    vmin, vmax = np.min(vals), np.max(vals)
    if vmin == vmax:
        x = np.array([vmin, vmax + 1e-8])
        y = np.array([1.0, 1.0])
        return x, y

    pos = vals[vals > 0]
    # If few positive values, use a smoothed histogram-like curve
    if pos.size < 5:
        x = np.linspace(vmin, vmax, 200)
        hist, edges = np.histogram(vals, bins=min(20, max(5, int(np.sqrt(vals.size)))) , density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        y = np.interp(x, centers, hist, left=0, right=0)
        return x, y

    if bandwidth is None:
        std = np.std(pos)
        n = pos.size
        bw = (n ** (-1.0 / 5.0)) * (std + 1e-12)  # Scott-like
        if not np.isfinite(bw) or bw <= 0:
            bw = max(1e-3, 0.1 * (np.max(pos) - np.min(pos)))
    else:
        bw = float(bandwidth)

    kde = KernelDensity(kernel="gaussian", bandwidth=bw)
    kde.fit(pos.reshape(-1, 1))

    x = np.linspace(np.min(pos), np.max(pos), 512)
    log_density = kde.score_samples(x.reshape(-1, 1))
    y = np.exp(log_density)
    return x, y


def _find_local_extrema(y: np.ndarray, min_peak_prominence: float, min_peak_distance: float) -> Tuple[np.ndarray, np.ndarray]:
    if y.size < 5:
        return np.array([], dtype=int), np.array([], dtype=int)
    ymax = float(np.max(y))
    if ymax <= 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    prominence = max(1e-12, float(min_peak_prominence) * ymax)
    distance = max(1, int(float(min_peak_distance) * y.size))
    peaks, _ = find_peaks(y, prominence=prominence, distance=distance)
    minima = argrelextrema(y, np.less, order=max(1, distance // 2))[0]
    return peaks, minima


def _maybe_inverse_transform(values: Union[float, np.ndarray], transform: str, original_scale: Optional[Dict[str, float]]):
    if transform == "none" or original_scale is None:
        return values
    if transform == "log1p":
        # original_scale can contain a multiplier, e.g., factor if needed
        if isinstance(values, np.ndarray):
            return np.expm1(values)
        else:
            return float(np.expm1(values))
    raise ValueError(f"Unknown transform '{transform}'")


# ----------------------------
# Core vector-based suggestion
# ----------------------------

def suggest_thresholds_for_vector(
    values: np.ndarray,
    name: str,
    bins: int = 100,
    zero_inflated: bool = False,
    min_peak_prominence: float = 0.01,
    min_peak_distance: float = 0.05,
    kde_bandwidth: Optional[float] = None,
    transform: str = "none",
    show: bool = True,
    save: Optional[str] = None,
) -> Dict[str, Any]:
    """Analyze a 1D array (one value per cell) and propose threshold candidates.

    - transform: 'none' or 'log1p'. The histogram/KDE/thresholds are computed on
      the transformed scale (if any). The return includes "original_scale" for
      converting thresholds back if needed.
    """
    vals = values.astype(float)
    finite_mask = np.isfinite(vals)
    vals = vals[finite_mask]

    if transform not in {"none", "log1p"}:
        raise ValueError("transform must be one of {'none','log1p'}")

    if transform == "log1p":
        vals_t = np.log1p(np.clip(vals, a_min=0, a_max=None))
    else:
        vals_t = vals

    # Histogram and KDE on transformed values
    hist_counts, hist_edges, centers, width = _prepare_hist(vals_t, bins=bins)
    x_kde, y_kde = _fit_kde(vals_t, bandwidth=kde_bandwidth)

    # Extrema on KDE
    peaks_idx, minima_idx = _find_local_extrema(y_kde, min_peak_prominence, min_peak_distance)
    peaks_x = x_kde[peaks_idx] if peaks_idx.size else np.array([])
    minima_x = x_kde[minima_idx] if minima_idx.size else np.array([])

    thresholds_t: List[float] = []
    if minima_idx.size > 0:
        thresholds_t = sorted(float(x) for x in minima_x)

    # Optionally include zero as threshold for zero-inflated metrics
    if zero_inflated and np.any(vals_t == 0) and np.any(vals_t > 0):
        thresholds_t = sorted(set(thresholds_t + [0.0]))

    # Plot
    if show or (save is not None):
        plt.figure(figsize=(8, 5))
        sns.set(style="whitegrid")
        plt.bar(centers, hist_counts, width=width, color="#c7d5eb", edgecolor="#5072a7", label="Cells per bin")
        if x_kde.size > 1:
            y_kde_scaled = y_kde * (np.sum(hist_counts) * np.mean(width))
            plt.plot(x_kde, y_kde_scaled, color="#d62728", lw=2, label="KDE (scaled)")
        if peaks_idx.size > 0:
            for px in peaks_x:
                plt.axvline(px, color="#ff7f0e", ls="--", lw=1.5, alpha=0.8)
        if len(thresholds_t) > 0:
            for t in thresholds_t:
                plt.axvline(t, color="#2ca02c", ls="-", lw=2, alpha=0.9)
        thr_text = ", ".join(f"{t:.3g}" for t in thresholds_t) if thresholds_t else "None"
        peak_text = ", ".join(f"{p:.3g}" for p in peaks_x) if peaks_idx.size else "None"
        tlabel = f" (transform={transform})" if transform != "none" else ""
        plt.title(f"{name}{tlabel}: histogram and KDE with suggested thresholds\nPeaks: {peak_text} | Thresholds: {thr_text}")
        plt.xlabel("Value" + (" (log1p)" if transform == "log1p" else ""))
        plt.ylabel("# Cells (per bin)")
        plt.legend(loc="best")
        plt.tight_layout()
        if save is not None:
            plt.savefig(save, dpi=200, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close()

    return {
        "name": name,
        "peaks": [float(x) for x in peaks_x],
        "valleys": [float(x) for x in minima_x],
        "thresholds": thresholds_t,
        "histogram_bins": hist_edges,
        "histogram_counts": hist_counts,
        "kde_x": x_kde,
        "kde_y": y_kde,
        "transform": transform,
        "original_scale": {"transform": transform},
    }


# -----------------------------------
# Public API for OBS (QC) metric usage
# -----------------------------------

def suggest_thresholds_for_obs(
    adata: AnnData,
    key: str,
    bins: int = 100,
    zero_inflated: bool = False,
    min_peak_prominence: float = 0.01,
    min_peak_distance: float = 0.05,
    kde_bandwidth: Optional[float] = None,
    transform: str = "none",
    show: bool = True,
    save: Optional[str] = None,
) -> Dict[str, Any]:
    """Suggest thresholds for a QC metric stored in adata.obs[key].

    - transform: 'none' or 'log1p'. Recommended to use 'log1p' for heavy-tailed
      metrics like total_counts or n_genes_by_counts.
    """
    if key not in adata.obs:
        raise KeyError(f"obs['{key}'] not found")
    col = adata.obs[key].to_numpy()
    return suggest_thresholds_for_vector(
        values=col,
        name=key,
        bins=bins,
        zero_inflated=zero_inflated,
        min_peak_prominence=min_peak_prominence,
        min_peak_distance=min_peak_distance,
        kde_bandwidth=kde_bandwidth,
        transform=transform,
        show=show,
        save=save,
    )


def apply_threshold_to_obs_metric(
    adata: AnnData,
    key: str,
    threshold: Union[float, Tuple[float, float]],
    direction: str = ">=",
    key_added: Optional[str] = None,
    original_scale: Optional[Dict[str, float]] = None,
) -> np.ndarray:
    """Apply a threshold to an obs metric and store a boolean mask in adata.obs.

    Parameters
    - key: obs column name (e.g., 'total_counts', 'pct_counts_mt')
    - threshold: a float if using '>', '>=', '<', '<='; a (low, high) tuple when
      using 'between' or 'outside'
    - direction: one of {'>', '>=', '<', '<=', 'between', 'outside'}
      - 'between' keeps values in [low, high]
      - 'outside' keeps values < low or > high
    - original_scale: if provided and contains transform info, will convert the
      threshold(s) back from log scale using inverse transform (currently log1p)

    Returns the boolean mask (True = keep/passes condition) and writes it to
    adata.obs[key_added].
    """
    if key not in adata.obs:
        raise KeyError(f"obs['{key}'] not found")
    vals = adata.obs[key].to_numpy().astype(float)

    # Inverse-transform thresholds if needed
    def inv(t: float) -> float:
        if original_scale and original_scale.get("transform") == "log1p":
            return float(np.expm1(t))
        return float(t)

    if direction in (">", ">=", "<", "<="):
        thr = inv(float(threshold))  # type: ignore
        if direction == ">":
            mask = vals > thr
        elif direction == ">=":
            mask = vals >= thr
        elif direction == "<":
            mask = vals < thr
        else:
            mask = vals <= thr
    elif direction in ("between", "outside"):
        if not (isinstance(threshold, (tuple, list)) and len(threshold) == 2):
            raise ValueError("threshold must be (low, high) for 'between'/'outside'")
        low, high = float(threshold[0]), float(threshold[1])
        low_i, high_i = inv(low), inv(high)
        if direction == "between":
            mask = (vals >= low_i) & (vals <= high_i)
        else:
            mask = (vals < low_i) | (vals > high_i)
    else:
        raise ValueError("direction must be one of {'>', '>=', '<', '<=', 'between', 'outside'}")

    if key_added is None:
        key_added = f"{key}_passes_{direction}"
    adata.obs[key_added] = mask.astype(bool)
    return mask


def batch_suggest_thresholds_for_obs(
    adata: AnnData,
    keys: Iterable[str],
    bins: int = 100,
    zero_inflated: bool = False,
    min_peak_prominence: float = 0.01,
    min_peak_distance: float = 0.05,
    kde_bandwidth: Optional[float] = None,
    transform: str = "none",
    show: bool = False,
    save_dir: Optional[str] = None,
) -> Dict[str, ThresholdSuggestion]:
    results: Dict[str, ThresholdSuggestion] = {}
    for k in keys:
        save_path = f"{save_dir.rstrip('/')}/{k}_thresholds.png" if save_dir is not None else None
        out = suggest_thresholds_for_obs(
            adata,
            key=k,
            bins=bins,
            zero_inflated=zero_inflated,
            min_peak_prominence=min_peak_prominence,
            min_peak_distance=min_peak_distance,
            kde_bandwidth=kde_bandwidth,
            transform=transform,
            show=show,
            save=save_path,
        )
        results[k] = ThresholdSuggestion(
            name=out["name"],
            peaks=list(out["peaks"]),
            valleys=list(out["valleys"]),
            thresholds=list(out["thresholds"]),
            histogram_bins=np.asarray(out["histogram_bins"]),
            histogram_counts=np.asarray(out["histogram_counts"]),
            kde_x=np.asarray(out["kde_x"]),
            kde_y=np.asarray(out["kde_y"]),
            transform=out.get("transform", "none"),
            original_scale=out.get("original_scale"),
        )
    return results


# -----------------------------------
# Gene expression API (kept for parity)
# -----------------------------------

def _get_gene_vector(
    adata: AnnData,
    gene: str | int,
    layer: Optional[str] = None,
    use_raw: bool = False,
) -> Tuple[str, np.ndarray]:
    """Return gene name and expression vector across cells.

    - gene can be a gene name (adata.var_names) or integer index.
    - layer: if provided, use adata.layers[layer]. Otherwise use adata.X.
    - use_raw: if True and adata.raw exists, use adata.raw.X/var_names.
    """
    if use_raw and getattr(adata, "raw", None) is not None:
        X = adata.raw.X
        var_names = np.array(adata.raw.var_names)
    else:
        X = adata.layers[layer] if layer is not None else adata.X
        var_names = np.array(adata.var_names)

    if isinstance(gene, int):
        idx = gene
        gene_name = str(var_names[idx])
    else:
        matches = np.where(var_names == gene)[0]
        if matches.size == 0:
            raise ValueError(f"Gene '{gene}' not found in var_names (use_raw={use_raw}).")
        idx = int(matches[0])
        gene_name = gene

    if sp.issparse(X):
        vec = X[:, idx].toarray().ravel()
    else:
        vec = np.asarray(X[:, idx]).ravel()
    return gene_name, vec


def suggest_thresholds_for_gene(
    adata: AnnData,
    gene: str | int,
    layer: Optional[str] = None,
    use_raw: bool = False,
    bins: int = 100,
    zero_inflated: bool = True,
    min_peak_prominence: float = 0.01,
    min_peak_distance: float = 0.05,
    kde_bandwidth: Optional[float] = None,
    show: bool = True,
    save: Optional[str] = None,
) -> Dict[str, Any]:
    gene_name, vec = _get_gene_vector(adata, gene=gene, layer=layer, use_raw=use_raw)
    out = suggest_thresholds_for_vector(
        values=vec,
        name=gene_name,
        bins=bins,
        zero_inflated=zero_inflated,
        min_peak_prominence=min_peak_prominence,
        min_peak_distance=min_peak_distance,
        kde_bandwidth=kde_bandwidth,
        transform="none",
        show=show,
        save=save,
    )
    return out


def apply_threshold_to_obs(
    adata: AnnData,
    gene: str | int,
    threshold: float,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    strictly_greater: bool = True,
) -> np.ndarray:
    gene_name, vec = _get_gene_vector(adata, gene=gene, layer=layer, use_raw=use_raw)
    mask = (vec > threshold) if strictly_greater else (vec >= threshold)
    if key_added is None:
        key_added = f"{gene_name}_positive"
    adata.obs[key_added] = mask.astype(bool)
    return mask


# ----------------
# CLI helper
# ----------------

def interactive_select_threshold(
    thresholds: List[float],
    prompt: str = "Select threshold index: ",
) -> Optional[float]:
    """Simple CLI selector. Returns chosen threshold value or None if not selected.

    Keep the human-in-the-loop aspect: the code proposes candidates, and a
    person chooses which to apply. Useful inside notebooks or scripts.
    """
    if not thresholds:
        print("No thresholds available.")
        return None
    print("Available threshold candidates (sorted):")
    for i, t in enumerate(sorted(thresholds)):
        print(f"  [{i}] {t:.6g}")
    try:
        idx = int(input(prompt))
    except Exception:
        print("Invalid input.")
        return None
    if 0 <= idx < len(thresholds):
        return sorted(thresholds)[idx]
    print("Index out of range.")
    return None
