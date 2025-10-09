import os
import glob
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import harmonypy as hm
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from scipy.stats import entropy, spearmanr
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix
from collections import Counter

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

'''
## Functions

'''

from statsmodels.stats.multitest import multipletests
def check_pc_correlation(adata, n_pcs=None, pc_col_name="X_pca"):
    """Check correlation between PCs and QC metrics to identify potential confounders."""
    if n_pcs is None:
        n_pcs_check = adata.obsm[pc_col_name].shape[1]
    else: n_pcs_check = min(n_pcs, adata.obsm[pc_col_name].shape[1])
    pcs = pd.DataFrame(
        adata.obsm[pc_col_name][:, :n_pcs_check],
        index=adata.obs_names,
        columns=["PC"+str(i+1).zfill(2) for i in range(n_pcs_check)]
    )

    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", 'pct_counts_ribo']
    qc = adata.obs[qc_cols].astype(float)
    qc['species'] = adata.obs['species'].astype('category').cat.codes
    qc_cols.append('species')

    rows = []
    for pc in pcs.columns:
        for cov in qc_cols:
            # Spearman es más robusto; cambia a pearsonr si prefieres lineal
            r, p = spearmanr(pcs[pc].values, qc[cov].values, nan_policy="omit")
            rows.append({"PC": pc, "covariate": cov, "r": r, "p": p})

    corr = pd.DataFrame(rows)
    corr["q"] = multipletests(corr["p"], method="fdr_bh")[1]


    return corr

def plot_pc_qc_correlation(corr, save_path=None):
    """Plot heatmap of PC vs QC metric correlations."""
    mat = corr.pivot(index="covariate", columns="PC", values="r")
    plt.figure(figsize=(10,3))
    im = plt.imshow(mat.values, aspect="auto", interpolation="nearest")
    plt.yticks(range(mat.shape[0]), mat.index)
    plt.xticks(range(mat.shape[1]), mat.columns, rotation=90)
    plt.colorbar(im, label="Spearman r")
    plt.title("Correlation PC vs QC")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()
    

'''
## Pre-processing
- Start with raw_merged data
- Filter
- Normalize
- Log transform
- Find HVG
- Combat (batch correction)
'''

#Indicate file paths
data_path = '/home/raquelcr/scanpy/'
fig_path = '/home/raquelcr/scanpy/figures'
out_path = '/home/raquelcr/scanpy/predicted'
os.makedirs(fig_path, exist_ok=True)
os.makedirs(out_path, exist_ok=True)

rev_n = 7
filename = os.path.join(data_path, f'all_adata_outer{rev_n}.h5ad')
adata = sc.read(filename)
print(f"Succesfully loaded {filename}. Shape: {adata.shape}")

print(adata.obs['dataset_name'].value_counts())
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", 'pct_counts_ribo'],
    jitter=0.4,
    multi_panel=True,
    save= os.path.join(fig_path,f'qc_violin{rev_n}.png')
)
sc.pp.scale(adata)
sc.pp.pca(adata, svd_solver='arpack', use_highly_variable=False)
sc.pl.pca(
    adata,
    color=["cell_supertype","cell_supertype", "dataset_name", "dataset_name"],
    dimensions=[(0,1),(2,3), (0,1),(2,3)],
    ncols=2,
    size=2,
    save=f'_by_supertype{rev_n}.png',
    legend_loc='none' 
)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f'{rev_n}.png')