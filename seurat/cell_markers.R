library(tibble)
library(dplyr)
library(tidyr)

investigado <- tibble(
  tipo_celular = c(
    rep("Neuronas", 12),
    rep("Glía", 8),
    rep("Progenitoras", 3),
    rep("Vasculares", 3),
    rep("Ependimarias", 2)
  ),
  subtipo_celular = c(
    rep("Dopaminérgicas (SNc)", 6),
    rep("Dopaminérgicas (VTA)", 3),
    rep("GABAérgicas", 3),
    rep("Astrocitos", 6),
    rep("Oligodendrocitos", 3),
    rep("Microglía", 1),
    rep("Progenitoras neurales", 3),
    rep("Endoteliales", 3),
    rep("Ependimarias", 2)
  ),
  marcador = c(
    # Neuronas SNc
    "TH", "SLC6A3", "SLC18A2", "FOXA2", "CALB1", "SOX6",
    # Neuronas VTA
    "TH", "SLC6A3", "NR4A2",
    # GABAérgicas
    "GAD1", "GAD2", "SLC32A1",
    # Astrocitos
    "GFAP", "AQP4", "SLC1A2", "ALDH1L1", 'GJB6', 'GJA1'
    # Oligodendrocitos
    "OLIG1", "OLIG2", "SOX10",
    # Microglía
    "CX3CR1",
    # Progenitoras
    "NES", "SOX2", "ASCL1",
    # Vasculares
    "CLDN5", "PECAM1", "FLT1",
    # Ependimarias
    "FOXJ1", "SPDEF"
  )
  neg_marker = c(
) %>% arrange(tipo_celular, subtipo_celular)

# 1. Mesencéfalo (expandido)
mesencefalo <- tibble(
  tipo_celular = c(
    rep("Neuronas", 12),
    rep("Glía", 8),
    rep("Progenitoras", 3),
    rep("Vasculares", 3),
    rep("Ependimarias", 2)
  ),
  subtipo_celular = c(
    rep("Dopaminérgicas (SNc)", 6),
    rep("Dopaminérgicas (VTA)", 3),
    rep("GABAérgicas", 3),
    rep("Astrocitos", 4),
    rep("Oligodendrocitos", 3),
    rep("Microglía", 1),
    rep("Progenitoras neurales", 3),
    rep("Endoteliales", 3),
    rep("Ependimarias", 2)
  ),
  marcador = c(
    # Neuronas SNc
    "TH", "SLC6A3", "SLC18A2", "FOXA2", "CALB1", "SOX6",
    # Neuronas VTA
    "TH", "SLC6A3", "NR4A2",
    # GABAérgicas
    "GAD1", "GAD2", "SLC32A1",
    # Astrocitos
    "GFAP", "AQP4", "SLC1A2", "ALDH1L1",
    # Oligodendrocitos
    "OLIG1", "OLIG2", "SOX10",
    # Microglía
    "CX3CR1",
    # Progenitoras
    "NES", "SOX2", "ASCL1",
    # Vasculares
    "CLDN5", "PECAM1", "FLT1",
    # Ependimarias
    "FOXJ1", "SPDEF"
  )
) %>% arrange(tipo_celular, subtipo_celular)

# 2. Hipocampo (expandido)
hipocampo <- tibble(
  tipo_celular = c(
    rep("Neuronas", 9),
    rep("Interneuronas", 9),
    rep("Glía", 8),
    rep("Progenitoras", 3),
    rep("Vasculares", 3)
  ),
  subtipo_celular = c(
    rep("Piramidales CA1", 3),
    rep("Piramidales CA3", 2),
    rep("Granulares DG", 4),
    rep("PVALB+", 3),
    rep("SST+", 3),
    rep("VIP+", 3),
    rep("Astrocitos", 4),
    rep("Oligodendrocitos", 3),
    rep("Microglía", 1),
    rep("Progenitoras", 3),
    rep("Endoteliales", 3)
  ),
  marcador = c(
    # CA1
    "SLC17A7", "CAMK2A", "SATB2",
    # CA3
    "SLC17A7", "PCP4",
    # DG
    "PROX1", "GRIA1", "CALB2", "TBR1",
    # PVALB+
    "PVALB", "GAD1", "SLC32A1",
    # SST+
    "SST", "NPY", "LHX6",
    # VIP+
    "VIP", "CCK", "CALB2",
    # Astrocitos
    "GFAP", "AQP4", "SLC1A2", "FABP7",
    # Oligodendrocitos
    "OLIG1", "PLP1", "SOX10",
    # Microglía
    "TMEM119",
    # Progenitoras
    "SOX2", "NES", "DCX",
    # Vasculares
    "CLDN5", "PECAM1", "PDGFRB"
  )
) %>% arrange(tipo_celular, subtipo_celular)

# 3. Corteza Cerebral (expandido)
corteza <- tibble(
  tipo_celular = c(
    rep("Neuronas", 12),
    rep("Interneuronas", 12),
    rep("Glía", 9),
    rep("Progenitoras", 3),
    rep("Vasculares", 3)
  ),
  subtipo_celular = c(
    rep("Piramidales L2/3", 3),
    rep("Piramidales L4", 2),
    rep("Piramidales L5", 4),
    rep("Piramidales L6", 3),
    rep("PVALB+", 3),
    rep("SST+", 3),
    rep("VIP+", 3),
    rep("LAMP5+", 3),
    rep("Astrocitos", 4),
    rep("Oligodendrocitos", 3),
    rep("Microglía", 2),
    rep("Progenitoras", 3),
    rep("Endoteliales", 3)
  ),
  marcador = c(
    # L2/3
    "SLC17A7", "CUX1", "CUX2",
    # L4
    "RORB", "SYT6",
    # L5
    "SLC17A7", "FEZF2", "BCL11B", "TLE4",
    # L6
    "SLC17A7", "FOXP2", "NRP1",
    # PVALB+
    "PVALB", "GAD1", "SLC32A1",
    # SST+
    "SST", "NPY", "LHX6",
    # VIP+
    "VIP", "CCK", "CALB2",
    # LAMP5+
    "LAMP5", "SNCG", "RELN",
    # Astrocitos
    "GFAP", "AQP4", "SLC1A2", "APOE",
    # Oligodendrocitos
    "OLIG1", "PLP1", "MBP",
    # Microglía
    "CX3CR1", "P2RY12",
    # Progenitoras
    "SOX2", "NES", "ASCL1",
    # Vasculares
    "CLDN5", "PECAM1", "FLT1"
  ))
) %>% arrange(tipo_celular, subtipo_celular)

mesencefalo %>% 
  filter(subtipo_celular == "Dopaminérgicas (SNc)") %>%
  pull(marcador)  # Extraer los marcadores

write.csv(corteza,paste0(output_dir,'/cell_markers_cerebral_cortex.csv'),row.names = FALSE)
write.csv(hipocampo,paste0(output_dir,'/cell_markers_hippocampus.csv'),row.names = FALSE)
write.csv(mesencefalo,paste0(output_dir,'/cell_markers_midbrain.csv'),row.names = FALSE)

