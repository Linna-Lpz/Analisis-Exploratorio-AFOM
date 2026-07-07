# ==============================================================================
# instalar_paquetes.R
# Verifica e instala todas las dependencias del pipeline AFOM
# Ejecutar UNA SOLA VEZ al configurar el entorno (local o remoto)
# ==============================================================================

# 1. BiocManager (necesario para paquetes de Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Paquetes de CRAN
paquetes_cran <- c(
  # Infraestructura y rutas
  "here",
  # Árboles filogenéticos y distancias
  "ape", "TreeDist",
  # Clustering
  "cluster", "clValid", "mstknnclust", "igraph",
  # Reducción dimensional
  "uwot",       # UMAP
  "vegan",      # nMDS
  "Rtsne",      # t-SNE
  # Visualización
  "ggplot2", "patchwork", "pheatmap", "RColorBrewer", "viridisLite",
  # Manipulación de datos
  "dplyr", "tidyr", "stringr",
  # Lectura/escritura de Excel
  "openxlsx", "readxl",
  # APIs y web
  "httr", "jsonlite",
  # Enriquecimiento funcional (CRAN)
  "gprofiler2",
  # Utilidades
  "digest"
)

nuevos_cran <- paquetes_cran[!(paquetes_cran %in% installed.packages()[, "Package"])]
if (length(nuevos_cran) > 0) {
  message("Instalando desde CRAN: ", paste(nuevos_cran, collapse = ", "))
  install.packages(nuevos_cran)
} else {
  message("Todos los paquetes CRAN ya están instalados.")
}

# 3. Paquetes de Bioconductor
paquetes_bioc <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA"
)

nuevos_bioc <- paquetes_bioc[!(paquetes_bioc %in% installed.packages()[, "Package"])]
if (length(nuevos_bioc) > 0) {
  message("Instalando desde Bioconductor: ", paste(nuevos_bioc, collapse = ", "))
  BiocManager::install(nuevos_bioc, ask = FALSE)
} else {
  message("Todos los paquetes Bioconductor ya están instalados.")
}

# 4. Verificación final
todos <- c(paquetes_cran, paquetes_bioc)
faltantes <- todos[!(todos %in% installed.packages()[, "Package"])]

if (length(faltantes) == 0) {
  message("\n¡Todas las dependencias instaladas correctamente!")
} else {
  warning("\nFaltan por instalar: ", paste(faltantes, collapse = ", "))
}
