# =============================================================================
# 18_CALCULAR_BHI.R
# Cálculo del Índice de Homogeneidad Biológica (BHI)
# =============================================================================

source(here::here("Scripts", "config.R"))
library(openxlsx)
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
})

cat("=== CALCULANDO BHI PARA MST-kNN ===\n")

ruta_pca_excel <- file.path(DIR_RESULTS, paste0("pca_coordenadas_", NOMBRE_BDD, ".xlsx"))
if (!file.exists(ruta_pca_excel)) {
  stop("No se encontró el Excel PCA: ", ruta_pca_excel)
}

df_genes <- read.xlsx(ruta_pca_excel, sheet = "Genes_Agrupados", startRow = 2, colNames = TRUE)

algoritmo <- if(exists("ALGORITMO_DOWNSTREAM")) ALGORITMO_DOWNSTREAM else "AUTO"
COL_CLUSTER <- if (algoritmo == "MST-kNN") {
  "Cluster_MSTKNN_Iter"
} else if (algoritmo == "K-Means") {
  "Cluster_KMeans_Iter"
} else if (algoritmo == "CLARA") {
  "Cluster_CLARA_Iter"
} else if (algoritmo == "PAM") {
  "Cluster_PAM"
} else {
  "Cluster_MSTKNN_Iter"
}

if (!COL_CLUSTER %in% colnames(df_genes)) {
  cols_disponibles <- grep("^Cluster_", colnames(df_genes), value = TRUE)
  if (length(cols_disponibles) > 0) {
    COL_CLUSTER <- cols_disponibles[1]
  } else {
    stop("No hay columnas de cluster disponibles en Genes_Agrupados")
  }
}

df_genes$Cluster_Final <- df_genes[[COL_CLUSTER]]
df_validos <- df_genes[!is.na(df_genes$HGNC_Symbol) & !is.na(df_genes$Cluster_Final), ]

# Obtener mapeo de HGNC a GO
cat("Obteniendo anotaciones GO...\n")
# NOTA: dplyr y AnnotationDbi ambos exportan select(). Al cargar dplyr después
# de org.Hs.eg.db, dplyr::select() queda por encima en el search path y tapa
# a AnnotationDbi::select(), que es la que sabe operar sobre un objeto OrgDb.
# Se referencia el paquete explícitamente para evitar el conflicto sin
# depender del orden de los library().
hgnc_to_go <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = unique(df_validos$HGNC_Symbol),
    keytype = "SYMBOL",
    columns = "GO"
  )
)
hgnc_to_go <- hgnc_to_go[!is.na(hgnc_to_go$GO), ]

go_list <- split(hgnc_to_go$GO, hgnc_to_go$SYMBOL)

# Calcular BHI
cat("Calculando BHI por clúster...\n")
clusters <- unique(df_validos$Cluster_Final)
bhi_per_cluster <- c()

for (c in clusters) {
  genes <- df_validos$HGNC_Symbol[df_validos$Cluster_Final == c]
  genes <- genes[genes %in% names(go_list)]
  n <- length(genes)
  
  if (n < 2) next # BHI requiere pares
  
  matches <- 0
  total_pairs <- 0
  for (i in 1:(n-1)) {
    go_i <- go_list[[genes[i]]]
    for (j in (i+1):n) {
      go_j <- go_list[[genes[j]]]
      if (length(intersect(go_i, go_j)) > 0) {
        matches <- matches + 1
      }
      total_pairs <- total_pairs + 1
    }
  }
  bhi_per_cluster <- c(bhi_per_cluster, matches / total_pairs)
}

bhi_global <- mean(bhi_per_cluster)
cat(sprintf("-> BHI Global Obtenido: %.4f\n", bhi_global))

df_res <- data.frame(Metrica="BHI", Valor=bhi_global)
write.csv(df_res, file.path(DIR_RESULTS, "bhi_resultado.csv"), row.names=FALSE)
cat("Finalizado.\n")