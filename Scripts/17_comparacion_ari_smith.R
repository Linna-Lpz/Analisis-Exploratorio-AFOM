# =============================================================================
# 17_comparacion_ari_smith.R
# =============================================================================

library(here)
library(ape)
library(TreeDist)
library(cluster)
library(openxlsx)
library(mclust)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

adjusted_rand_index <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c) / (a + b + c + d)) /
    ((a + b + a + c) / 2 - (a + b) * (a + c) / (a + b + c + d))
  return(ARI)
}

cat("=== COMPARACION RF VS METRICAS DE SMITH (2020) ===\n")

ruta_medioide <- file.path(DIR_RESULTS, paste0("ranking_medioide_", NOMBRE_BDD, ".xlsx"))
ranking_df  <- read.xlsx(ruta_medioide, sheet = 1, startRow = 2)
ids_324     <- as.character(ranking_df$nombre_arbol)

archivos_bosque <- list.files(DIR_CACHE, pattern = "^bosque_.*\\.rds$", full.names = TRUE)
bosque_rds <- archivos_bosque[which.min(file.size(archivos_bosque))]
bosque_total <- readRDS(bosque_rds)
bosque_324 <- bosque_total[names(bosque_total) %in% ids_324]
bosque_324 <- bosque_324[!sapply(bosque_324, is.null)]
class(bosque_324) <- "multiPhylo"
rm(bosque_total); gc()

cat("Calculando ClusteringInfoDistance...\n")
matriz_smith <- ClusteringInfoDistance(bosque_324, normalize = TRUE)
mat_smith_324 <- as.matrix(matriz_smith)

ruta_cl_control <- file.path(DIR_RESULTS, paste0("control_324_clustering_", NOMBRE_BDD, ".xlsx"))
asig_control <- read.xlsx(ruta_cl_control, sheet = "Asignaciones", startRow = 2)

k_optimo <- length(unique(asig_control$Cluster))

set.seed(42)
res_clara_smith <- clara(mat_smith_324, k = k_optimo, metric = "euclidean", 
                         samples = 50, sampsize = min(nrow(mat_smith_324), 40 + 2 * k_optimo), 
                         keep.data = FALSE, rngR = TRUE)

asig_smith_df <- data.frame(
  Arbol = as.character(1:nrow(mat_smith_324)),
  Cluster_Smith = res_clara_smith$clustering,
  stringsAsFactors = FALSE
)

comparacion <- merge(asig_control, asig_smith_df, by="Arbol")
cat("Tabla de comparacion:\n")
print(table(comparacion$Cluster, comparacion$Cluster_Smith))
ari_score <- adjustedRandIndex(comparacion$Cluster, comparacion$Cluster_Smith)
cat("-> ARI entre RF y Smith (ClusteringInfoDistance):", round(ari_score, 4), "\n")

resumen_ari <- data.frame(Metrica = "ARI", Valor = ari_score)
wb <- createWorkbook()
addWorksheet(wb, "ARI")
writeData(wb, "ARI", resumen_ari)
addWorksheet(wb, "Asignaciones")
writeData(wb, "Asignaciones", comparacion)
saveWorkbook(wb, file.path(DIR_RESULTS, paste0("control_324_smith_ari_", NOMBRE_BDD, ".xlsx")), overwrite = TRUE)

cat("=== FIN ===\n")
