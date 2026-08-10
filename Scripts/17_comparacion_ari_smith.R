# =============================================================================
# 17_comparacion_ari_smith.R
# =============================================================================
library(here)
library(ape)
library(TreeDist)
library(cluster)
library(openxlsx)
library(mclust)
library(mstknnclust)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

cat("=== COMPARACION RF VS METRICAS DE SMITH (2020) ===\n")

# =============================================================================
# 1. CARGAR LOS ARBOLES DEL CORE SET DESDE CACHE
#    Generado en 02_filtro_arboles_medioide.R
# =============================================================================
ruta_cache_core_set <- file.path(DIR_CACHE, "conjunto_core_set.rds")
if (!file.exists(ruta_cache_core_set))
  stop("Falta: ", ruta_cache_core_set,
       "\nEjecuta primero 02_filtro_arboles_medioide.R")

bosque_324 <- readRDS(ruta_cache_core_set)
class(bosque_324) <- "multiPhylo"

cat(sprintf("Arboles cargados: %d\n", length(bosque_324)))

# =============================================================================
# 2. CALCULAR DISTANCIA DE SMITH (ClusteringInfoDistance)
# =============================================================================
cat("Calculando ClusteringInfoDistance...\n")

matriz_smith  <- ClusteringInfoDistance(bosque_324, normalize = TRUE)
mat_smith_324 <- as.matrix(matriz_smith)

# Forzamos SIEMPRE los nombres reales de los arboles como rownames/colnames,
# sin importar lo que haya devuelto ClusteringInfoDistance/as.matrix (que
# puede traer indices numericos "1","2",... en vez de los nombres reales).
# Esto asume que el orden de bosque_324 se preserva al calcular la distancia,
# igual que se asume en el resto del pipeline (p. ej. con RobinsonFoulds).
stopifnot(nrow(mat_smith_324) == length(bosque_324))
rownames(mat_smith_324) <- names(bosque_324)
colnames(mat_smith_324) <- names(bosque_324)

# =============================================================================
# 3. CLUSTERING CON MST-KNN SOBRE LA MATRIZ DE SMITH
# =============================================================================
cat("Clustering MST-KNN sobre matriz de Smith...\n")

set.seed(42)
mst_res_smith <- mstknnclust::mst.knn(mat_smith_324)
k_optimo_smith <- mst_res_smith$cnumber

cat(sprintf("K optimo (Smith, mst-knn): %d\n", k_optimo_smith))

asig_smith_df <- data.frame(
  Arbol         = names(mst_res_smith$cluster),
  Cluster_Smith = as.integer(mst_res_smith$cluster),
  stringsAsFactors = FALSE
)

# =============================================================================
# 4. COMPARAR CONTRA EL CLUSTERING RF (control_324_clustering_....xlsx)
# =============================================================================
ruta_cl_control <- file.path(DIR_RESULTS, paste0("control_324_clustering_", NOMBRE_BDD, ".xlsx"))
asig_control <- read.xlsx(ruta_cl_control, sheet = "Asignaciones", startRow = 2)

# Cruce correcto: por nombre real de arbol (Arbol), no por indice numerico
comparacion <- merge(asig_control, asig_smith_df, by = "Arbol")

n_perdidos <- nrow(asig_control) - nrow(comparacion)
if (n_perdidos != 0) {
  warning(sprintf(
    "%d arboles de asig_control no encontraron match en asig_smith_df. ",
    n_perdidos
  ))
}

cat("Tabla de comparacion:\n")
print(table(comparacion$Cluster, comparacion$Cluster_Smith))

ari_score <- adjustedRandIndex(comparacion$Cluster, comparacion$Cluster_Smith)
cat("-> ARI entre RF (clusters=", length(unique(comparacion$Cluster)),
    ") y Smith/mst-knn (clusters=", k_optimo_smith, "):",
    round(ari_score, 4), "\n", sep = "")

# =============================================================================
# 5. EXPORTAR RESULTADOS
# =============================================================================
resumen_ari <- data.frame(
  Metrica = c("ARI", "K_RF", "K_Smith_mstknn", "N_arboles_comparados"),
  Valor   = c(ari_score,
              length(unique(comparacion$Cluster)),
              k_optimo_smith,
              nrow(comparacion)),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "ARI")
writeData(wb, "ARI", resumen_ari)
addWorksheet(wb, "Asignaciones")
writeData(wb, "Asignaciones", comparacion)
saveWorkbook(wb, file.path(DIR_RESULTS, paste0("control_324_smith_ari_", NOMBRE_BDD, ".xlsx")), overwrite = TRUE)

cat("=== FIN ===\n")