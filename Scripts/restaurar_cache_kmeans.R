# ==============================================================================
# ARCHIVO: restaurar_cache_kmeans.R
# PROPÓSITO: Reconstruir el archivo de caché de K-Means (kmeans_calidad_res.rds)
#            en Resultados/cache/ a partir de la salida impresa en consola.
#            Esto evita volver a ejecutar las 58 horas de cálculo en el servidor.
# ==============================================================================

source(here::here("Scripts", "config.R"))

cat("=== RESTAURANDO CACHÉ DE CALIDAD K-MEANS ===\n")

# Reconstruir dataframe de calidad a partir del log impreso
res_calidad_clusters <- data.frame(
  Method       = rep("Kmeans", 14),
  k            = 2:15,
  Dunn         = c(0.322, 0.230, 0.066, 0.060, 0.140, 0.140, 0.077, 0.033, 0.060, 0.126, 0.126, 0.077, 0.060, 0.060),
  Connectivity = c(11805.25, 22419.83, 28230.92, 33535.61, 36852.56, 37808.31, 39173.17, 40662.31, 41788.44, 42063.74, 42899.94, 43093.05, 43524.91, 43685.64),
  Silhouette   = c(0.165, 0.037, -0.018, -0.064, -0.094, -0.102, -0.115, -0.131, -0.144, -0.147, -0.159, -0.160, -0.167, -0.170),
  stringsAsFactors = FALSE
)

# Convertir tipos
res_calidad_clusters$k            <- as.integer(res_calidad_clusters$k)
res_calidad_clusters$Dunn         <- as.numeric(res_calidad_clusters$Dunn)
res_calidad_clusters$Connectivity <- as.numeric(res_calidad_clusters$Connectivity)
res_calidad_clusters$Silhouette   <- as.numeric(res_calidad_clusters$Silhouette)

# Tiempos de ejecución originales (211245 segundos)
tiempo_kmeans_valores <- c(211245.0, 0.0, 211245.0)
names(tiempo_kmeans_valores) <- c("user.self", "sys.self", "elapsed")
tiempo_kmeans <- structure(tiempo_kmeans_valores, class = "proc_time")

# Estructurar caché
cache_kmeans <- list(
  res_calidad_clusters = res_calidad_clusters,
  tiempo_kmeans        = tiempo_kmeans
)

# Asegurar que existe el directorio de caché
if (!dir.exists(DIR_CACHE)) {
  dir.create(DIR_CACHE, recursive = TRUE)
  cat("Creado directorio de caché:", DIR_CACHE, "\n")
}

ruta_cache <- file.path(DIR_CACHE, "kmeans_calidad_res.rds")
saveRDS(cache_kmeans, file = ruta_cache)

cat("Caché de calidad K-Means guardado exitosamente en:", ruta_cache, "\n")
cat("Distribución cargada:\n")
print(res_calidad_clusters)
cat("\n=== COMPLETADO ===\n")
