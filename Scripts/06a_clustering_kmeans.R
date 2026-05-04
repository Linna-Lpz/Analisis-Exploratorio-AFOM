# =============================================================================
# CLUSTERING K-MEANS — LEYENDO MATRIZ RF DESDE CACHÉ .rds
# =============================================================================
library(clValid)
library(cluster)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# 1. LEER MATRIZ DESDE CACHÉ .rds
# =============================================================================
ruta_cache_matriz <- file.path(DIR_PROCESSED, "cache", "matriz_rf.rds")

if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché. Ejecuta primero el script de cálculo RF.")
}

cat("Cargando matriz RF desde caché...\n")
matriz_cuadrada <- readRDS(ruta_cache_matriz)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# 2. FUNCIÓN K-MEANS
# =============================================================================
metodo_KMEANS <- function(matriz_distancia_arboles, res_calidad_clusters) {
  
  posibles_k <- seq(2, 15)
  set.seed(2)
  
  for (i in posibles_k) {
    cat(sprintf("  Calculando k = %d...\n", i))
    
    clusters <- kmeans(matriz_distancia_arboles, centers = i, nstart = 25)
    
    puntaje_Dunn         <- dunn(distance = matriz_distancia_arboles, clusters$cluster)
    puntaje_Connectivity <- connectivity(distance = matriz_distancia_arboles, clusters$cluster)
    puntaje_Silhouette   <- mean(silhouette(clusters$cluster, matriz_distancia_arboles)[, 3])
    
    tmp <- c("Kmeans", i,
             round(puntaje_Dunn, 3),
             round(puntaje_Connectivity, 3),
             round(puntaje_Silhouette, 3))
    
    res_calidad_clusters <- rbind(res_calidad_clusters, tmp)
  }
  
  colnames(res_calidad_clusters) <- c("Method", "k", "Dunn", "Connectivity", "Silhouette")
  
  res_calidad_clusters$k            <- as.integer(res_calidad_clusters$k)
  res_calidad_clusters$Dunn         <- as.numeric(res_calidad_clusters$Dunn)
  res_calidad_clusters$Connectivity <- as.numeric(res_calidad_clusters$Connectivity)
  res_calidad_clusters$Silhouette   <- as.numeric(res_calidad_clusters$Silhouette)
  
  return(res_calidad_clusters)
}

# =============================================================================
# 3. EJECUTAR K-MEANS
# =============================================================================
cat("Ejecutando K-Means para k = 2 a 15...\n")
res_calidad_clusters <- data.frame()

tiempo_kmeans <- system.time({
  res_calidad_clusters <- metodo_KMEANS(
    matriz_distancia_arboles = matriz_cuadrada,
    res_calidad_clusters     = res_calidad_clusters
  )
})

print(res_calidad_clusters)

# =============================================================================
# 4. K ÓPTIMO Y ASIGNACIONES
# =============================================================================
k_optimo <- res_calidad_clusters$k[which.max(res_calidad_clusters$Silhouette)]
cat(sprintf("K óptimo según Silhouette: %d\n", k_optimo))

set.seed(2)
clusters_optimos <- kmeans(matriz_cuadrada, centers = k_optimo, nstart = 25)

asignaciones_df <- data.frame(
  Arbol   = rownames(matriz_cuadrada),
  Cluster = clusters_optimos$cluster
)

# Guardar asignaciones también como .rds para scripts posteriores
DIR_CACHE <- file.path(DIR_PROCESSED, "cache")
saveRDS(asignaciones_df,   file = file.path(DIR_CACHE, "kmeans_asignaciones.rds"))
saveRDS(clusters_optimos,  file = file.path(DIR_CACHE, "kmeans_modelo_optimo.rds"))

# =============================================================================
# 5. EXPORTAR CSVs
# =============================================================================
write.table(res_calidad_clusters,
            file = file.path(DIR_RESULTS, "kmeans_calidad_clusters.csv"),
            sep = ";", row.names = FALSE, quote = FALSE)

write.table(asignaciones_df,
            file = file.path(DIR_RESULTS, "kmeans_asignaciones_k_optimo.csv"),
            sep = ";", row.names = FALSE, quote = FALSE)

# =============================================================================
# 6. EXCEL CON TIEMPOS + RESULTADOS
# =============================================================================
t_kmeans <- as.numeric(tiempo_kmeans)

tiempos_df <- data.frame(
  Proceso          = c("K-Means (k=2 a 15)", "TOTAL"),
  Tiempo_Usuario_s = c(t_kmeans[1], t_kmeans[1]),
  Tiempo_Sistema_s = c(t_kmeans[2], t_kmeans[2]),
  Tiempo_Total_s   = c(t_kmeans[3], t_kmeans[3])
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Tiempos",
                              titulo_tabla = "Reporte de Tiempos de Ejecución (Segundos)",
                              datos        = tiempos_df,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Calidad_Clusters",
                              titulo_tabla = paste0("Índices de Calidad K-Means — K óptimo: ", k_optimo),
                              datos        = res_calidad_clusters,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Asignaciones_K_Optimo",
                              titulo_tabla = paste0("Asignaciones — K óptimo = ", k_optimo),
                              datos        = asignaciones_df,
                              anchos_col   = "auto")

ruta_excel <- file.path(DIR_RESULTS, "kmeans_resultados.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Resultados guardados en:", ruta_excel, "\n")
cat(sprintf("Resumen: k óptimo = %d | Silhouette = %.3f\n",
            k_optimo, max(res_calidad_clusters$Silhouette)))