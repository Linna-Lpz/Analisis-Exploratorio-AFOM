# =============================================================================
# CLUSTERING K-MEANS DESDE COORDENADAS (UMAP / t-SNE / PCA)
# =============================================================================
library(clValid)
library(cluster)
library(openxlsx)
library(readxl)
library(ggplot2)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# 1. PARÁMETROS
# =============================================================================
METODO_REDUCCION <- "pca"    # "umap", "tsne" o "pca"

ruta_coords <- file.path(DIR_RESULTS,
                         sprintf("%s_coordenadas_%s.xlsx", METODO_REDUCCION, NOMBRE_BDD))

# =============================================================================
# 2. LECTURA DE COORDENADAS DESDE EXCEL
# =============================================================================
cat("Cargando coordenadas desde Excel...\n")

if (!file.exists(ruta_coords)) {
  stop("No se encontró el archivo: ", ruta_coords)
}

coords_df       <- read_excel(ruta_coords, skip = 1)
nombres_arboles <- as.character(coords_df[[1]])
cols_coords     <- names(coords_df)[-1]

# Limpiar comas decimales (Excel en español)
coords_df[cols_coords] <- lapply(coords_df[cols_coords], function(x) {
  if (is.character(x)) as.numeric(gsub(",", ".", x)) else as.numeric(x)
})

coords <- as.matrix(coords_df[, -1])
rownames(coords) <- nombres_arboles

# Para PCA: usar solo las 2 primeras componentes
if (METODO_REDUCCION == "pca") {
  coords <- coords[, 1:2]
  cat("PCA: usando solo PC1 y PC2\n")
}

cat("Coordenadas cargadas:", nrow(coords), "árboles ×", ncol(coords), "dimensiones\n")
cat("NAs en coords:", sum(is.na(coords)), "\n")

# =============================================================================
# 3. FUNCIÓN K-MEANS PARA COORDENADAS
# =============================================================================
metodo_KMEANS_coords <- function(coords, res_calidad_clusters) {
  
  matriz_dist <- dist(coords, method = "euclidean")
  
  posibles_k <- seq(2, 15)
  set.seed(2)
  
  for (i in posibles_k) {
    cat(sprintf("  Calculando k = %d...\n", i))
    
    clusters <- kmeans(coords, centers = i, nstart = 25)
    
    puntaje_Dunn         <- dunn(distance = matriz_dist, clusters$cluster)
    puntaje_Connectivity <- connectivity(distance = matriz_dist, clusters$cluster)
    puntaje_Silhouette   <- mean(silhouette(clusters$cluster, matriz_dist)[, 3])
    
    tmp <- c(METODO_REDUCCION, i,
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
# 4. EJECUTAR K-MEANS
# =============================================================================
cat(sprintf("Ejecutando K-Means sobre coordenadas %s para k = 2 a 15...\n", METODO_REDUCCION))
res_calidad_clusters <- data.frame()

tiempo_kmeans <- system.time({
  res_calidad_clusters <- metodo_KMEANS_coords(
    coords               = coords,
    res_calidad_clusters = res_calidad_clusters
  )
})

# =============================================================================
# 5. K ÓPTIMO
# =============================================================================
k_optimo <- res_calidad_clusters$k[which.max(res_calidad_clusters$Silhouette)]
cat(sprintf("K óptimo según Silhouette: %d\n", k_optimo))

# =============================================================================
# 6. GRÁFICO DE SILHOUETTE POR K
# =============================================================================
cat("\nTabla completa de métricas:\n")
print(res_calidad_clusters)

p_silhouette <- ggplot(res_calidad_clusters, aes(x = k, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = k_optimo, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 2:15) +
  labs(
    title    = "Silhouette por número de clusters",
    subtitle = "Línea roja = k óptimo seleccionado automáticamente"
  ) +
  theme_bw()

print(p_silhouette)

# =============================================================================
# 7. ASIGNACIONES CON K ÓPTIMO
# =============================================================================
set.seed(2)
clusters_optimos <- kmeans(coords, centers = k_optimo, nstart = 25)

asignaciones_df <- data.frame(
  Arbol   = rownames(coords),
  Cluster = clusters_optimos$cluster
)

# =============================================================================
# 8. GUARDAR CACHÉ .rds
# =============================================================================
DIR_CACHE <- file.path(DIR_PROCESSED, "cache")

saveRDS(asignaciones_df,
        file = file.path(DIR_CACHE, sprintf("kmeans_%s_asignaciones.rds", METODO_REDUCCION)))
saveRDS(clusters_optimos,
        file = file.path(DIR_CACHE, sprintf("kmeans_%s_modelo_optimo.rds", METODO_REDUCCION)))

# =============================================================================
# 9. EXPORTAR CSVs
# =============================================================================
write.table(res_calidad_clusters,
            file = file.path(DIR_RESULTS, sprintf("kmeans_%s_calidad_clusters.csv", METODO_REDUCCION)),
            sep = ";", row.names = FALSE, quote = FALSE)

write.table(asignaciones_df,
            file = file.path(DIR_RESULTS, sprintf("kmeans_%s_asignaciones_k_optimo.csv", METODO_REDUCCION)),
            sep = ";", row.names = FALSE, quote = FALSE)

# =============================================================================
# 10. EXCEL CON TIEMPOS + RESULTADOS
# =============================================================================
t_kmeans <- as.numeric(tiempo_kmeans)

tiempos_df <- data.frame(
  Proceso          = c(sprintf("K-Means %s (k=2 a 15)", METODO_REDUCCION), "TOTAL"),
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
                              titulo_tabla = sprintf("Índices de Calidad K-Means %s — K óptimo: %d",
                                                     METODO_REDUCCION, k_optimo),
                              datos        = res_calidad_clusters,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Asignaciones_K_Optimo",
                              titulo_tabla = sprintf("Asignaciones %s — K óptimo = %d",
                                                     METODO_REDUCCION, k_optimo),
                              datos        = asignaciones_df,
                              anchos_col   = "auto")

ruta_excel <- file.path(DIR_RESULTS, sprintf("kmeans_%s_resultados.xlsx", METODO_REDUCCION))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)

cat("Resultados guardados en:", ruta_excel, "\n")
cat(sprintf("Resumen: método = %s | k óptimo = %d | Silhouette = %.3f\n",
            METODO_REDUCCION, k_optimo, max(res_calidad_clusters$Silhouette)))

# =============================================================================
# 11. GRÁFICO DE CLUSTERS — exportar como .png
# =============================================================================
col_x <- colnames(coords)[1]
col_y <- colnames(coords)[2]

plot_df <- data.frame(
  coords,
  Cluster = as.factor(clusters_optimos$cluster),
  Arbol   = rownames(coords)
)

p_clusters <- ggplot(plot_df, aes(x = .data[[col_x]],
                                  y = .data[[col_y]],
                                  color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(group = Cluster), type = "norm",
               linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = sprintf("K-Means sobre %s — K óptimo = %d", METODO_REDUCCION, k_optimo),
    subtitle = sprintf("Silhouette = %.3f", max(res_calidad_clusters$Silhouette)),
    x        = col_x,
    y        = col_y,
    color    = "Cluster"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right"
  )

ruta_png <- file.path(DIR_RESULTS, sprintf("kmeans_%s_clusters.png", METODO_REDUCCION))

ggsave(ruta_png,
       plot   = p_clusters,
       width  = 8,
       height = 6,
       dpi    = 300)

cat("Gráfico guardado en:", ruta_png, "\n")