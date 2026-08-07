# =============================================================================
# CLUSTERING K-MEANS — LEYENDO MATRIZ RF DESDE CACHÉ .rds
# =============================================================================
library(clValid)
library(cluster)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# LEER MATRIZ DESDE CACHÉ .rds
# =============================================================================
ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché. Ejecuta primero el script de cálculo RF.")
}

cat("Cargando matriz RF desde caché...\n")
matriz_cuadrada <- readRDS(ruta_cache_matriz)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# EMBEDDING MDS (k=20) PARA ESPACIO MÉTRICO CORRECTO
# =============================================================================
ruta_cache_mds <- file.path(DIR_CACHE, "mds_coords_k20.rds")
N_DIMS <- 20

if (file.exists(ruta_cache_mds)) {
  cat(sprintf("Cargando embedding MDS (k=%d) desde caché...\n", N_DIMS))
  mds_coords <- readRDS(ruta_cache_mds)
} else {
  cat(sprintf("Calculando embedding MDS (k=%d) para K-Means...\n", N_DIMS))
  tiempo_mds <- system.time({
    mds_coords <- cmdscale(as.dist(matriz_cuadrada), k = N_DIMS)
  })
  saveRDS(mds_coords, file = ruta_cache_mds)
  cat(sprintf("MDS calculado en %.1f segundos y guardado en caché.\n", tiempo_mds[3]))
}

# Calcular varianza explicada (aproximación rápida)
var_total <- sum(matriz_cuadrada^2) / (2 * nrow(matriz_cuadrada)^2)
var_mds   <- sum(apply(mds_coords, 2, var))
var_exp   <- (var_mds / var_total) * 100
cat(sprintf("Varianza explicada por %d dimensiones MDS: %.2f%%\n\n", N_DIMS, var_exp))

# =============================================================================
# FUNCIÓN K-MEANS
# =============================================================================
metodo_KMEANS <- function(matriz_distancia_arboles, mds_coords, res_calidad_clusters) {
  
  posibles_k <- seq(2, 15)
  set.seed(2)
  
  for (i in posibles_k) {
    cat(sprintf("  Calculando k = %d...\n", i))
    
    # K-Means opera sobre las coordenadas MDS (espacio euclidiano representativo de RF)
    clusters <- kmeans(mds_coords, centers = i, nstart = 25)
    
    # Las métricas de validación se evalúan sobre la matriz de disimilitud RF original
    puntaje_Dunn         <- dunn(distance = matriz_distancia_arboles, clusters$cluster)
    puntaje_Connectivity <- connectivity(distance = matriz_distancia_arboles, clusters$cluster)
    puntaje_Silhouette   <- mean(silhouette(clusters$cluster, as.dist(matriz_distancia_arboles))[, 3])
    
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
# EJECUTAR K-MEANS (con caché condicional)
# =============================================================================
ruta_cache_calidad <- file.path(DIR_CACHE, "kmeans_calidad_res.rds")

if (file.exists(ruta_cache_calidad)) {
  cat("Resultados de calidad de clusters K-Means encontrados en caché. Cargando...\n")
  cache_kmeans <- readRDS(ruta_cache_calidad)
  res_calidad_clusters <- cache_kmeans$res_calidad_clusters
  tiempo_kmeans        <- cache_kmeans$tiempo_kmeans
} else {
  cat("Ejecutando K-Means para k = 2 a 15...\n")
  res_calidad_clusters <- data.frame()
  
  tiempo_kmeans <- system.time({
    res_calidad_clusters <- metodo_KMEANS(
      matriz_distancia_arboles = matriz_cuadrada,
      mds_coords               = mds_coords,
      res_calidad_clusters     = res_calidad_clusters
    )
  })
  
  # Guardar en caché para evitar tener que repetir si ocurre algún error posterior
  cache_kmeans <- list(
    res_calidad_clusters = res_calidad_clusters,
    tiempo_kmeans        = tiempo_kmeans
  )
  saveRDS(cache_kmeans, file = ruta_cache_calidad)
  cat("Resultados de K-Means guardados en caché:", ruta_cache_calidad, "\n")
}

print(res_calidad_clusters)


# =============================================================================
# K ÓPTIMO Y ASIGNACIONES
# =============================================================================
k_optimo <- res_calidad_clusters$k[which.max(res_calidad_clusters$Silhouette)]
cat(sprintf("K óptimo según Silhouette (sobre dist RF): %d\n", k_optimo))

set.seed(2)
clusters_optimos <- kmeans(mds_coords, centers = k_optimo, nstart = 25)

asignaciones_df <- data.frame(
  Arbol   = rownames(matriz_cuadrada),
  Cluster = clusters_optimos$cluster
)

# Guardar asignaciones también como .rds para scripts posteriores
saveRDS(asignaciones_df,   file = file.path(DIR_CACHE, "kmeans_asignaciones.rds"))
saveRDS(clusters_optimos,  file = file.path(DIR_CACHE, "kmeans_modelo_optimo.rds"))

# =============================================================================
# ASIGNACIONES ADICIONALES PARA k = 10 Y k = 15
# =============================================================================
k_extra <- c(10, 15)
asignaciones_extra <- list()

for (ke in k_extra) {
  cat(sprintf("Generando asignaciones para k = %d...\n", ke))
  set.seed(2)
  cl_extra <- kmeans(mds_coords, centers = ke, nstart = 25)

  asig_extra <- data.frame(
    Arbol   = rownames(matriz_cuadrada),
    Cluster = cl_extra$cluster
  )

  # Agregar columna con tamaño de cada cluster para facilitar revisión
  tamanos <- as.integer(table(cl_extra$cluster))
  asig_extra$Tamano_Cluster <- tamanos[asig_extra$Cluster]

  asignaciones_extra[[as.character(ke)]] <- asig_extra

  # Resumen de tamaños por cluster
  resumen_k <- data.frame(
    Cluster = seq_len(ke),
    Tamano  = tamanos
  )
  cat(sprintf("  k=%d — Tamaños de clusters: min=%d, max=%d\n",
              ke, min(tamanos), max(tamanos)))
}

# =============================================================================
# EXCEL CON TIEMPOS + RESULTADOS
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
                              anchos_col   = "auto")# =============================================================================
# CLUSTERING K-MEANS — LEYENDO MATRIZ RF DESDE CACHÉ .rds
# =============================================================================
library(clValid)
library(cluster)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# LEER MATRIZ DESDE CACHÉ .rds
# =============================================================================
ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché. Ejecuta primero el script de cálculo RF.")
}

cat("Cargando matriz RF desde caché...\n")
matriz_cuadrada <- readRDS(ruta_cache_matriz)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# FUNCIÓN K-MEANS
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
# EJECUTAR K-MEANS (con caché condicional)
# =============================================================================
ruta_cache_calidad <- file.path(DIR_CACHE, "kmeans_calidad_res.rds")

if (file.exists(ruta_cache_calidad)) {
  cat("Resultados de calidad de clusters K-Means encontrados en caché. Cargando...\n")
  cache_kmeans <- readRDS(ruta_cache_calidad)
  res_calidad_clusters <- cache_kmeans$res_calidad_clusters
  tiempo_kmeans        <- cache_kmeans$tiempo_kmeans
} else {
  cat("Ejecutando K-Means para k = 2 a 15...\n")
  res_calidad_clusters <- data.frame()
  
  tiempo_kmeans <- system.time({
    res_calidad_clusters <- metodo_KMEANS(
      matriz_distancia_arboles = matriz_cuadrada,
      res_calidad_clusters     = res_calidad_clusters
    )
  })
  
  # Guardar en caché para evitar tener que repetir si ocurre algún error posterior
  cache_kmeans <- list(
    res_calidad_clusters = res_calidad_clusters,
    tiempo_kmeans        = tiempo_kmeans
  )
  saveRDS(cache_kmeans, file = ruta_cache_calidad)
  cat("Resultados de K-Means guardados en caché:", ruta_cache_calidad, "\n")
}

print(res_calidad_clusters)


# =============================================================================
# K ÓPTIMO Y ASIGNACIONES
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
saveRDS(asignaciones_df,   file = file.path(DIR_CACHE, "kmeans_asignaciones.rds"))
saveRDS(clusters_optimos,  file = file.path(DIR_CACHE, "kmeans_modelo_optimo.rds"))

# =============================================================================
# ASIGNACIONES ADICIONALES PARA k = 10 Y k = 15
# =============================================================================
k_extra <- c(10, 15)
asignaciones_extra <- list()

for (ke in k_extra) {
  cat(sprintf("Generando asignaciones para k = %d...\n", ke))
  set.seed(2)
  cl_extra <- kmeans(matriz_cuadrada, centers = ke, nstart = 25)
  
  asig_extra <- data.frame(
    Arbol   = rownames(matriz_cuadrada),
    Cluster = cl_extra$cluster
  )
  
  # Agregar columna con tamaño de cada cluster para facilitar revisión
  tamanos <- as.integer(table(cl_extra$cluster))
  asig_extra$Tamano_Cluster <- tamanos[asig_extra$Cluster]
  
  asignaciones_extra[[as.character(ke)]] <- asig_extra
  
  # Resumen de tamaños por cluster
  resumen_k <- data.frame(
    Cluster = seq_len(ke),
    Tamano  = tamanos
  )
  cat(sprintf("  k=%d — Tamaños de clusters: min=%d, max=%d\n",
              ke, min(tamanos), max(tamanos)))
}

# =============================================================================
# EXCEL CON TIEMPOS + RESULTADOS
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

# Hojas adicionales para k = 10 y k = 15
for (ke in k_extra) {
  nombre_hoja <- paste0("Asignaciones_K_", ke)
  titulo      <- paste0("Asignaciones — K = ", ke)
  wb <- agregar_hoja_formateada(wb           = wb,
                                nombre_hoja  = nombre_hoja,
                                titulo_tabla = titulo,
                                datos        = asignaciones_extra[[as.character(ke)]],
                                anchos_col   = "auto")
}

ruta_excel <- file.path(DIR_RESULTS, "kmeans_resultados.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Resultados guardados en:", ruta_excel, "\n")
cat(sprintf("Resumen: k óptimo = %d | Silhouette = %.3f\n",
            k_optimo, max(res_calidad_clusters$Silhouette)))

# Hojas adicionales para k = 10 y k = 15
for (ke in k_extra) {
  nombre_hoja <- paste0("Asignaciones_K_", ke)
  titulo      <- paste0("Asignaciones — K = ", ke)
  wb <- agregar_hoja_formateada(wb           = wb,
                                nombre_hoja  = nombre_hoja,
                                titulo_tabla = titulo,
                                datos        = asignaciones_extra[[as.character(ke)]],
                                anchos_col   = "auto")
}

ruta_excel <- file.path(DIR_RESULTS, "kmeans_resultados.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Resultados guardados en:", ruta_excel, "\n")
cat(sprintf("Resumen: k óptimo = %d | Silhouette = %.3f\n",
            k_optimo, max(res_calidad_clusters$Silhouette)))