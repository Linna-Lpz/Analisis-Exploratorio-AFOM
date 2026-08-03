# =============================================================================
# CLUSTERING PAM (K-MEDOIDS) SOBRE MATRIZ RF
# =============================================================================
library(cluster)    # pam(), silhouette()
library(clValid)    # dunn(), connectivity()
library(openxlsx)

# --- Cargar configuración y funciones globales ---
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
dist_rf <- as.dist(matriz_cuadrada)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# FUNCIÓN PAM
# =============================================================================
metodo_PAM <- function(dist_arboles, matriz_arboles, res_calidad_clusters) {
  
  posibles_k <- seq(2, 15)
  
  for (i in posibles_k) {
    cat(sprintf("  Calculando k = %d...\n", i))
    
    # PAM acepta directamente objetos dist
    pam_resultado <- pam(dist_arboles, k = i, diss = TRUE)
    
    puntaje_Dunn         <- dunn(distance = matriz_arboles, pam_resultado$clustering)
    puntaje_Connectivity <- connectivity(distance = matriz_arboles, pam_resultado$clustering)
    # PAM ya calcula silhouette internamente; extraemos el promedio directamente
    puntaje_Silhouette   <- pam_resultado$silinfo$avg.width
    
    tmp <- c("PAM", i,
             round(puntaje_Dunn, 3),
             round(puntaje_Connectivity, 3),
             round(puntaje_Silhouette, 3))
    
    res_calidad_clusters <- rbind(res_calidad_clusters, tmp)
  }
  
  colnames(res_calidad_clusters) <- c("Method", "k", "Dunn", "Connectivity", "Silhouette")
  
  # Convertir columnas numéricas
  res_calidad_clusters$k            <- as.integer(res_calidad_clusters$k)
  res_calidad_clusters$Dunn         <- as.numeric(res_calidad_clusters$Dunn)
  res_calidad_clusters$Connectivity <- as.numeric(res_calidad_clusters$Connectivity)
  res_calidad_clusters$Silhouette   <- as.numeric(res_calidad_clusters$Silhouette)
  
  return(res_calidad_clusters)
}

# =============================================================================
# EJECUTAR PAM Y MEDIR TIEMPOS
# =============================================================================
cat("Ejecutando PAM para k = 2 a 15...\n")

res_calidad_clusters <- data.frame()
tiempo_inicio        <- proc.time()

res_calidad_clusters <- metodo_PAM(
  dist_arboles         = dist_rf,
  matriz_arboles       = matriz_cuadrada,
  res_calidad_clusters = res_calidad_clusters
)

tiempo_pam <- proc.time() - tiempo_inicio
cat("PAM completado.\n")
print(res_calidad_clusters)

# =============================================================================
# IDENTIFICAR K ÓPTIMO Y RE-EJECUTAR
# =============================================================================
k_optimo <- res_calidad_clusters$k[which.max(res_calidad_clusters$Silhouette)]
cat(sprintf("\nK óptimo según Silhouette: %d\n", k_optimo))
cat(sprintf("Silhouette promedio: %.3f\n", max(res_calidad_clusters$Silhouette)))

pam_optimo <- pam(dist_rf, k = k_optimo, diss = TRUE)

# Tabla de asignaciones: árbol → cluster → medoide del cluster
medoides_por_cluster <- pam_optimo$medoids[pam_optimo$clustering]

asignaciones_df <- data.frame(
  Arbol   = rownames(matriz_cuadrada),
  Cluster = pam_optimo$clustering,
  Medoide = medoides_por_cluster
)

# Tabla de medoides (representantes de cada cluster)
medoides_df <- data.frame(
  Cluster = seq_len(k_optimo),
  Medoide = pam_optimo$medoids,
  Tamano  = as.integer(table(pam_optimo$clustering))
)

cat("\nMedoides (árboles representativos por cluster):\n")
print(medoides_df)

# =============================================================================
# REPORTE EXCEL (tiempos + calidad + asignaciones + medoides)
# =============================================================================
cat("Generando reporte Excel...\n")

t_pam <- as.numeric(tiempo_pam)

tiempos_df <- data.frame(
  Proceso          = c("PAM (k=2 a 15)", "TOTAL"),
  Tiempo_Usuario_s = c(t_pam[1], t_pam[1]),
  Tiempo_Sistema_s = c(t_pam[2], t_pam[2]),
  Tiempo_Total_s   = c(t_pam[3], t_pam[3])
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Tiempos",
                              titulo_tabla = "Reporte de Tiempos de Ejecución (Segundos)",
                              datos        = tiempos_df,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Calidad_Clusters",
                              titulo_tabla = paste0("Índices de Calidad PAM (k=2 a 15) — K óptimo: ", k_optimo),
                              datos        = res_calidad_clusters,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Asignaciones_K_Optimo",
                              titulo_tabla = paste0("Asignaciones de Árboles — K óptimo = ", k_optimo),
                              datos        = asignaciones_df,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Medoides",
                              titulo_tabla = paste0("Medoides por Cluster — K óptimo = ", k_optimo),
                              datos        = medoides_df,
                              anchos_col   = "auto")

ruta_excel <- file.path(DIR_RESULTS, "pam_resultados.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)

cat("Resultados completos guardados en:", ruta_excel, "\n")
cat(sprintf("Resumen final: k óptimo = %d | Silhouette = %.3f\n",
            k_optimo,
            max(res_calidad_clusters$Silhouette)))