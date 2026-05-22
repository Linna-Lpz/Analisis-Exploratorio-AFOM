# =============================================================================
# CLUSTERING CLARA (CLUSTERING LARGE APPLICATIONS) SOBRE MATRIZ RF
# =============================================================================
library(cluster)    # clara(), silhouette()
library(clValid)    # dunn(), connectivity()
library(openxlsx)

# --- Cargar configuración y funciones globales ---
source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# LEER MATRIZ DESDE CACHÉ .rds
# =============================================================================
ruta_cache_matriz <- file.path(DIR_PROCESSED, "cache", "matriz_rf.rds")
if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché. Ejecuta primero el script de cálculo RF.")
}
cat("Cargando matriz RF desde caché...\n")
matriz_cuadrada <- readRDS(ruta_cache_matriz)
dist_rf <- as.dist(matriz_cuadrada)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# FUNCIÓN CLARA
# =============================================================================
# NOTA: clara() no acepta objetos dist directamente — requiere la matriz
# numérica y metric = "euclidean" o "manhattan". Para trabajar con distancias
# RF precomputadas se usa el argumento keep.data = FALSE y se pasa la
# matriz cuadrada directamente como insumo de distancias.
# dunn() y connectivity() siguen recibiendo la matriz original.
# =============================================================================
metodo_CLARA <- function(matriz_arboles, res_calidad_clusters,
                         muestras = 50, sampsize = NULL) {
  
  posibles_k <- seq(2, 15)
  n          <- nrow(matriz_arboles)
  
  # sampsize por defecto: min(n, 40 + 2*k) por iteración — valor recomendado
  # en Kaufman & Rousseeuw. Si el usuario pasa un valor fijo se respeta.
  usar_sampsize_fijo <- !is.null(sampsize)
  
  set.seed(2)
  
  for (i in posibles_k) {
    cat(sprintf("  Calculando k = %d...\n", i))
    
    ss <- if (usar_sampsize_fijo) sampsize else min(n, 40 + 2 * i)
    
    clara_resultado <- clara(
      x         = matriz_arboles,
      k         = i,
      metric    = "euclidean",   # distancias ya precomputadas en la matriz
      samples   = muestras,      # número de submuestras aleatorias
      sampsize  = ss,            # tamaño de cada submuestra
      keep.data = FALSE,         # no duplicar la matriz en memoria
      rngR      = TRUE           # usar RNG de R (respeta set.seed)
    )
    
    puntaje_Dunn         <- dunn(distance = matriz_arboles, clara_resultado$clustering)
    puntaje_Connectivity <- connectivity(distance = matriz_arboles, clara_resultado$clustering)
    puntaje_Silhouette   <- clara_resultado$silinfo$avg.width
    
    tmp <- c("CLARA", i,
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
# EJECUTAR CLARA Y MEDIR TIEMPOS
# =============================================================================
cat("Ejecutando CLARA para k = 2 a 15...\n")
cat(sprintf("Configuración: 50 submuestras | sampsize = 40 + 2*k (dinámico)\n\n"))

res_calidad_clusters <- data.frame()
tiempo_inicio        <- proc.time()

res_calidad_clusters <- metodo_CLARA(
  matriz_arboles       = matriz_cuadrada,
  res_calidad_clusters = res_calidad_clusters,
  muestras             = 50
)

tiempo_clara <- proc.time() - tiempo_inicio
cat("\nCLARA completado.\n")
print(res_calidad_clusters)

# =============================================================================
# IDENTIFICAR K ÓPTIMO Y RE-EJECUTAR
# =============================================================================
k_optimo <- res_calidad_clusters$k[which.max(res_calidad_clusters$Silhouette)]
cat(sprintf("\nK óptimo según Silhouette: %d\n", k_optimo))
cat(sprintf("Silhouette promedio:        %.3f\n", max(res_calidad_clusters$Silhouette)))

set.seed(2)
n <- nrow(matriz_cuadrada)

clara_optimo <- clara(
  x         = matriz_cuadrada,
  k         = k_optimo,
  metric    = "euclidean",
  samples   = 50,
  sampsize  = min(n, 40 + 2 * k_optimo),
  keep.data = FALSE,
  rngR      = TRUE
)

# Tabla de asignaciones: árbol → cluster → medoide del cluster
medoides_por_arbol <- clara_optimo$medoids[clara_optimo$clustering, ]

asignaciones_df <- data.frame(
  Arbol   = rownames(matriz_cuadrada),
  Cluster = clara_optimo$clustering,
  Medoide = rownames(clara_optimo$medoids)[clara_optimo$clustering]
)

# Tabla de medoides con tamaño de cada cluster
medoides_df <- data.frame(
  Cluster  = seq_len(k_optimo),
  Medoide  = rownames(clara_optimo$medoids),
  Tamano   = as.integer(table(clara_optimo$clustering)),
  Silhouette_Cluster = round(
    tapply(clara_optimo$silinfo$widths[, "sil_width"],
           clara_optimo$silinfo$widths[, "cluster"],
           mean), 3
  )
)

cat("\nMedoides (árboles representativos por cluster):\n")
print(medoides_df)

# =============================================================================
# TABLA COMPARATIVA K-MEANS vs PAM vs CLARA
# =============================================================================
# Leer resultados previos si existen
ruta_kmeans <- file.path(DIR_RESULTS, "kmeans_calidad_clusters.csv")
ruta_pam    <- file.path(DIR_RESULTS, "pam_calidad_clusters.csv")

comparativa_disponible <- file.exists(ruta_kmeans) && file.exists(ruta_pam)

if (comparativa_disponible) {
  cat("\nGenerando tabla comparativa de métodos...\n")
  
  kmeans_df <- read.table(ruta_kmeans, sep = ";", header = TRUE)
  pam_df    <- read.table(ruta_pam,    sep = ";", header = TRUE)
  
  comparativa_df <- rbind(kmeans_df, pam_df, res_calidad_clusters)
  
  # Resumen: mejor k por método según Silhouette
  resumen_df <- do.call(rbind, lapply(
    split(comparativa_df, comparativa_df$Method),
    function(df) {
      idx <- which.max(df$Silhouette)
      data.frame(
        Metodo             = df$Method[idx],
        K_Optimo           = df$k[idx],
        Dunn_Optimo        = df$Dunn[idx],
        Connectivity_Optimo = df$Connectivity[idx],
        Silhouette_Optimo  = df$Silhouette[idx]
      )
    }
  ))
  rownames(resumen_df) <- NULL
  
  cat("\nResumen comparativo (mejor k por método):\n")
  print(resumen_df)
  
} else {
  cat("\nNo se encontraron resultados previos de K-Means o PAM para comparar.\n")
  comparativa_df <- res_calidad_clusters
  resumen_df     <- NULL
}


# =============================================================================
# REPORTE EXCEL
# =============================================================================
cat("Generando reporte Excel...\n")

t_clara <- as.numeric(tiempo_clara)

tiempos_df <- data.frame(
  Proceso          = c("CLARA (k=2 a 15)", "TOTAL"),
  Tiempo_Usuario_s = c(t_clara[1], t_clara[1]),
  Tiempo_Sistema_s = c(t_clara[2], t_clara[2]),
  Tiempo_Total_s   = c(t_clara[3], t_clara[3])
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Tiempos",
                              titulo_tabla = "Reporte de Tiempos de Ejecución (Segundos)",
                              datos        = tiempos_df,
                              anchos_col   = "auto")

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Calidad_Clusters",
                              titulo_tabla = paste0("Índices de Calidad CLARA (k=2 a 15) — K óptimo: ", k_optimo),
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

# Hoja comparativa si está disponible
if (comparativa_disponible) {
  wb <- agregar_hoja_formateada(wb           = wb,
                                nombre_hoja  = "Comparativa_Metodos",
                                titulo_tabla = "Comparativa K-Means vs PAM vs CLARA (todos los k)",
                                datos        = comparativa_df,
                                anchos_col   = "auto")
  
  wb <- agregar_hoja_formateada(wb           = wb,
                                nombre_hoja  = "Resumen_Comparativo",
                                titulo_tabla = "Mejor K por Método — Resumen Ejecutivo",
                                datos        = resumen_df,
                                anchos_col   = "auto")
}

ruta_excel <- file.path(DIR_RESULTS, "clara_resultados.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)

cat("Resultados completos guardados en:", ruta_excel, "\n")
cat(sprintf("Resumen final: k óptimo = %d | Silhouette = %.3f\n",
            k_optimo,
            max(res_calidad_clusters$Silhouette)))