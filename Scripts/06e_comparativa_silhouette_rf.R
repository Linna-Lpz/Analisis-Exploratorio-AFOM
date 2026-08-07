# =============================================================================
# 06e_comparativa_silhouette_rf.R
# PROPÓSITO: Recalcular Silhouette de los cuatro algoritmos de clustering
#            sobre la MISMA disimilitud RF, para hacer la comparación
#            de la Tabla 3.2 metodológicamente válida.
#
# PROBLEMA IDENTIFICADO:
#   - CLARA: silhouette calculada internamente con metric="euclidean" sobre
#             las filas de la matriz RF → espacio euclidiano, NO disimilitud RF.
#   - K-Means: silhouette calculada sobre la matriz cuadrada como si fuera
#               distancia, pero kmeans() operó en espacio euclidiano.
#   - PAM: silhouette calculada directamente sobre dist_rf. ✓
#   - MST-kNN: silhouette calculada sobre submatriz RF. ✓
#
# SOLUCIÓN:
#   Tomar las asignaciones de k óptimo de cada algoritmo y recalcular
#   silhouette() sobre dist_rf para todos por igual.
# =============================================================================
library(cluster)    # silhouette()
library(clValid)    # dunn(), connectivity()
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# 1. CARGAR MATRIZ RF Y CONSTRUIR OBJETO dist
# =============================================================================
cat("Cargando matriz RF desde caché...\n")
ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")
if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché.")
}
matriz_cuadrada <- readRDS(ruta_cache_matriz)
dist_rf         <- as.dist(matriz_cuadrada)
n               <- nrow(matriz_cuadrada)
cat(sprintf("Matriz cargada: %d x %d\n", n, n))

# =============================================================================
# 2. LEER ASIGNACIONES ÓPTIMAS DE CADA ALGORITMO
#    (usamos las ya guardadas; k óptimo según Silhouette original = k=2 en todos)
# =============================================================================

# ---- 2a. K-Means (desde caché .rds) ----------------------------------------
cat("\nLeyendo asignaciones K-Means...\n")
ruta_kmeans_asig <- file.path(DIR_CACHE, "kmeans_asignaciones.rds")
if (!file.exists(ruta_kmeans_asig)) {
  stop("No se encontró kmeans_asignaciones.rds. Ejecuta primero 06a_clustering_kmeans.R")
}
kmeans_asig_df <- readRDS(ruta_kmeans_asig)
# Asegurarse de que el orden de filas coincida con la matriz
kmeans_asig_df <- kmeans_asig_df[match(rownames(matriz_cuadrada), kmeans_asig_df$Arbol), ]
asig_kmeans    <- kmeans_asig_df$Cluster
k_kmeans       <- length(unique(asig_kmeans))
cat(sprintf("  K-Means: k=%d, n=%d\n", k_kmeans, sum(!is.na(asig_kmeans))))

# ---- 2b. PAM (desde Excel de resultados) ------------------------------------
cat("Leyendo asignaciones PAM...\n")
ruta_pam_xlsx <- file.path(DIR_RESULTS, "pam_resultados.xlsx")
if (!file.exists(ruta_pam_xlsx)) {
  stop("No se encontró pam_resultados.xlsx. Ejecuta primero 06b_clustering_pam.R")
}
pam_asig_df <- read.xlsx(ruta_pam_xlsx, sheet = "Asignaciones_K_Optimo", startRow = 2)
pam_asig_df <- pam_asig_df[match(rownames(matriz_cuadrada), pam_asig_df$Arbol), ]
asig_pam    <- pam_asig_df$Cluster
k_pam       <- length(unique(asig_pam[!is.na(asig_pam)]))
cat(sprintf("  PAM: k=%d, n=%d\n", k_pam, sum(!is.na(asig_pam))))

# ---- 2c. CLARA (desde Excel de resultados) ----------------------------------
cat("Leyendo asignaciones CLARA...\n")
ruta_clara_xlsx <- file.path(DIR_RESULTS, "clara_resultados.xlsx")
if (!file.exists(ruta_clara_xlsx)) {
  stop("No se encontró clara_resultados.xlsx. Ejecuta primero 06c_clustering_clara.R")
}
clara_asig_df <- read.xlsx(ruta_clara_xlsx, sheet = "Asignaciones_K_Optimo", startRow = 2)
clara_asig_df <- clara_asig_df[match(rownames(matriz_cuadrada), clara_asig_df$Arbol), ]
asig_clara    <- clara_asig_df$Cluster
k_clara       <- length(unique(asig_clara[!is.na(asig_clara)]))
cat(sprintf("  CLARA: k=%d, n=%d\n", k_clara, sum(!is.na(asig_clara))))

# ---- 2d. MST-kNN (desde Excel de resultados) --------------------------------
cat("Leyendo asignaciones MST-kNN...\n")
ruta_mst_xlsx <- file.path(DIR_RESULTS, "mstknn_resultados2.xlsx")
if (!file.exists(ruta_mst_xlsx)) {
  stop("No se encontró mstknn_resultados2.xlsx. Ejecuta primero 06d_clustering_mstknn.R")
}
mst_asig_df <- read.xlsx(ruta_mst_xlsx, sheet = "Asignaciones", startRow = 2)

# MST-kNN puede no cubrir todos los nodos; alinear con la matriz
mst_asig_df <- mst_asig_df[match(rownames(matriz_cuadrada), mst_asig_df$Arbol), ]
asig_mst    <- mst_asig_df$Cluster
n_mst_validos <- sum(!is.na(asig_mst))
k_mst         <- length(unique(asig_mst[!is.na(asig_mst)]))
cat(sprintf("  MST-kNN: k=%d, n=%d (de %d totales)\n", k_mst, n_mst_validos, n))

# =============================================================================
# 3. RECALCULAR SILHOUETTE SOBRE dist_rf PARA CADA ALGORITMO
# =============================================================================
cat("\n=== Recalculando Silhouette sobre dist_rf ===\n\n")

calcular_sil_rf <- function(nombre_alg, asignaciones, dist_obj, matriz_cuad) {
  
  # Filtrar NAs (para MST-kNN que puede no cubrir todos los nodos)
  mask      <- !is.na(asignaciones)
  asig_val  <- as.integer(asignaciones[mask])
  n_val     <- sum(mask)
  k_val     <- length(unique(asig_val))
  
  cat(sprintf("Calculando %s: k=%d, n=%d\n", nombre_alg, k_val, n_val))
  
  if (k_val < 2) {
    cat(sprintf("  AVISO: k=%d, silhouette no aplicable.\n", k_val))
    return(list(
      Algoritmo   = nombre_alg,
      k           = k_val,
      n           = n_val,
      Silhouette  = NA,
      Dunn        = NA,
      Connectivity = NA
    ))
  }
  
  # Si todos los nodos están presentes, usar dist_obj directamente
  if (n_val == attr(dist_obj, "Size")) {
    dist_uso <- dist_obj
  } else {
    # Subconjunto de la matriz cuadrada para nodos con asignación
    idx_val  <- which(mask)
    sub_mat  <- matriz_cuad[idx_val, idx_val]
    dist_uso <- as.dist(sub_mat)
  }
  
  # Silhouette
  sil_obj <- tryCatch(
    silhouette(asig_val, dist_uso),
    error = function(e) { cat("  ERROR silhouette:", e$message, "\n"); NULL }
  )
  sil_avg <- if (!is.null(sil_obj)) round(mean(sil_obj[, 3]), 3) else NA
  cat(sprintf("  Silhouette (sobre dist_rf): %.3f\n", sil_avg))
  
  # Dunn (sobre la matriz cuadrada correspondiente)
  if (n_val == attr(dist_obj, "Size")) {
    mat_uso <- matriz_cuad
  } else {
    mat_uso <- sub_mat
  }
  puntaje_Dunn <- tryCatch(
    round(dunn(distance = mat_uso, asig_val), 3),
    error = function(e) { cat("  ERROR Dunn:", e$message, "\n"); NA }
  )
  cat(sprintf("  Dunn (sobre dist_rf):      %.3f\n", puntaje_Dunn))
  
  # Connectivity (sobre la matriz cuadrada)
  puntaje_Conn <- tryCatch(
    round(connectivity(distance = mat_uso, asig_val), 3),
    error = function(e) { cat("  ERROR Connectivity:", e$message, "\n"); NA }
  )
  cat(sprintf("  Connectivity:              %.3f\n\n", puntaje_Conn))
  
  list(
    Algoritmo    = nombre_alg,
    k            = k_val,
    n            = n_val,
    Silhouette   = sil_avg,
    Dunn         = puntaje_Dunn,
    Connectivity = puntaje_Conn
  )
}

res_kmeans <- calcular_sil_rf("K-Means",  asig_kmeans, dist_rf, matriz_cuadrada)
res_pam    <- calcular_sil_rf("PAM",      asig_pam,    dist_rf, matriz_cuadrada)
res_clara  <- calcular_sil_rf("CLARA",    asig_clara,  dist_rf, matriz_cuadrada)
res_mst    <- calcular_sil_rf("MST-kNN",  asig_mst,    dist_rf, matriz_cuadrada)

# =============================================================================
# 4. CONSTRUIR TABLA COMPARATIVA FINAL
# =============================================================================
tabla_comparativa <- do.call(rbind, lapply(
  list(res_clara, res_mst, res_kmeans, res_pam),
  function(r) as.data.frame(r, stringsAsFactors = FALSE)
))
rownames(tabla_comparativa) <- NULL

cat("=== TABLA COMPARATIVA (Silhouette sobre disimilitud RF) ===\n")
print(tabla_comparativa)

# Identificar mejor algoritmo según Silhouette sobre RF
idx_mejor <- which.max(tabla_comparativa$Silhouette)
cat(sprintf("\nMejor algoritmo según Silhouette (sobre dist_rf): %s (Sil=%.3f)\n",
            tabla_comparativa$Algoritmo[idx_mejor],
            tabla_comparativa$Silhouette[idx_mejor]))

# =============================================================================
# 5. COMPARACIÓN CON VALORES ORIGINALES REPORTADOS EN EL LATEX
# =============================================================================
originales_df <- data.frame(
  Algoritmo         = c("CLARA (k=2)", "MST-kNN (k=5)", "K-Means (k=2)", "PAM (k=2)"),
  Sil_Original      = c(0.577, 0.234, 0.165, 0.007),
  Espacio_Original  = c("Euclidiano (filas RF)", "RF (submatriz)", 
                        "Euclidiano (filas RF)", "RF directo"),
  stringsAsFactors  = FALSE
)

cat("\n=== COMPARACIÓN CON VALORES ORIGINALMENTE REPORTADOS ===\n")
print(originales_df)
cat("\nNota: Los valores originales de CLARA y K-Means NO son comparables\n")
cat("      con PAM porque usan espacios métricos distintos.\n")
cat("      Los valores recalculados arriba son comparables entre sí.\n")

# =============================================================================
# 6. GUARDAR RESULTADOS
# =============================================================================
# Guardar como .rds para referencia posterior
ruta_rds <- file.path(DIR_CACHE, "comparativa_silhouette_rf.rds")
saveRDS(tabla_comparativa, ruta_rds)
cat(sprintf("\nResultados guardados en: %s\n", ruta_rds))

# Guardar como Excel
wb <- createWorkbook()

wb <- agregar_hoja_formateada(
  wb           = wb,
  nombre_hoja  = "Silhouette_RF_Comparable",
  titulo_tabla = "Coeficiente Silhouette recalculado sobre disimilitud RF (comparación válida)",
  datos        = tabla_comparativa,
  anchos_col   = "auto"
)

wb <- agregar_hoja_formateada(
  wb           = wb,
  nombre_hoja  = "Valores_Originales",
  titulo_tabla = "Valores originalmente reportados (NO comparables entre sí)",
  datos        = originales_df,
  anchos_col   = "auto"
)

ruta_excel <- file.path(DIR_RESULTS, "comparativa_silhouette_rf.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat(sprintf("Reporte Excel guardado en: %s\n", ruta_excel))

cat("\n=== FIN: 06e_comparativa_silhouette_rf.R ===\n")
