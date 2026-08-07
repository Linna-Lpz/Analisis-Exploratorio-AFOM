# =============================================================================
# 07_iterar_mstknn.R
# SUBDIVISIÓN ITERATIVA DE CLUSTERS MST-kNN CON CRITERIO BASADO EN SILHOUETTE
#
# DISEÑO:
#   - Punto de partida: asignaciones de MST-kNN (mejor escenario, k=5)
#   - Algoritmo de subdivisión: PAM sobre submatriz RF (no CLARA euclidiana)
#   - Criterio de parada por cluster: la Silhouette promedio del resultado
#     PAM(k=2) sobre la submatriz RF debe ser >= DELTA_SIL para aceptar el
#     corte. Si no hay ganancia real, el cluster se conserva intacto aunque
#     supere TAMANO_MAX. Esto garantiza que el número final de clusters refleja
#     estructura de datos, no aritmética de tamaño.
#   - TAMANO_MAX actúa como límite operacional adicional de seguridad.
#   - TAMANO_MIN_DURO: clusters con n < TAMANO_MIN_DURO no se subdividen y
#     se marcan como "insuficientes" para el enriquecimiento funcional.
# =============================================================================
library(cluster)    # pam(), silhouette()
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# PARÁMETROS
# =============================================================================
TAMANO_MAX      <- 200    # umbral operacional: clusters > N se intentan dividir
TAMANO_MIN_DURO <- 15     # clusters con n < N no se subdividen (ni se enriquecen)
DELTA_SIL       <- -999   # Forzamos subdivision ignorando el silhouette
K_SUBDIVISION   <- 2      # cortes binarios (k=2 en cada paso)
MAX_ITERACIONES <- 30     # tope de seguridad
SEED            <- 2

# =============================================================================
# 1. CARGAR INSUMOS
# =============================================================================
cat("=== CARGANDO INSUMOS ===\n")

# Matriz RF
ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")
if (!file.exists(ruta_cache_matriz)) stop("No se encontró matriz_rf.rds")
matriz_cuadrada <- readRDS(ruta_cache_matriz)
dist_rf         <- as.dist(matriz_cuadrada)
cat(sprintf("Matriz RF: %d x %d\n", nrow(matriz_cuadrada), ncol(matriz_cuadrada)))

# Asignaciones MST-kNN (mejor escenario)
ruta_mst_xlsx <- file.path(DIR_RESULTS, "mstknn_resultados2.xlsx")
if (!file.exists(ruta_mst_xlsx)) stop("No se encontró mstknn_resultados2.xlsx")

asig_inicial <- read.xlsx(ruta_mst_xlsx, sheet = "Asignaciones", startRow = 2,
                          colNames = TRUE)
asig_inicial$Arbol   <- as.character(asig_inicial$Arbol)
asig_inicial$Cluster <- as.integer(asig_inicial$Cluster)

# Filtrar NAs (nodos excluidos por MST-kNN)
asig_inicial <- asig_inicial[!is.na(asig_inicial$Cluster), ]

cat(sprintf("Árboles con asignación MST-kNN: %d de %d\n",
            nrow(asig_inicial), nrow(matriz_cuadrada)))
cat("Distribución inicial de clusters:\n")
print(sort(table(asig_inicial$Cluster), decreasing = TRUE))

# =============================================================================
# 2. FUNCIÓN DE SUBDIVISIÓN CON PAM SOBRE dist_rf
#    Devuelve: list(asig = vector nombrado, sil_avg = silhouette promedio)
#    sil_avg = NA si no se pudo subdividir
# =============================================================================
subdividir_con_pam <- function(ids_arboles, matriz_completa,
                               k = 2, seed = 2, min_n = 3) {

  n <- length(ids_arboles)

  if (n < max(k + 1, min_n)) {
    return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA,
                razon = "n_insuficiente"))
  }

  # Verificar que los ids estén en la matriz
  ids_validos <- ids_arboles[ids_arboles %in% rownames(matriz_completa)]
  if (length(ids_validos) < k + 1) {
    return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA,
                razon = "ids_no_en_matriz"))
  }

  submatriz <- matriz_completa[ids_validos, ids_validos]
  dist_sub  <- as.dist(submatriz)

  set.seed(seed)
  pam_res <- tryCatch(
    pam(dist_sub, k = k, diss = TRUE),
    error = function(e) {
      cat(sprintf("    ERROR pam: %s\n", e$message))
      NULL
    }
  )

  if (is.null(pam_res)) {
    return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA,
                razon = "pam_error"))
  }

  sil_obj <- silhouette(pam_res$clustering, dist_sub)
  sil_avg <- mean(sil_obj[, 3])

  list(
    asig    = setNames(as.integer(pam_res$clustering), ids_validos),
    sil_avg = round(sil_avg, 4),
    razon   = "ok"
  )
}

# =============================================================================
# 3. BUCLE ITERATIVO DE SUBDIVISIÓN
# =============================================================================
cat("\n=== INICIANDO SUBDIVISIÓN ITERATIVA CON CRITERIO DE SILHOUETTE ===\n")
cat(sprintf("TAMANO_MAX     : %d (umbral operacional)\n", TAMANO_MAX))
cat(sprintf("TAMANO_MIN_DURO: %d (mínimo absoluto para subdivisión)\n", TAMANO_MIN_DURO))
cat(sprintf("DELTA_SIL      : %.3f (ganancia mínima de Silhouette para aceptar corte)\n", DELTA_SIL))
cat(sprintf("k por corte    : %d\n\n", K_SUBDIVISION))

# Estado inicial
estado_actual <- setNames(asig_inicial$Cluster, asig_inicial$Arbol)

# Log completo de subdivisiones
log_subdivisiones <- data.frame()
clusters_rechazados <- data.frame()   # clusters grandes que no se subdividieron por criterio Sil
iteracion         <- 0L
tiempo_total      <- proc.time()

repeat {
  iteracion <- iteracion + 1L

  if (iteracion > MAX_ITERACIONES) {
    cat(sprintf("AVISO: Se alcanzó el límite de %d iteraciones.\n", MAX_ITERACIONES))
    break
  }

  tamanos_actuales <- table(estado_actual)
  # Candidatos: clusters con n > TAMANO_MAX Y n >= TAMANO_MIN_DURO
  clusters_grandes <- as.integer(names(tamanos_actuales[
    tamanos_actuales > TAMANO_MAX & tamanos_actuales >= TAMANO_MIN_DURO
  ]))

  if (length(clusters_grandes) == 0) {
    cat(sprintf("Iteración %d: ningún cluster supera TAMANO_MAX = %d. Finalizado.\n",
                iteracion, TAMANO_MAX))
    break
  }

  cat(sprintf("--- Iteración %d: %d cluster(s) candidatos a subdividir ---\n",
              iteracion, length(clusters_grandes)))

  proximo_id <- max(estado_actual) + 1L
  nuevo_estado <- estado_actual
  hubo_cambio <- FALSE

  for (cid in clusters_grandes) {

    arboles_cid    <- names(estado_actual[estado_actual == cid])
    n_cid          <- length(arboles_cid)

    cat(sprintf("  Cluster %d (n=%d):\n", cid, n_cid))

    # Calcular Silhouette del cluster padre sobre su submatriz
    # (como referencia: si es un cluster cohesivo k=1, sil=0)
    # Para k=2, la Silhouette del sub-resultado debe ser >= DELTA_SIL
    resultado_sub <- subdividir_con_pam(
      ids_arboles     = arboles_cid,
      matriz_completa = matriz_cuadrada,
      k               = K_SUBDIVISION,
      seed            = SEED,
      min_n           = TAMANO_MIN_DURO
    )

    if (is.na(resultado_sub$sil_avg) || resultado_sub$razon != "ok") {
      cat(sprintf("    → No subdividible (%s). Conservando.\n", resultado_sub$razon))
      clusters_rechazados <- rbind(clusters_rechazados, data.frame(
        Iteracion      = iteracion,
        Cluster        = cid,
        N              = n_cid,
        Sil_Propuesta  = NA,
        Decision       = resultado_sub$razon
      ))
      next
    }

    sil_propuesta <- resultado_sub$sil_avg

    if (sil_propuesta < DELTA_SIL) {
      cat(sprintf("    → Silhouette propuesta = %.4f < DELTA_SIL = %.3f. Conservando.\n",
                  sil_propuesta, DELTA_SIL))
      clusters_rechazados <- rbind(clusters_rechazados, data.frame(
        Iteracion      = iteracion,
        Cluster        = cid,
        N              = n_cid,
        Sil_Propuesta  = sil_propuesta,
        Decision       = "sil_insuficiente"
      ))
      next
    }

    # Aceptar subdivisión
    cat(sprintf("    → Silhouette propuesta = %.4f >= %.3f. ACEPTADA.\n",
                sil_propuesta, DELTA_SIL))

    asig_nueva <- resultado_sub$asig
    sub_ids    <- sort(unique(asig_nueva))

    mapa_ids <- setNames(
      c(cid, seq(proximo_id, proximo_id + length(sub_ids) - 2L)),
      sub_ids
    )
    asig_reetiquetada <- mapa_ids[as.character(asig_nueva)]

    tamanos_nuevos <- table(asig_reetiquetada)
    cat(sprintf("    → %d sub-clusters: %s\n",
                length(tamanos_nuevos),
                paste(paste0("n=", sort(as.integer(tamanos_nuevos), decreasing = TRUE)),
                      collapse = ", ")))

    ids_validos_cid <- names(asig_nueva)
    nuevo_estado[ids_validos_cid] <- asig_reetiquetada
    proximo_id <- proximo_id + length(sub_ids) - 1L
    hubo_cambio <- TRUE

    log_subdivisiones <- rbind(log_subdivisiones, data.frame(
      Iteracion           = iteracion,
      Cluster_Original    = cid,
      N_Original          = n_cid,
      Sil_Sub             = round(sil_propuesta, 4),
      N_SubClusters       = length(tamanos_nuevos),
      Tamanos_SubClusters = paste(sort(as.integer(tamanos_nuevos), decreasing = TRUE),
                                  collapse = " | ")
    ))
  }

  estado_actual <- nuevo_estado

  if (!hubo_cambio) {
    cat(sprintf("Iteración %d: ningún cluster aceptó subdivisión (todos bajo DELTA_SIL). Finalizando.\n",
                iteracion))
    break
  }
}

tiempo_total <- as.numeric(proc.time() - tiempo_total)[3]

# =============================================================================
# 4. TABLA FINAL DE ASIGNACIONES (reetiquetado secuencial)
# =============================================================================
cat("\n=== CONSTRUYENDO TABLA FINAL ===\n")

ids_originales <- sort(unique(estado_actual))
mapa_final     <- setNames(seq_along(ids_originales), ids_originales)
cluster_final  <- mapa_final[as.character(estado_actual)]

asig_final_df <- data.frame(
  Arbol         = names(estado_actual),
  Cluster_Final = as.integer(cluster_final),
  stringsAsFactors = FALSE
)
asig_final_df <- asig_final_df[order(asig_final_df$Cluster_Final), ]
rownames(asig_final_df) <- NULL

# Distribución de tamaños
tamanos_finales <- as.data.frame(table(Cluster = asig_final_df$Cluster_Final))
tamanos_finales$Cluster <- as.integer(as.character(tamanos_finales$Cluster))
colnames(tamanos_finales)[2] <- "Tamano"
tamanos_finales <- tamanos_finales[order(tamanos_finales$Tamano, decreasing = TRUE), ]

n_clusters_final  <- length(unique(asig_final_df$Cluster_Final))
n_bajo_15         <- sum(tamanos_finales$Tamano < 15)
n_bajo_50         <- sum(tamanos_finales$Tamano < 50)

cat(sprintf("Clusters iniciales (MST-kNN)      : %d\n",
            length(unique(asig_inicial$Cluster))))
cat(sprintf("Clusters finales                  : %d\n", n_clusters_final))
cat(sprintf("Árboles cubiertos                 : %d de %d\n",
            nrow(asig_final_df), nrow(matriz_cuadrada)))
cat(sprintf("Iteraciones realizadas            : %d\n", iteracion - 1L))
cat(sprintf("Tiempo total                      : %.1f s\n", tiempo_total))
cat(sprintf("Clusters con n < 15 (insuficientes): %d de %d (%.1f%%)\n",
            n_bajo_15, n_clusters_final,
            100 * n_bajo_15 / n_clusters_final))
cat(sprintf("Clusters con n < 50               : %d de %d (%.1f%%)\n",
            n_bajo_50, n_clusters_final,
            100 * n_bajo_50 / n_clusters_final))
cat("\nResumen estadístico de tamaños:\n")
print(summary(tamanos_finales$Tamano))

# =============================================================================
# 5. SILHOUETTE GLOBAL DE LA PARTICIÓN FINAL (sobre dist_rf)
# =============================================================================
cat("\n=== CALCULANDO SILHOUETTE GLOBAL DE LA PARTICIÓN FINAL ===\n")
cat("(sobre disimilitud RF — puede tardar varios minutos)\n")

# Usar solo los nodos con asignación
ids_con_asig   <- asig_final_df$Arbol
ids_con_asig   <- ids_con_asig[ids_con_asig %in% rownames(matriz_cuadrada)]
submat_final   <- matriz_cuadrada[ids_con_asig, ids_con_asig]
dist_final     <- as.dist(submat_final)
asig_vec_final <- asig_final_df$Cluster_Final[
  match(ids_con_asig, asig_final_df$Arbol)
]

tiempo_sil <- proc.time()
sil_final  <- tryCatch(
  silhouette(asig_vec_final, dist_final),
  error = function(e) { cat("ERROR silhouette global:", e$message, "\n"); NULL }
)
tiempo_sil <- as.numeric(proc.time() - tiempo_sil)[3]

sil_global_avg <- if (!is.null(sil_final)) round(mean(sil_final[, 3]), 4) else NA
cat(sprintf("Silhouette global (partición final, %d clusters): %.4f\n",
            n_clusters_final, sil_global_avg))
cat(sprintf("Tiempo cálculo Silhouette: %.1f s\n", tiempo_sil))

# Silhouette por cluster
if (!is.null(sil_final)) {
  sil_por_cluster <- tapply(sil_final[, 3], sil_final[, 1], mean)
  tamanos_finales$Silhouette_Cluster <- round(
    sil_por_cluster[as.character(tamanos_finales$Cluster)], 4
  )
} else {
  tamanos_finales$Silhouette_Cluster <- NA
}

# Etiqueta de suficiencia para enriquecimiento
tamanos_finales$Apto_Enriquecimiento <- tamanos_finales$Tamano >= 15

# =============================================================================
# 6. GUARDAR ASIGNACIONES COMO .rds PARA SCRIPTS POSTERIORES
# =============================================================================
ruta_cache_asig <- file.path(DIR_CACHE, "mstknn_iter_asignaciones.rds")
saveRDS(asig_final_df, ruta_cache_asig)
cat(sprintf("Asignaciones guardadas en: %s\n", ruta_cache_asig))

# =============================================================================
# 7. EXPORTAR EXCEL COMPLETO
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

# Parámetros
parametros_df <- data.frame(
  Parametro = c("tamano_max", "tamano_min_duro", "delta_sil", "k_subdivision",
                "max_iteraciones", "iteraciones_realizadas", "seed",
                "clusters_iniciales_mstknn", "clusters_finales",
                "arboles_cubiertos", "arboles_totales",
                "n_clusters_bajo_15", "n_clusters_bajo_50",
                "silhouette_global_final", "tiempo_total_s"),
  Valor     = c(TAMANO_MAX, TAMANO_MIN_DURO, DELTA_SIL, K_SUBDIVISION,
                MAX_ITERACIONES, iteracion - 1L, SEED,
                length(unique(asig_inicial$Cluster)), n_clusters_final,
                nrow(asig_final_df), nrow(matriz_cuadrada),
                n_bajo_15, n_bajo_50,
                sil_global_avg, round(tiempo_total, 1))
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(
  wb           = wb,
  nombre_hoja  = "Asignaciones_Finales",
  titulo_tabla = paste0("Asignaciones MST-kNN Iterativo — ",
                        n_clusters_final, " clusters"),
  datos        = asig_final_df,
  anchos_col   = "auto"
)

wb <- agregar_hoja_formateada(
  wb           = wb,
  nombre_hoja  = "Tamanos_y_Silhouette",
  titulo_tabla = paste0("Tamaño y Silhouette por Cluster — Sil global = ",
                        sil_global_avg),
  datos        = tamanos_finales,
  anchos_col   = "auto"
)

if (nrow(log_subdivisiones) > 0) {
  wb <- agregar_hoja_formateada(
    wb           = wb,
    nombre_hoja  = "Log_Subdivisiones_Aceptadas",
    titulo_tabla = "Log de Subdivisiones Aceptadas (Sil >= DELTA_SIL)",
    datos        = log_subdivisiones,
    anchos_col   = "auto"
  )
}

if (nrow(clusters_rechazados) > 0) {
  wb <- agregar_hoja_formateada(
    wb           = wb,
    nombre_hoja  = "Log_Rechazados",
    titulo_tabla = "Clusters grandes conservados (Sil < DELTA_SIL o n insuficiente)",
    datos        = clusters_rechazados,
    anchos_col   = "auto"
  )
}

wb <- agregar_hoja_formateada(
  wb           = wb,
  nombre_hoja  = "Parametros",
  titulo_tabla = "Parámetros del Proceso de Subdivisión Iterativa MST-kNN",
  datos        = parametros_df,
  anchos_col   = c(35, 25)
)

ruta_excel <- file.path(DIR_RESULTS, "mstknn_subdivision_iterativa.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Excel guardado:", ruta_excel, "\n")

cat(sprintf("\n=== RESUMEN FINAL ===\n"))
cat(sprintf("  Clusters MST-kNN iniciales : %d\n", length(unique(asig_inicial$Cluster))))
cat(sprintf("  Clusters finales           : %d\n", n_clusters_final))
cat(sprintf("  Silhouette global (RF)     : %.4f\n", sil_global_avg))
cat(sprintf("  Clusters aptos (n >= 15)   : %d\n",
            sum(tamanos_finales$Apto_Enriquecimiento)))
cat(sprintf("  Clusters insuficientes     : %d\n", n_bajo_15))
cat("\n=== COMPLETADO: 07_iterar_mstknn.R ===\n")
