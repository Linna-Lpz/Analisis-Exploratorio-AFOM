# =============================================================================
# SUBDIVISIÓN ITERATIVA DE CLUSTERS GRANDES — K-MEANS (subdivisión con PAM sobre RF)
#
# ACTUALIZACIóN: La subdivisión interna ya no usa K-Means euclidiano sino PAM
# directamente sobre la submatriz de disimilitud Robinson-Foulds. Se añade un
# criterio de parada por ganancia de Silhouette (DELTA_SIL): un cluster solo
# se subdivide si la Silhouette promedio del resultado PAM(k=2) sobre la
# submatriz RF es >= DELTA_SIL. Esto evita fraccionar clusters sin estructura
# real en el espacio de referencia correcto.
# =============================================================================
library(cluster)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# PARÁMETROS AJUSTABLES
# =============================================================================
TAMANO_MAX      <- 200    # clusters con más de N árboles son candidatos a subdividir
TAMANO_MIN_DURO <- 15     # clusters con n < N no se subdividen (ni se enriquecen)
TAMANO_MIN_META <- 50     # tamaño mínimo deseado (informativo)
DELTA_SIL       <- 0.01   # ganancia mínima de Silhouette (sobre RF) para aceptar corte
K_SUBDIVISION   <- 2      # k fijo para cada corte binario
MAX_ITERACIONES <- 30     # tope de seguridad contra loops infinitos
SEED            <- 2

# Fuente de asignaciones K-Means: número (10, 15) o "optimo"
K_FUENTE <- "optimo"

# =============================================================================
# 1. CARGAR INSUMOS
# =============================================================================
cat("=== CARGANDO INSUMOS ===\n")

ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")
if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché. Ejecuta primero el script de cálculo RF.")
}
matriz_cuadrada <- readRDS(ruta_cache_matriz)
cat("Matriz RF:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# dist_rf global (se reutilizará para la Silhouette final)
dist_rf <- as.dist(matriz_cuadrada)

ruta_kmeans <- file.path(DIR_RESULTS, "kmeans_resultados.xlsx")
if (!file.exists(ruta_kmeans)) {
  stop("No se encontró kmeans_resultados.xlsx. Ejecuta primero el script de K-Means.")
}

if (identical(K_FUENTE, "optimo")) {
  hoja_fuente   <- "Asignaciones_K_Optimo"
  titulo_fuente <- "k óptimo (K-Means)"
} else {
  hoja_fuente   <- paste0("Asignaciones_K_", K_FUENTE)
  titulo_fuente <- paste0("K-Means k = ", K_FUENTE)
}

cat(sprintf("Leyendo asignaciones K-Means desde hoja '%s'...\n", hoja_fuente))
asig_inicial <- read.xlsx(ruta_kmeans,
                          sheet    = hoja_fuente,
                          startRow = 2,
                          colNames = TRUE)
asig_inicial$Arbol   <- as.character(asig_inicial$Arbol)
asig_inicial$Cluster <- as.integer(asig_inicial$Cluster)

cat(sprintf("Árboles cargados: %d\n", nrow(asig_inicial)))
cat("Distribución inicial:\n")
print(sort(table(asig_inicial$Cluster), decreasing = TRUE))

# =============================================================================
# 2. FUNCIÓN DE SUBDIVISIÓN CON PAM SOBRE SUBMATRIZ RF
#    Devuelve: list(asig = vector nombrado, sil_avg = Silhouette promedio)
# =============================================================================
subdividir_con_pam <- function(ids_arboles, matriz_completa,
                               k = 2, seed = 2, min_n = 3) {

  n <- length(ids_arboles)

  if (n < max(k + 1, min_n)) {
    return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA,
                razon = "n_insuficiente"))
  }

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
# 3. BUCLE ITERATIVO
# =============================================================================
cat("\n=== INICIANDO SUBDIVISIÓN ITERATIVA (K-Means → PAM sobre RF) ===\n")
cat(sprintf("Umbral subdivisión : n > %d\n", TAMANO_MAX))
cat(sprintf("TAMANO_MIN_DURO   : n < %d no se subdividen\n", TAMANO_MIN_DURO))
cat(sprintf("DELTA_SIL         : %.3f (Sil mínima para aceptar corte)\n", DELTA_SIL))
cat(sprintf("k por corte       : %d\n\n", K_SUBDIVISION))

estado_actual       <- setNames(asig_inicial$Cluster, asig_inicial$Arbol)
log_subdivisiones   <- data.frame()
clusters_rechazados <- data.frame()
iteracion           <- 0L
tiempo_total        <- proc.time()

repeat {
  iteracion <- iteracion + 1L
  
  if (iteracion > MAX_ITERACIONES) {
    cat(sprintf("AVISO: límite de %d iteraciones alcanzado.\n", MAX_ITERACIONES))
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

  cat(sprintf("--- Iteración %d: %d cluster(s) candidatos ---\n",
              iteracion, length(clusters_grandes)))

  proximo_id   <- max(estado_actual) + 1L
  nuevo_estado <- estado_actual
  hubo_cambio  <- FALSE

  for (cid in clusters_grandes) {
    arboles_cid <- names(estado_actual[estado_actual == cid])
    n_cid       <- length(arboles_cid)

    cat(sprintf("\n  Cluster %d (n=%d):\n", cid, n_cid))

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
        Iteracion     = iteracion, Cluster = cid, N = n_cid,
        Sil_Propuesta = NA, Decision = resultado_sub$razon
      ))
      next
    }

    sil_propuesta <- resultado_sub$sil_avg

    if (sil_propuesta < DELTA_SIL) {
      cat(sprintf("    → Sil propuesta = %.4f < DELTA_SIL = %.3f. Conservando.\n",
                  sil_propuesta, DELTA_SIL))
      clusters_rechazados <- rbind(clusters_rechazados, data.frame(
        Iteracion     = iteracion, Cluster = cid, N = n_cid,
        Sil_Propuesta = sil_propuesta, Decision = "sil_insuficiente"
      ))
      next
    }

    cat(sprintf("    → Sil propuesta = %.4f >= %.3f. ACEPTADA.\n",
                sil_propuesta, DELTA_SIL))

    asig_nueva <- resultado_sub$asig
    sub_ids    <- sort(unique(asig_nueva))
    mapa_ids   <- setNames(
      c(cid, seq(proximo_id, proximo_id + length(sub_ids) - 2L)),
      sub_ids
    )
    asig_reetiquetada <- mapa_ids[as.character(asig_nueva)]

    ids_validos_cid               <- names(asig_nueva)
    nuevo_estado[ids_validos_cid] <- asig_reetiquetada
    proximo_id                    <- proximo_id + length(sub_ids) - 1L
    hubo_cambio                   <- TRUE

    tamanos_nuevos <- sort(as.integer(table(asig_reetiquetada)), decreasing = TRUE)
    log_subdivisiones <- rbind(log_subdivisiones, data.frame(
      Iteracion           = iteracion,
      Cluster_Original    = cid,
      N_Original          = n_cid,
      Sil_Sub             = round(sil_propuesta, 4),
      N_SubClusters       = length(tamanos_nuevos),
      Tamanos_SubClusters = paste(tamanos_nuevos, collapse = " | "),
      stringsAsFactors    = FALSE
    ))
  }

  estado_actual <- nuevo_estado

  if (!hubo_cambio) {
    cat(sprintf("Iteración %d: ningún cluster aceptó subdivisión (Sil < DELTA_SIL). Finalizando.\n",
                iteracion))
    break
  }
}

tiempo_total <- as.numeric(proc.time() - tiempo_total)[3]

# =============================================================================
# 4. TABLA FINAL
# =============================================================================
cat("\n=== CONSTRUYENDO TABLA FINAL ===\n")

ids_unicos  <- sort(unique(estado_actual))
mapa_final  <- setNames(seq_along(ids_unicos), ids_unicos)
cluster_fin <- mapa_final[as.character(estado_actual)]

asig_final_df <- data.frame(
  Arbol         = names(estado_actual),
  Cluster_Final = as.integer(cluster_fin),
  stringsAsFactors = FALSE
)
asig_final_df <- asig_final_df[order(asig_final_df$Cluster_Final), ]
rownames(asig_final_df) <- NULL

tamanos_finales <- as.data.frame(table(Cluster = asig_final_df$Cluster_Final))
tamanos_finales$Cluster <- as.integer(as.character(tamanos_finales$Cluster))
colnames(tamanos_finales)[2] <- "Tamano"
tamanos_finales <- tamanos_finales[order(tamanos_finales$Tamano, decreasing = TRUE), ]

n_clusters_final <- length(unique(asig_final_df$Cluster_Final))
n_bajo_15        <- sum(tamanos_finales$Tamano < 15)
n_bajo_50        <- sum(tamanos_finales$Tamano < 50)

cat(sprintf("Clusters iniciales (K-Means) : %d\n", length(unique(asig_inicial$Cluster))))
cat(sprintf("Clusters finales             : %d\n", n_clusters_final))
cat(sprintf("Clusters con n < 15         : %d (%.1f%%)\n",
            n_bajo_15, 100 * n_bajo_15 / n_clusters_final))
cat(sprintf("Clusters con n < 50         : %d (%.1f%%)\n",
            n_bajo_50, 100 * n_bajo_50 / n_clusters_final))
cat(sprintf("Árboles cubiertos            : %d / %d\n",
            nrow(asig_final_df), nrow(matriz_cuadrada)))
cat(sprintf("Iteraciones                  : %d\n", iteracion - 1L))
cat(sprintf("Tiempo total                 : %.1f s\n", tiempo_total))
cat("\nResumen tamaños:\n")
print(summary(tamanos_finales$Tamano))

# =============================================================================
# 4b. SILHOUETTE GLOBAL DE LA PARTICIÓN FINAL (sobre dist_rf)
# =============================================================================
cat("\n=== CALCULANDO SILHOUETTE GLOBAL DE LA PARTICIÓN FINAL ===\n")
cat("(sobre disimilitud RF — puede tardar varios minutos)\n")

ids_con_asig   <- asig_final_df$Arbol[asig_final_df$Arbol %in% rownames(matriz_cuadrada)]
submat_final   <- matriz_cuadrada[ids_con_asig, ids_con_asig]
dist_final     <- as.dist(submat_final)
asig_vec_final <- asig_final_df$Cluster_Final[match(ids_con_asig, asig_final_df$Arbol)]

tiempo_sil <- proc.time()
sil_final  <- tryCatch(
  silhouette(asig_vec_final, dist_final),
  error = function(e) { cat("ERROR silhouette global:", e$message, "\n"); NULL }
)
tiempo_sil <- as.numeric(proc.time() - tiempo_sil)[3]

sil_global_avg <- if (!is.null(sil_final)) round(mean(sil_final[, 3]), 4) else NA
cat(sprintf("Silhouette global (%d clusters, dist RF): %.4f\n",
            n_clusters_final, sil_global_avg))
cat(sprintf("Tiempo cálculo Silhouette: %.1f s\n", tiempo_sil))

if (!is.null(sil_final)) {
  sil_por_cluster <- tapply(sil_final[, 3], sil_final[, 1], mean)
  tamanos_finales$Silhouette_Cluster <- round(
    sil_por_cluster[as.character(tamanos_finales$Cluster)], 4
  )
} else {
  tamanos_finales$Silhouette_Cluster <- NA
}
tamanos_finales$Apto_Enriquecimiento <- tamanos_finales$Tamano >= 15

# =============================================================================
# 5. EXPORTAR
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

parametros_df <- data.frame(
  Parametro = c("algoritmo_origen", "algoritmo_subdivision", "k_fuente",
                "hoja_fuente", "tamano_max", "tamano_min_duro", "delta_sil",
                "k_subdivision", "max_iteraciones", "iteraciones_realizadas",
                "clusters_iniciales", "clusters_finales",
                "n_clusters_bajo_15", "n_clusters_bajo_50",
                "silhouette_global_final", "arboles_totales", "tiempo_total_s"),
  Valor = c("K-Means", "PAM sobre submatriz RF", as.character(K_FUENTE),
            hoja_fuente, TAMANO_MAX, TAMANO_MIN_DURO, DELTA_SIL,
            K_SUBDIVISION, MAX_ITERACIONES, iteracion - 1L,
            length(unique(asig_inicial$Cluster)), n_clusters_final,
            n_bajo_15, n_bajo_50,
            sil_global_avg,
            nrow(asig_final_df), round(tiempo_total, 1))
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb, "Asignaciones_K_Optimo",
                              paste0("Asignaciones Finales — ", n_clusters_final,
                                     " clusters (K-Means → PAM/RF desde ", titulo_fuente, ")"),
                              asig_final_df, anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Tamanos_y_Silhouette",
                              paste0("Tamaño y Silhouette por Cluster (Sil global = ",
                                     sil_global_avg, ")"),
                              tamanos_finales, anchos_col = "auto")

if (nrow(log_subdivisiones) > 0) {
  wb <- agregar_hoja_formateada(wb, "Log_Subdivisiones_Aceptadas",
                                "Log de Subdivisiones Aceptadas (Sil >= DELTA_SIL)",
                                log_subdivisiones, anchos_col = "auto")
}

if (nrow(clusters_rechazados) > 0) {
  wb <- agregar_hoja_formateada(wb, "Log_Rechazados",
                                "Clusters conservados (Sil < DELTA_SIL)",
                                clusters_rechazados, anchos_col = "auto")
}

wb <- agregar_hoja_formateada(wb, "Parametros",
                              "Parámetros del Proceso",
                              parametros_df, anchos_col = c(35, 25))

ruta_excel <- file.path(DIR_RESULTS, "kmeans_subdivision_iterativa.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Excel guardado:", ruta_excel, "\n")
cat("\n=== COMPLETADO ===\n")