# =============================================================================
# SUBDIVISIÓN ITERATIVA DE CLUSTERS GRANDES — CLARA
# Subdivide clusters con n > TAMANO_MAX hasta que todos tengan n <= TAMANO_MAX
# =============================================================================
library(cluster)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# PARÁMETROS AJUSTABLES
# =============================================================================
TAMANO_MAX      <- 200   # clusters con más de N árboles se subdividen
TAMANO_MIN_META <- 50    # tamaño mínimo deseado (informativo, no es un corte duro)
K_SUBDIVISION   <- 2     # k fijo para cada corte — divide en 2 en cada iteración
MAX_ITERACIONES <- 20    # tope de seguridad para evitar loops infinitos
MUESTRAS_CLARA  <- 50
SEED            <- 2

# k_extra disponibles en clara_resultados.xlsx (deben existir como hojas)
K_FUENTE        <- 10    # usar asignaciones de k=10 como punto de partida
# alternativa: "optimo" para usar k óptimo

# =============================================================================
# 1. CARGAR INSUMOS
# =============================================================================
cat("=== CARGANDO INSUMOS ===\n")

# Matriz RF
ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")
matriz_cuadrada <- if (file.exists(ruta_cache_matriz)) {
  readRDS(ruta_cache_matriz)
} else {
  as.matrix(read.table(file.path(DIR_RESULTS, "matriz_rf_conjunto.csv"),
                       sep = ";", header = TRUE, row.names = 1,
                       check.names = FALSE))
}
cat("Matriz RF:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# Asignaciones CLARA iniciales
ruta_clara <- file.path(DIR_RESULTS, "clara_resultados.xlsx")

if (K_FUENTE == "optimo") {
  hoja_fuente  <- "Asignaciones_K_Optimo"
  titulo_fuente <- "k óptimo"
} else {
  hoja_fuente  <- paste0("Asignaciones_K_", K_FUENTE)
  titulo_fuente <- paste0("k = ", K_FUENTE)
}

cat(sprintf("Leyendo asignaciones iniciales desde hoja '%s'...\n", hoja_fuente))

asig_inicial <- read.xlsx(ruta_clara,
                          sheet    = hoja_fuente,
                          startRow = 2,
                          colNames = TRUE)

asig_inicial$Arbol   <- as.character(asig_inicial$Arbol)
asig_inicial$Cluster <- as.integer(asig_inicial$Cluster)

cat(sprintf("Árboles cargados: %d\n", nrow(asig_inicial)))
cat("Distribución inicial:\n")
print(sort(table(asig_inicial$Cluster), decreasing = TRUE))

# =============================================================================
# 2. FUNCIÓN DE SUBDIVISIÓN DE UN CLUSTER
# =============================================================================
subdividir_con_clara <- function(ids_arboles, submatriz, k, muestras, seed) {
  
  n          <- length(ids_arboles)
  k_real     <- min(k, n - 1)
  
  if (k_real < 2) {
    cat(sprintf("    n=%d demasiado pequeño para subdividir. Conservando.\n", n))
    return(setNames(rep(1L, n), ids_arboles))
  }
  
  set.seed(seed)
  
  resultado <- tryCatch(
    clara(submatriz, k = k_real, metric = "euclidean",
          samples   = muestras,
          sampsize  = min(n, 40 + 2 * k_real),
          keep.data = FALSE,
          rngR      = TRUE),
    error = function(e) {
      cat(sprintf("    ERROR clara: %s\n", e$message))
      NULL
    }
  )
  
  if (is.null(resultado)) return(setNames(rep(1L, n), ids_arboles))
  
  setNames(as.integer(resultado$clustering), ids_arboles)
}

# =============================================================================
# 3. BUCLE ITERATIVO DE SUBDIVISIÓN
# =============================================================================
cat("\n=== INICIANDO SUBDIVISIÓN ITERATIVA ===\n")
cat(sprintf("Umbral de subdivisión : n > %d\n", TAMANO_MAX))
cat(sprintf("Objetivo              : clusters entre %d y %d árboles\n",
            TAMANO_MIN_META, TAMANO_MAX))
cat(sprintf("k por corte           : %d\n\n", K_SUBDIVISION))

# Estado inicial: vector nombrado Arbol -> Cluster_Final
estado_actual <- setNames(asig_inicial$Cluster, asig_inicial$Arbol)

# Log de todas las subdivisiones realizadas
log_subdivisiones <- data.frame()
iteracion         <- 0L
tiempo_total      <- proc.time()

repeat {
  iteracion <- iteracion + 1L
  
  if (iteracion > MAX_ITERACIONES) {
    cat(sprintf("AVISO: Se alcanzó el límite de %d iteraciones.\n", MAX_ITERACIONES))
    break
  }
  
  tamanos_actuales  <- table(estado_actual)
  clusters_grandes  <- as.integer(names(tamanos_actuales[tamanos_actuales > TAMANO_MAX]))
  
  if (length(clusters_grandes) == 0) {
    cat(sprintf("Iteración %d: todos los clusters tienen n <= %d. Finalizado.\n",
                iteracion, TAMANO_MAX))
    break
  }
  
  cat(sprintf("--- Iteración %d: %d cluster(s) a subdividir ---\n",
              iteracion, length(clusters_grandes)))
  
  # Reetiquetado global: el próximo cluster nuevo recibe el id max+1
  proximo_id <- max(estado_actual) + 1L
  nuevo_estado <- estado_actual
  
  for (cid in clusters_grandes) {
    
    arboles_cid     <- names(estado_actual[estado_actual == cid])
    n_cid           <- length(arboles_cid)
    arboles_validos <- arboles_cid[arboles_cid %in% rownames(matriz_cuadrada)]
    
    cat(sprintf("  Cluster %d (n=%d)...\n", cid, n_cid))
    
    submatriz   <- matriz_cuadrada[arboles_validos, arboles_validos]
    asig_nueva  <- subdividir_con_clara(arboles_validos, submatriz,
                                        k        = K_SUBDIVISION,
                                        muestras = MUESTRAS_CLARA,
                                        seed     = SEED)
    
    # Sub-cluster 1 conserva el id original; sub-cluster 2+ recibe ids nuevos
    sub_ids         <- sort(unique(asig_nueva))
    mapa_ids        <- setNames(c(cid, seq(proximo_id, proximo_id + length(sub_ids) - 2L)),
                                sub_ids)
    asig_reetiquetada <- mapa_ids[as.character(asig_nueva)]
    
    tamanos_nuevos <- table(asig_reetiquetada)
    cat(sprintf("    → %d sub-clusters: %s\n",
                length(tamanos_nuevos),
                paste(paste0("n=", sort(as.integer(tamanos_nuevos),
                                        decreasing = TRUE)),
                      collapse = ", ")))
    
    nuevo_estado[arboles_validos] <- asig_reetiquetada
    proximo_id <- proximo_id + length(sub_ids) - 1L
    
    # Registrar en log
    log_subdivisiones <- rbind(log_subdivisiones, data.frame(
      Iteracion          = iteracion,
      Cluster_Original   = cid,
      N_Original         = n_cid,
      N_SubClusters      = length(tamanos_nuevos),
      Tamanos_SubClusters = paste(sort(as.integer(tamanos_nuevos),
                                       decreasing = TRUE), collapse = " | ")
    ))
  }
  
  estado_actual <- nuevo_estado
}

tiempo_total <- as.numeric(proc.time() - tiempo_total)[3]

# =============================================================================
# 4. CONSTRUIR TABLA FINAL DE ASIGNACIONES
# =============================================================================
cat("\n=== CONSTRUYENDO TABLA FINAL ===\n")

# Reetiquetado secuencial limpio (1, 2, 3, ...) para el resultado final
ids_originales <- sort(unique(estado_actual))
mapa_final     <- setNames(seq_along(ids_originales), ids_originales)
cluster_final  <- mapa_final[as.character(estado_actual)]

asig_final_df <- data.frame(
  Arbol          = names(estado_actual),
  Cluster_Final  = as.integer(cluster_final),
  stringsAsFactors = FALSE
)
asig_final_df <- asig_final_df[order(asig_final_df$Cluster_Final), ]
rownames(asig_final_df) <- NULL

tamanos_finales <- as.data.frame(table(Cluster = asig_final_df$Cluster_Final))
tamanos_finales$Cluster <- as.integer(as.character(tamanos_finales$Cluster))
colnames(tamanos_finales)[2] <- "Tamano"
tamanos_finales <- tamanos_finales[order(tamanos_finales$Tamano, decreasing = TRUE), ]

cat(sprintf("Clusters iniciales   : %d\n", length(unique(asig_inicial$Cluster))))
cat(sprintf("Clusters finales     : %d\n", length(unique(asig_final_df$Cluster_Final))))
cat(sprintf("Árboles cubiertos    : %d de %d\n",
            nrow(asig_final_df), nrow(matriz_cuadrada)))
cat(sprintf("Iteraciones realizadas: %d\n", iteracion - 1L))
cat(sprintf("Tiempo total         : %.1f segundos\n", tiempo_total))
cat("\nResumen de tamaños finales:\n")
print(summary(tamanos_finales$Tamano))
cat(sprintf("Clusters dentro del rango [%d, %d]: %d de %d (%.1f%%)\n",
            TAMANO_MIN_META, TAMANO_MAX,
            sum(tamanos_finales$Tamano >= TAMANO_MIN_META &
                  tamanos_finales$Tamano <= TAMANO_MAX),
            nrow(tamanos_finales),
            100 * mean(tamanos_finales$Tamano >= TAMANO_MIN_META &
                         tamanos_finales$Tamano <= TAMANO_MAX)))

# =============================================================================
# 5. EXPORTAR
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

# Parámetros del proceso
parametros_df <- data.frame(
  Parametro = c("k_fuente", "hoja_fuente", "tamano_max_subdivision",
                "tamano_min_meta", "k_por_corte", "max_iteraciones",
                "iteraciones_realizadas", "muestras_clara",
                "clusters_iniciales", "clusters_finales",
                "arboles_totales", "tiempo_total_s"),
  Valor     = c(as.character(K_FUENTE), hoja_fuente, TAMANO_MAX,
                TAMANO_MIN_META, K_SUBDIVISION, MAX_ITERACIONES,
                iteracion - 1L, MUESTRAS_CLARA,
                length(unique(asig_inicial$Cluster)),
                length(unique(asig_final_df$Cluster_Final)),
                nrow(asig_final_df), round(tiempo_total, 1))
)

tiempos_df <- data.frame(
  Proceso  = c("Subdivisión iterativa CLARA", "TOTAL"),
  Tiempo_s = c(round(tiempo_total, 2), round(tiempo_total, 2))
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb, "Asignaciones_K_Optimo",
                              paste0("Asignaciones Finales — ",
                                     length(unique(asig_final_df$Cluster_Final)),
                                     " clusters (subdivisión iterativa desde ",
                                     titulo_fuente, ")"),
                              asig_final_df,
                              anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Tamanos_Finales",
                              paste0("Tamaño de Clusters Finales (n <= ", TAMANO_MAX, ")"),
                              tamanos_finales,
                              anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Log_Subdivisiones",
                              "Log de Subdivisiones por Iteración",
                              log_subdivisiones,
                              anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Parametros",
                              "Parámetros del Proceso",
                              parametros_df,
                              anchos_col = c(35, 25))

wb <- agregar_hoja_formateada(wb, "Tiempos",
                              "Tiempos de Ejecución",
                              tiempos_df,
                              anchos_col = "auto")

ruta_excel <- file.path(DIR_RESULTS, "clara_subdivision_iterativa.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Excel guardado:", ruta_excel, "\n")
cat("\n=== COMPLETADO ===\n")