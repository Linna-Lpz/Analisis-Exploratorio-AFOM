# =============================================================================
# SUBDIVISIÓN ITERATIVA DE CLUSTERS GRANDES — K-MEANS
# Lee asignaciones de K-Means y subdivide con K-Means los clusters con
# n > TAMANO_MAX hasta que todos tengan n <= TAMANO_MAX
# =============================================================================
library(cluster)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# PARÁMETROS AJUSTABLES
# =============================================================================
TAMANO_MAX      <- 200   # clusters con más de N árboles se subdividen
TAMANO_MIN_META <- 50    # tamaño mínimo deseado (informativo)
TAM_MIN_CORTE   <- 20    # sub-cluster con MENOS de este nº invalida el corte
K_INICIAL       <- 2     # k con el que se intenta primero cada corte
K_MAXIMO        <- 6     # k máximo a probar si los intentos con k menor fallan
NSTART          <- 25    # nstart de kmeans (igual que en el script principal)
MAX_ITERACIONES <- 20    # tope de seguridad contra loops infinitos
SEEDS           <- c(2L, 13L, 42L, 77L, 123L)

# Fuente de asignaciones K-Means: número (10, 15) o "optimo"
K_FUENTE <- 15

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
# 2. FUNCIÓN DE SUBDIVISIÓN CON K-MEANS
#    Prueba k = K_INICIAL..K_MAXIMO × SEEDS.
#    Acepta el primer corte donde todos los sub-clusters >= TAM_MIN_CORTE.
#    Devuelve NULL si no encuentra ningún corte válido.
# =============================================================================
subdividir <- function(ids_arboles, submatriz, k_inicial, k_maximo,
                       tam_min_corte, nstart, seeds) {
  n <- length(ids_arboles)
  
  for (k in seq(k_inicial, k_maximo)) {
    if (k >= n) break
    
    for (s in seeds) {
      set.seed(s)
      res <- tryCatch(
        kmeans(submatriz, centers = k, nstart = nstart),
        error = function(e) NULL
      )
      if (is.null(res)) next
      
      tam_sub <- as.integer(table(res$cluster))
      
      if (all(tam_sub >= tam_min_corte)) {
        cat(sprintf("    Corte aceptado: k=%d, seed=%d → sub-clusters: %s\n",
                    k, s,
                    paste(sort(tam_sub, decreasing = TRUE), collapse = " | ")))
        return(setNames(as.integer(res$cluster), ids_arboles))
      } else {
        cat(sprintf("    k=%d, seed=%d rechazado — mínimo sub-cluster: %d\n",
                    k, s, min(tam_sub)))
      }
    }
  }
  
  return(NULL)
}

# =============================================================================
# 3. BUCLE ITERATIVO
# =============================================================================
cat("\n=== INICIANDO SUBDIVISIÓN ITERATIVA (K-Means) ===\n")
cat(sprintf("Umbral subdivisión     : n > %d\n", TAMANO_MAX))
cat(sprintf("Mínimo por sub-cluster : %d\n", TAM_MIN_CORTE))
cat(sprintf("k a probar             : %d .. %d\n", K_INICIAL, K_MAXIMO))
cat(sprintf("nstart                 : %d\n", NSTART))
cat(sprintf("Seeds por intento      : %s\n\n", paste(SEEDS, collapse = ", ")))

estado_actual         <- setNames(asig_inicial$Cluster, asig_inicial$Arbol)
log_subdivisiones     <- data.frame()
clusters_irresolubles <- integer(0)
iteracion             <- 0L
tiempo_total          <- proc.time()

repeat {
  iteracion <- iteracion + 1L
  
  if (iteracion > MAX_ITERACIONES) {
    cat(sprintf("AVISO: límite de %d iteraciones alcanzado.\n", MAX_ITERACIONES))
    break
  }
  
  tamanos_actuales <- table(estado_actual)
  clusters_grandes <- setdiff(
    as.integer(names(tamanos_actuales[tamanos_actuales > TAMANO_MAX])),
    clusters_irresolubles
  )
  
  if (length(clusters_grandes) == 0) {
    cat(sprintf("Iteración %d: no quedan clusters subdividibles. Finalizado.\n", iteracion))
    break
  }
  
  cat(sprintf("--- Iteración %d: %d cluster(s) a subdividir ---\n",
              iteracion, length(clusters_grandes)))
  
  proximo_id   <- max(estado_actual) + 1L
  nuevo_estado <- estado_actual
  
  for (cid in clusters_grandes) {
    arboles_cid     <- names(estado_actual[estado_actual == cid])
    n_cid           <- length(arboles_cid)
    arboles_validos <- arboles_cid[arboles_cid %in% rownames(matriz_cuadrada)]
    
    cat(sprintf("\n  Cluster %d (n=%d)...\n", cid, n_cid))
    
    submatriz  <- matriz_cuadrada[arboles_validos, arboles_validos]
    asig_nueva <- subdividir(
      ids_arboles   = arboles_validos,
      submatriz     = submatriz,
      k_inicial     = K_INICIAL,
      k_maximo      = K_MAXIMO,
      tam_min_corte = TAM_MIN_CORTE,
      nstart        = NSTART,
      seeds         = SEEDS
    )
    
    if (is.null(asig_nueva)) {
      cat(sprintf("    Cluster %d marcado como irresoluble (conservado intacto).\n", cid))
      clusters_irresolubles <- c(clusters_irresolubles, cid)
      
      log_subdivisiones <- rbind(log_subdivisiones, data.frame(
        Iteracion           = iteracion,
        Cluster_Original    = cid,
        N_Original          = n_cid,
        N_SubClusters       = 0L,
        Tamanos_SubClusters = "NO_SUBDIVIDIDO",
        stringsAsFactors    = FALSE
      ))
      next
    }
    
    sub_ids  <- sort(unique(asig_nueva))
    mapa_ids <- setNames(
      c(cid, seq(proximo_id, proximo_id + length(sub_ids) - 2L)),
      sub_ids
    )
    asig_reetiquetada             <- mapa_ids[as.character(asig_nueva)]
    nuevo_estado[arboles_validos] <- asig_reetiquetada
    proximo_id                    <- proximo_id + length(sub_ids) - 1L
    
    tamanos_nuevos <- sort(as.integer(table(asig_reetiquetada)), decreasing = TRUE)
    log_subdivisiones <- rbind(log_subdivisiones, data.frame(
      Iteracion           = iteracion,
      Cluster_Original    = cid,
      N_Original          = n_cid,
      N_SubClusters       = length(tamanos_nuevos),
      Tamanos_SubClusters = paste(tamanos_nuevos, collapse = " | "),
      stringsAsFactors    = FALSE
    ))
  }
  
  estado_actual <- nuevo_estado
  cat("\n")
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

cat(sprintf("Clusters iniciales (K-Means) : %d\n", length(unique(asig_inicial$Cluster))))
cat(sprintf("Clusters finales             : %d\n", nrow(tamanos_finales)))
cat(sprintf("  — dentro [%d, %d]         : %d\n",
            TAMANO_MIN_META, TAMANO_MAX,
            sum(tamanos_finales$Tamano >= TAMANO_MIN_META &
                  tamanos_finales$Tamano <= TAMANO_MAX)))
cat(sprintf("  — irresolubles             : %d\n", length(clusters_irresolubles)))
cat(sprintf("Árboles cubiertos            : %d / %d\n",
            nrow(asig_final_df), nrow(matriz_cuadrada)))
cat(sprintf("Iteraciones                  : %d\n", iteracion - 1L))
cat(sprintf("Tiempo total                 : %.1f s\n", tiempo_total))
cat("\nResumen tamaños:\n")
print(summary(tamanos_finales$Tamano))

if (length(clusters_irresolubles) > 0) {
  ids_irr <- mapa_final[as.character(clusters_irresolubles)]
  ids_irr <- ids_irr[!is.na(ids_irr)]
  cat(sprintf("\nAVISO: %d cluster(s) no subdivididos (> %d árboles):\n",
              length(ids_irr), TAMANO_MAX))
  print(tamanos_finales[tamanos_finales$Cluster %in% ids_irr, ])
}

# =============================================================================
# 5. EXPORTAR
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

parametros_df <- data.frame(
  Parametro = c("algoritmo_origen", "algoritmo_subdivision", "k_fuente",
                "hoja_fuente", "tamano_max_subdivision", "tamano_min_meta",
                "tam_min_corte", "k_inicial", "k_maximo", "nstart",
                "seeds", "max_iteraciones", "iteraciones_realizadas",
                "clusters_iniciales", "clusters_finales",
                "clusters_irresolubles", "arboles_totales", "tiempo_total_s"),
  Valor = c("K-Means", "K-Means", as.character(K_FUENTE),
            hoja_fuente, TAMANO_MAX, TAMANO_MIN_META, TAM_MIN_CORTE,
            K_INICIAL, K_MAXIMO, NSTART,
            paste(SEEDS, collapse = ","), MAX_ITERACIONES, iteracion - 1L,
            length(unique(asig_inicial$Cluster)), nrow(tamanos_finales),
            length(clusters_irresolubles),
            nrow(asig_final_df), round(tiempo_total, 1))
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb, "Asignaciones_K_Optimo",
                              paste0("Asignaciones Finales — ", nrow(tamanos_finales),
                                     " clusters (subdivisión iterativa K-Means desde ", titulo_fuente, ")"),
                              asig_final_df, anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Tamanos_Finales",
                              paste0("Tamaño de Clusters Finales (objetivo <= ", TAMANO_MAX, ")"),
                              tamanos_finales, anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Log_Subdivisiones",
                              "Log de Subdivisiones por Iteración",
                              log_subdivisiones, anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Parametros",
                              "Parámetros del Proceso",
                              parametros_df, anchos_col = c(35, 25))

ruta_excel <- file.path(DIR_RESULTS, "kmeans_subdivision_iterativa.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Excel guardado:", ruta_excel, "\n")
cat("\n=== COMPLETADO ===\n")