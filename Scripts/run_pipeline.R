# ==============================================================================
# run_pipeline.R
# PIPELINE MAESTRO DINÁMICO — Ejecuta los scripts en orden secuencial
# ==============================================================================
library(here)
library(openxlsx)

source(here("Scripts", "config.R"))
source(here("Scripts", "00_funciones_globales.R"))

# ==============================================================================
# 1. CONFIGURACIÓN DEL PIPELINE
# ==============================================================================

# Variable de Control Dinámico:
# "AUTO" : El pipeline evaluará empíricamente los algoritmos y orquestará 
#          los pasos siguientes (downstream) para el ganador (mejor Silhouette).
# Alternativamente, puedes forzar el algoritmo de downstream usando uno de los siguientes:
# "MST-kNN", "CLARA", "K-Means", "PAM"
ALGORITMO_DOWNSTREAM <- "AUTO"

# -- Pasos a ejecutar (TRUE = ejecutar, FALSE = saltar) --
PASOS_EJECUTAR <- list(
  # PASOS BASE
  paso_01  = TRUE,  # 01 — Análisis por árbol
  paso_02  = TRUE,  # 02 — Árbol medioide
  paso_03  = TRUE,  # 03 — Comparar medioide vs árboles
  paso_04  = TRUE,  # 04 — Matriz Robinson-Foulds
  paso_05a = TRUE,  # 05a — Clustering K-Means
  paso_05b = TRUE,  # 05b — Clustering PAM
  paso_05c = TRUE,  # 05c — Clustering CLARA
  paso_05d = TRUE,  # 05d — Clustering MST-kNN
  paso_05e = TRUE,  # 05e — Comparativa Silhouette
  
  # PASOS DOWNSTREAM (Dependen del algoritmo seleccionado o AUTO)
  paso_07  = TRUE,  # 07 — Iterar clustering
  paso_06a = TRUE,  # 06a — UMAP (Reducc. Dim.)
  paso_06b = TRUE,  # 06b — t-SNE (Reducc. Dim.)
  paso_06c = TRUE,  # 06c — PCA (Reducc. Dim.)
  paso_08  = TRUE,  # 08 — Etiquetas HGNC
  paso_09  = TRUE,  # 09 — Enriquecimiento
  paso_10  = TRUE,  # 10 — Burbujas Enriquecimiento
  paso_11  = TRUE,  # 11 — Evaluación Sesgo
  paso_12  = TRUE   # 12 — Validación Nula (GO)
)

# ¿Detener el pipeline si un paso falla?
DETENER_EN_ERROR <- TRUE

# ==============================================================================
# 2. DEFINICIÓN DE PASOS BASE (Siempre se ejecutan)
# ==============================================================================
# Estos pasos construyen el espacio filogenético, evalúan todos los candidatos
# y realizan la reducción dimensional inicial.
PASOS_BASE <- list(
  paso_01  = list("01 — Análisis por árbol",           "Scripts/01_analisis_por_arbol.R"),
  paso_02  = list("02 — Árbol medioide",                "Scripts/02_generar_arbol_medioide.R"),
  paso_03  = list("03 — Comparar medioide vs árboles", "Scripts/03_comparar_medioide_vs_arboles.R"),
  paso_04  = list("04 — Matriz Robinson-Foulds",       "Scripts/04_calcular_matriz_rf.R"),
  paso_05a = list("05a — Clustering K-Means",          "Scripts/06a_clustering_kmeans.R"),
  paso_05b = list("05b — Clustering PAM",              "Scripts/06b_clustering_pam.R"),
  paso_05c = list("05c — Clustering CLARA",            "Scripts/06c_clustering_clara.R"),
  paso_05d = list("05d — Clustering MST-kNN",          "Scripts/06d_clustering_mstknn.R"),
  paso_05e = list("05e — Comparativa Silhouette",      "Scripts/06e_comparativa_silhouette_rf.R")
)

# ==============================================================================
# 3. INFRAESTRUCTURA INTERNA
# ==============================================================================

log_pipeline <- data.frame(
  Paso    = character(),
  Script  = character(),
  Estado  = character(),
  Tiempo  = numeric(),
  Mensaje = character(),
  stringsAsFactors = FALSE
)

log_memoria <- NULL

ejecutar_paso <- function(nombre_paso, ruta_script, clave_paso = "Paso") {
  cat(sprintf("\n%s\n>>> INICIANDO: %s\n%s\n", strrep("=", 70), nombre_paso, strrep("=", 70)))
  
  if (!file.exists(ruta_script)) {
    msg <- paste("Script no encontrado:", ruta_script)
    cat("  [OMITIDO]", msg, "\n")
    return(list(estado = "OMITIDO", tiempo = 0, mensaje = msg))
  }
  
  t_inicio <- proc.time()
  
  resultado <- tryCatch({
    env_aislado <- new.env(parent = globalenv())
    source(ruta_script, echo = FALSE, local = env_aislado)
    rm(env_aislado)
    list(estado = "OK", mensaje = "")
  }, error = function(e) {
    list(estado = "ERROR", mensaje = conditionMessage(e))
  })
  
  t_fin     <- proc.time()
  t_elapsed <- round((t_fin - t_inicio)[["elapsed"]], 1)
  
  if (resultado$estado == "OK") {
    cat(sprintf("\n  [OK] %s completado en %.1f s\n", nombre_paso, t_elapsed))
  } else {
    cat(sprintf("\n  [ERROR] %s falló después de %.1f s\n", nombre_paso, t_elapsed))
    cat("  Mensaje:", resultado$mensaje, "\n")
  }
  
  return(list(estado = resultado$estado, tiempo = t_elapsed, mensaje = resultado$mensaje, clave = clave_paso, script = ruta_script))
}

registrar_ejecucion <- function(log_df, res) {
  rbind(log_df, data.frame(
    Paso    = res$clave,
    Script  = basename(res$script),
    Estado  = res$estado,
    Tiempo  = res$tiempo,
    Mensaje = res$mensaje,
    stringsAsFactors = FALSE
  ))
}

# ==============================================================================
# 4. FASE 1: EJECUCIÓN BASE Y EVALUACIÓN
# ==============================================================================

t_total_inicio <- proc.time()
cat("\n", strrep("*", 70), "\n")
cat("   PIPELINE DINÁMICO AFOM — OrthoMaM\n")
cat("   Inicio:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("*", 70), "\n")

pipeline_interrumpido <- FALSE

for (clave in names(PASOS_BASE)) {
  nombre <- PASOS_BASE[[clave]][[1]]
  script <- here(PASOS_BASE[[clave]][[2]])
  
  # Verificar si se debe ejecutar el paso
  if (!is.null(PASOS_EJECUTAR[[clave]]) && !PASOS_EJECUTAR[[clave]]) {
    cat(sprintf("\n>>> SALTANDO: %s (Desactivado en configuración)\n", nombre))
    log_pipeline <- registrar_ejecucion(log_pipeline, list(
      clave = clave, script = script, estado = "SALTADO", tiempo = 0, mensaje = "Saltado por configuración"
    ))
    next
  }

  log_memoria <- registrar_memoria(log_memoria, nombre, "ANTES")
  res <- ejecutar_paso(nombre, script, clave_paso = clave)
  limpiar_memoria()
  log_memoria <- registrar_memoria(log_memoria, nombre, "DESPUÉS")
  
  log_pipeline <- registrar_ejecucion(log_pipeline, res)
  
  if (res$estado == "ERROR" && DETENER_EN_ERROR) {
    cat("\n", strrep("!", 70), "\n  PIPELINE DETENIDO en:", nombre, "\n", strrep("!", 70), "\n")
    pipeline_interrumpido <- TRUE
    break
  }
}

# ==============================================================================
# 5. FASE 2: DECISIÓN DINÁMICA DE DOWNSTREAM
# ==============================================================================

if (!pipeline_interrumpido) {
  cat("\n", strrep("=", 70), "\n")
  cat(">>> FASE DE DECISIÓN DINÁMICA\n")
  cat(strrep("=", 70), "\n")
  
  if (toupper(ALGORITMO_DOWNSTREAM) == "AUTO") {
    cat("Evaluando empíricamente el mejor algoritmo desde caché (06e)...\n")
    ruta_rds <- here("Resultados", "cache", "comparativa_silhouette_rf.rds")
    
    if (file.exists(ruta_rds)) {
      tabla_comp <- readRDS(ruta_rds)
      mejor_alg <- tabla_comp$Algoritmo[which.max(tabla_comp$Silhouette)]
      ALGORITMO_DOWNSTREAM <- mejor_alg
      cat(sprintf("=> GANADOR AUTOMÁTICO: %s\n", mejor_alg))
    } else {
      stop("[ERROR] No se encontró el archivo 'comparativa_silhouette_rf.rds'. Ejecute los pasos base primero.")
    }
  } else {
    cat(sprintf("=> ALGORITMO FORZADO POR USUARIO: %s\n", ALGORITMO_DOWNSTREAM))
  }
  
  # Seleccionar los scripts de downstream correspondientes al ganador
  # El script iterar se infiere con un switch para mantener el flujo unificado
  script_iterar <- switch(ALGORITMO_DOWNSTREAM,
                          "CLARA" = "Scripts/07_iterar_clara.R",
                          "K-Means" = "Scripts/07_iterar_k-means.R",
                          "PAM" = "Scripts/07_iterar_pam.R",
                          "MST-kNN" = "Scripts/07_iterar_mstknn.R",
                          "Scripts/07_iterar.R") # fallback por si acaso
  
  PASOS_DOWNSTREAM <- list(
    paso_07 = list(paste("07 — Iterar", ALGORITMO_DOWNSTREAM),   script_iterar),
    paso_06a = list("06a — UMAP (Reducc. Dim.)",         "Scripts/05a_reduccion_dimensional_umap.R"),
    paso_06b = list("06b — t-SNE (Reducc. Dim.)",        "Scripts/05c_reduccion_dimensional_tSNE.R"),
    paso_06c = list("06c — PCA (Reducc. Dim.)",          "Scripts/05d_reduccion_dimensional_PCA.R"),
    paso_08 = list(paste("08 — Etiquetas HGNC", ALGORITMO_DOWNSTREAM), "Scripts/08_etiquetas.R"),
    paso_09 = list(paste("09 — Enriquecimiento", ALGORITMO_DOWNSTREAM), "Scripts/09_enriquecimiento.R"),
    paso_10 = list(paste("10 — Burbujas Enriquecimiento", ALGORITMO_DOWNSTREAM), "Scripts/10_grafico_enriquecimiento.R")
  )
  
  # Añadir Evaluación de Sesgo y Modelo Nulo Global si es que hay un downstream válido
  if (length(PASOS_DOWNSTREAM) > 0) {
    PASOS_DOWNSTREAM$paso_11 <- list("11 — Evaluación Sesgo", "Scripts/evaluacion_sesgo.R")
    PASOS_DOWNSTREAM$paso_12 <- list("12 — Validación Nula (GO)", "Scripts/12_validacion_nula_GO.R")
  }
  
  # Ejecutar Downstream
  for (clave in names(PASOS_DOWNSTREAM)) {
    nombre <- PASOS_DOWNSTREAM[[clave]][[1]]
    script <- here(PASOS_DOWNSTREAM[[clave]][[2]])
    
    # Verificar si se debe ejecutar el paso
    if (!is.null(PASOS_EJECUTAR[[clave]]) && !PASOS_EJECUTAR[[clave]]) {
      cat(sprintf("\n>>> SALTANDO: %s (Desactivado en configuración)\n", nombre))
      log_pipeline <- registrar_ejecucion(log_pipeline, list(
        clave = clave, script = script, estado = "SALTADO", tiempo = 0, mensaje = "Saltado por configuración"
      ))
      next
    }

    log_memoria <- registrar_memoria(log_memoria, nombre, "ANTES")
    res <- ejecutar_paso(nombre, script, clave_paso = clave)
    limpiar_memoria()
    log_memoria <- registrar_memoria(log_memoria, nombre, "DESPUÉS")
    
    log_pipeline <- registrar_ejecucion(log_pipeline, res)
    
    if (res$estado == "ERROR" && DETENER_EN_ERROR) {
      cat("\n", strrep("!", 70), "\n  PIPELINE DETENIDO en:", nombre, "\n", strrep("!", 70), "\n")
      pipeline_interrumpido <- TRUE
      break
    }
  }
}

# ==============================================================================
# 6. INFORME FINAL Y EXPORTACIÓN
# ==============================================================================

t_total <- round((proc.time() - t_total_inicio)[["elapsed"]], 1)

cat("\n", strrep("*", 70), "\n")
cat("   RESUMEN DEL PIPELINE\n")
cat("   Fin:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("   Tiempo total: %.1f s  (%.1f min)\n", t_total, t_total / 60))
cat(strrep("*", 70), "\n\n")

cat(sprintf("  %-8s  %-45s  %-8s  %s\n", "Paso", "Script", "Estado", "Tiempo(s)"))
cat("  ", strrep("-", 70), "\n")
for (i in seq_len(nrow(log_pipeline))) {
  fila  <- log_pipeline[i, ]
  icono <- switch(fila$Estado, "OK" = "[OK]   ", "ERROR" = "[ERROR]", "OMITIDO" = "[miss] ", "SALTADO" = "[skip] ", "[?]    ")
  cat(sprintf("  %s  %-45s  %-8s  %.1f\n", icono, fila$Script, fila$Estado, fila$Tiempo))
}

errores <- log_pipeline[log_pipeline$Estado == "ERROR", ]
if (nrow(errores) > 0) {
  cat("\n  ERRORES ENCONTRADOS:\n")
  for (i in seq_len(nrow(errores))) {
    cat(sprintf("  - %s:\n    %s\n", errores$Script[i], errores$Mensaje[i]))
  }
}

if (!pipeline_interrumpido) {
  cat("\n  Pipeline completado exitosamente!\n")
}

# Exportar Logs
timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
ruta_log <- here("Resultados", paste0("pipeline_log_", timestamp_str, ".csv"))
write.csv(log_pipeline, ruta_log, row.names = FALSE)
cat("\n  Log general guardado en:", ruta_log, "\n")

if (!is.null(log_memoria) && nrow(log_memoria) > 0) {
  ruta_mem_csv <- here("Resultados", paste0("pipeline_memoria_", timestamp_str, ".csv"))
  write.csv(log_memoria, ruta_mem_csv, row.names = FALSE)
  
  ruta_mem_xlsx <- here("Resultados", paste0("pipeline_memoria_", timestamp_str, ".xlsx"))
  tryCatch({
    wb_mem <- createWorkbook()
    wb_mem <- agregar_hoja_formateada(wb_mem, "Log_Memoria", "Monitoreo RAM", log_memoria, c(22, 35, 12, 14, 16, 16, 10, 50))
    wb_mem <- agregar_hoja_formateada(wb_mem, "Resumen_Pipeline", "Resumen", log_pipeline, c(12, 45, 10, 12, 50))
    saveWorkbook(wb_mem, ruta_mem_xlsx, overwrite = TRUE)
    cat("  Log memoria (Excel) guardado en:", ruta_mem_xlsx, "\n")
  }, error = function(e) cat("  Aviso: no se pudo guardar memoria Excel:", e$message, "\n"))
  
  cat("\n", strrep("-", 70), "\n   RESUMEN DE MEMORIA\n", strrep("-", 70), "\n")
  cat(sprintf("  Pico RAM R     : %.1f MB\n", max(log_memoria$RAM_R_MB, na.rm = TRUE)))
  if (any(!is.na(log_memoria$RAM_Sistema_MB))) {
    cat(sprintf("  Pico RAM Sistema: %.1f MB (%.1f%%)\n", max(log_memoria$RAM_Sistema_MB, na.rm = TRUE), max(log_memoria$Pct_Uso, na.rm = TRUE)))
  }
  cat(strrep("-", 70), "\n\n")
}
