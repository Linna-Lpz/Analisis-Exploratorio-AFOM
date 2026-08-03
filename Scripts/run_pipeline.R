# ==============================================================================
# run_pipeline.R
# PIPELINE MAESTRO — Ejecuta todos los scripts en orden secuencial
#
# USO:
#   Opción A (RStudio Server): Abrir este archivo y pulsar "Source"
#   Opción B (Terminal/Rscript): Rscript Scripts/run_pipeline.R
#
# REQUISITOS PREVIOS:
#   - El proyecto debe abrirse desde el .Rproj (para que here() funcione)
#   - Los datos originales deben estar en Datos/Original/archiveTreesV12/
#   - Los paquetes deben estar instalados (ver Scripts/instalar_paquetes.R)
#
# CONTROL DE EJECUCIÓN:
#   - Cambia a FALSE los pasos que NO quieras re-ejecutar (ej. si ya tienes caché)
#   - Si un paso falla, el pipeline detiene la ejecución e informa el error
#
# GESTIÓN DE MEMORIA:
#   - Cada script se ejecuta en un entorno local aislado para evitar
#     acumulación de objetos en RAM (previene crash OOM en el servidor)
#   - Se registra uso de RAM antes/después de cada paso
#   - El log de memoria se exporta a Excel y CSV en Resultados/
# ==============================================================================

library(here)
library(openxlsx)

# Cargar funciones globales (incluye funciones de monitoreo de RAM)
source(here("Scripts", "config.R"))
source(here("Scripts", "00_funciones_globales.R"))

# ==============================================================================
# 1. CONFIGURACIÓN DEL PIPELINE
# ==============================================================================

# -- Pasos a ejecutar (TRUE = ejecutar, FALSE = saltar) --
PASOS <- list(
  paso_01  = TRUE,   # Análisis de especies por árbol
  paso_02  = TRUE,   # Generar árbol medioide
  paso_03  = TRUE,   # Comparar medioide vs árboles (genera arboles_injertados/)
  paso_04  = TRUE,   # Calcular matriz Robinson-Foulds
  paso_06a = TRUE,   # Clustering — K-Means
  paso_06i = FALSE,  # Restaurar cache K-Means
  paso_07a = TRUE,   # Iterar - K-means
  paso_06b = TRUE,   # Clustering — PAM
  paso_06c = TRUE,   # Clustering — CLARA
  paso_07b = TRUE,   # Iterar - CLARA
  paso_06d = TRUE,   # Clustering — Mstknn
  paso_05a = TRUE,   # Reducción dimensional — UMAP
  paso_05b = TRUE,   # Reducción dimensional — t-SNE
  paso_05c = TRUE,   # Reducción dimensional — PCA
  paso_05d = TRUE,   # Reducción dimensional — nMDS
  paso_08  = TRUE,   # Etiquetas HGNC (requiere internet: mygene.info)
  paso_09  = TRUE,   # Enriquecimiento funcional (requiere internet: g:Profiler)
  paso_10  = TRUE,   # Gráfico de burbujas de enriquecimiento
  paso_11  = TRUE,   # Mapa de calor funcional
  paso_12  = TRUE    # Evaluación sesgo
)

# -- ¿Detener el pipeline si un paso falla? --
DETENER_EN_ERROR <- TRUE

# ==============================================================================
# 2. INFRAESTRUCTURA INTERNA
# ==============================================================================

# Registro de resultados por paso
log_pipeline <- data.frame(
  Paso    = character(),
  Script  = character(),
  Estado  = character(),
  Tiempo  = numeric(),
  Mensaje = character(),
  stringsAsFactors = FALSE
)

# Registro de memoria (se exporta a Excel/CSV al final)
log_memoria <- NULL

# Función auxiliar: ejecuta un script en entorno AISLADO y captura resultado
ejecutar_paso <- function(nombre_paso, ruta_script) {

  cat(sprintf(
    "\n%s\n>>> INICIANDO: %s\n%s\n",
    strrep("=", 70), nombre_paso, strrep("=", 70)
  ))

  if (!file.exists(ruta_script)) {
    msg <- paste("Script no encontrado:", ruta_script)
    cat("  [OMITIDO]", msg, "\n")
    return(list(estado = "OMITIDO", tiempo = 0, mensaje = msg))
  }

  t_inicio <- proc.time()

  resultado <- tryCatch({
    # ============================================================
    # CLAVE: ejecutar en entorno local aislado
    # Las variables creadas dentro del script se destruyen al terminar,
    # evitando acumulación de RAM entre pasos.
    # parent = globalenv() permite acceso a funciones globales y config.
    # ============================================================
    env_aislado <- new.env(parent = globalenv())
    source(ruta_script, echo = FALSE, local = env_aislado)
    
    # Destruir el entorno aislado explícitamente
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

  return(list(
    estado  = resultado$estado,
    tiempo  = t_elapsed,
    mensaje = resultado$mensaje
  ))
}

# Función para registrar en el log
registrar <- function(log_df, clave, nombre_script, res) {
  rbind(log_df, data.frame(
    Paso    = clave,
    Script  = nombre_script,
    Estado  = res$estado,
    Tiempo  = res$tiempo,
    Mensaje = res$mensaje,
    stringsAsFactors = FALSE
  ))
}

# ==============================================================================
# 3. EJECUCIÓN SECUENCIAL
# ==============================================================================

t_total_inicio <- proc.time()
cat("\n", strrep("*", 70), "\n")
cat("   PIPELINE AFOM — OrthoMaM\n")
cat("   Inicio:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("*", 70), "\n")

# Mapa: clave -> (nombre legible, ruta script relativa al proyecto)
definicion_pasos <- list(
  paso_01  = list("01 — Análisis por árbol",           "Scripts/01_analisis_por_arbol.R"),
  paso_02  = list("02 — Árbol medioide",                "Scripts/02_generar_arbol_medioide.R"),
  paso_03  = list("03 — Comparar medioide vs árboles", "Scripts/03_comparar_medioide_vs_arboles.R"),
  paso_04  = list("04 — Matriz Robinson-Foulds",       "Scripts/04_calcular_matriz_rf.R"),
  paso_06a = list("06a — Clustering K-Means",          "Scripts/06a_clustering_kmeans.R"),
  paso_06i = list("06i - Restaurar cache K-Means",     "Scripts/restaurar_cache_kmeans.R"),
  paso_07a = list("07a — Iterar K-Means",              "Scripts/07_iterar_k-means.R"),
  paso_06b = list("06b — Clustering PAM",              "Scripts/06b_clustering_pam.R"),
  paso_06c = list("06c — Clustering CLARA",            "Scripts/06c_clustering_clara.R"),
  paso_07b = list("07b — Iterar CLARA",                "Scripts/07_iterar_clara.R"),
  paso_06d = list("06c — Clustering Mstknn",           "Scripts/06d_clustering_mstknn.R"),
  paso_05a = list("05a — UMAP",                        "Scripts/05a_reduccion_dimensional_umap.R"),
  paso_05b = list("05b — t-SNE",                       "Scripts/05c_reduccion_dimensional_tSNE.R"),
  paso_05c = list("05c — PCA",                         "Scripts/05d_reduccion_dimensional_PCA.R"),
  paso_05d = list("05d — nMDS",                        "Scripts/05b_reduccion_dimensional_nMDS.R"),
  paso_08  = list("08 — Etiquetas HGNC",               "Scripts/08_etiquetas.R"),
  paso_09  = list("09 — Enriquecimiento funcional",    "Scripts/09_enriquecimiento.R"),
  paso_10  = list("10 — Gráfico de burbujas",          "Scripts/10_grafico_enriquecimiento.R"),
  paso_11  = list("11 — Mapa de calor",                "Scripts/11_mapa_de_calor.R"),
  paso_12  = list("12 - Evaluacion sesgo",             "Scripts/evaluacion_sesgo.R")
)

# Ejecutar cada paso según la configuración
pipeline_interrumpido <- FALSE

for (clave in names(definicion_pasos)) {

  if (!isTRUE(PASOS[[clave]])) {
    cat(sprintf("\n  [SALTADO] %s\n", definicion_pasos[[clave]][[1]]))
    log_pipeline <- registrar(
      log_pipeline, clave,
      basename(definicion_pasos[[clave]][[2]]),
      list(estado = "SALTADO", tiempo = 0, mensaje = "Desactivado en PASOS")
    )
    next
  }

  nombre <- definicion_pasos[[clave]][[1]]
  script <- here(definicion_pasos[[clave]][[2]])
  
  # --- Registrar RAM ANTES del paso ---
  log_memoria <- registrar_memoria(log_memoria, nombre, "ANTES")
  
  # --- Ejecutar paso ---
  res <- ejecutar_paso(nombre, script)
  
  # --- Forzar limpieza de memoria entre pasos ---
  limpiar_memoria()
  
  # --- Registrar RAM DESPUÉS del paso (post-gc) ---
  log_memoria <- registrar_memoria(log_memoria, nombre, "DESPUÉS")

  log_pipeline <- registrar(log_pipeline, clave, basename(script), res)

  if (res$estado == "ERROR" && DETENER_EN_ERROR) {
    cat("\n", strrep("!", 70), "\n")
    cat("  PIPELINE DETENIDO en:", nombre, "\n")
    cat("  Razón:", res$mensaje, "\n")
    cat(strrep("!", 70), "\n")
    pipeline_interrumpido <- TRUE
    break
  }
}

# ==============================================================================
# 4. INFORME FINAL
# ==============================================================================

t_total <- round((proc.time() - t_total_inicio)[["elapsed"]], 1)

cat("\n", strrep("*", 70), "\n")
cat("   RESUMEN DEL PIPELINE\n")
cat("   Fin:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("   Tiempo total: %.1f s  (%.1f min)\n", t_total, t_total / 60))
cat(strrep("*", 70), "\n\n")

# Tabla de resultados
cat(sprintf("  %-8s  %-45s  %-8s  %s\n", "Paso", "Script", "Estado", "Tiempo(s)"))
cat("  ", strrep("-", 70), "\n")
for (i in seq_len(nrow(log_pipeline))) {
  fila  <- log_pipeline[i, ]
  icono <- switch(fila$Estado,
    "OK"      = "[OK]   ",
    "ERROR"   = "[ERROR]",
    "SALTADO" = "[skip] ",
    "OMITIDO" = "[miss] ",
    "[?]    "
  )
  cat(sprintf("  %s  %-45s  %-8s  %.1f\n",
              icono, fila$Script, fila$Estado, fila$Tiempo))
}

# Mostrar errores si los hubo
errores <- log_pipeline[log_pipeline$Estado == "ERROR", ]
if (nrow(errores) > 0) {
  cat("\n  ERRORES ENCONTRADOS:\n")
  for (i in seq_len(nrow(errores))) {
    cat(sprintf("  - %s:\n    %s\n", errores$Script[i], errores$Mensaje[i]))
  }
}

if (!pipeline_interrumpido && all(log_pipeline$Estado %in% c("OK", "SALTADO", "OMITIDO"))) {
  cat("\n  Pipeline completado exitosamente!\n")
}

# ==============================================================================
# 5. EXPORTAR LOG DE EJECUCIÓN (CSV)
# ==============================================================================
ruta_log <- here("Resultados", paste0("pipeline_log_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
tryCatch(
  write.csv(log_pipeline, ruta_log, row.names = FALSE),
  error = function(e) cat("  Aviso: no se pudo guardar log CSV:", e$message, "\n")
)
cat("  Log guardado en:", ruta_log, "\n")

# ==============================================================================
# 6. EXPORTAR LOG DE MEMORIA (Excel + CSV)
# ==============================================================================
if (!is.null(log_memoria) && nrow(log_memoria) > 0) {
  
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # --- CSV ---
  ruta_mem_csv <- here("Resultados",
                       paste0("pipeline_memoria_", timestamp_str, ".csv"))
  tryCatch(
    write.csv(log_memoria, ruta_mem_csv, row.names = FALSE),
    error = function(e) cat("  Aviso: no se pudo guardar memoria CSV:", e$message, "\n")
  )
  
  # --- Excel con formato ---
  ruta_mem_xlsx <- here("Resultados",
                        paste0("pipeline_memoria_", timestamp_str, ".xlsx"))
  tryCatch({
    wb_mem <- createWorkbook()
    
    wb_mem <- agregar_hoja_formateada(
      wb           = wb_mem,
      nombre_hoja  = "Log_Memoria",
      titulo_tabla = paste0("Monitoreo de RAM — Pipeline AFOM (", NOMBRE_BDD, ")"),
      datos        = log_memoria,
      anchos_col   = c(22, 35, 12, 14, 16, 16, 10, 50)
    )
    
    # Agregar también el resumen del pipeline
    wb_mem <- agregar_hoja_formateada(
      wb           = wb_mem,
      nombre_hoja  = "Resumen_Pipeline",
      titulo_tabla = paste0("Resumen de Ejecución — ", NOMBRE_BDD,
                            " (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")"),
      datos        = log_pipeline,
      anchos_col   = c(12, 45, 10, 12, 50)
    )
    
    saveWorkbook(wb_mem, ruta_mem_xlsx, overwrite = TRUE)
    rm(wb_mem)
    
    cat("  Log de memoria (Excel):", ruta_mem_xlsx, "\n")
    cat("  Log de memoria (CSV)  :", ruta_mem_csv, "\n")
    
  }, error = function(e) {
    cat("  Aviso: no se pudo guardar memoria Excel:", e$message, "\n")
  })
  
  # --- Resumen de RAM en consola ---
  cat("\n", strrep("-", 70), "\n")
  cat("   RESUMEN DE MEMORIA\n")
  cat(strrep("-", 70), "\n")
  
  pasos_antes  <- log_memoria[log_memoria$Momento == "ANTES", ]
  pasos_despues <- log_memoria[log_memoria$Momento == "DESPUÉS", ]
  
  cat(sprintf("  Pico RAM R     : %.1f MB\n", max(log_memoria$RAM_R_MB, na.rm = TRUE)))
  
  if (any(!is.na(log_memoria$RAM_Sistema_MB))) {
    cat(sprintf("  Pico RAM Sistema: %.1f / %.1f MB (%.1f%%)\n",
                max(log_memoria$RAM_Sistema_MB, na.rm = TRUE),
                max(log_memoria$RAM_Total_MB, na.rm = TRUE),
                max(log_memoria$Pct_Uso, na.rm = TRUE)))
  }
  
  cat(strrep("-", 70), "\n\n")
  
} else {
  cat("  No se registraron datos de memoria.\n")
}
