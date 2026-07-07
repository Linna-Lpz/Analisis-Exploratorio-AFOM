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
# ==============================================================================

library(here)

# ==============================================================================
# 1. CONFIGURACIÓN DEL PIPELINE
# ==============================================================================

# -- Pasos a ejecutar (TRUE = ejecutar, FALSE = saltar) --
PASOS <- list(
  paso_01  = TRUE,   # Análisis de especies por árbol
  paso_02  = TRUE,   # Generar árbol medioide
  paso_03  = TRUE,   # Comparar medioide vs árboles (genera arboles_podados/)
  paso_04  = TRUE,   # Calcular matriz Robinson-Foulds
  paso_05a = TRUE,   # Reducción dimensional — UMAP
  paso_05b = TRUE,   # Reducción dimensional — nMDS
  paso_05c = TRUE,   # Reducción dimensional — t-SNE
  paso_05d = TRUE,   # Reducción dimensional — PCA  <- proceso más pesado
  paso_06a = TRUE,   # Clustering — K-Means (versión principal)
  paso_08  = TRUE,   # Etiquetas HGNC (requiere internet: mygene.info)
  paso_09  = TRUE,   # Enriquecimiento funcional (requiere internet: g:Profiler)
  paso_10  = TRUE,   # Gráfico de burbujas de enriquecimiento
  paso_11  = TRUE    # Mapa de calor funcional
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

# Función auxiliar: ejecuta un script y captura resultado
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
    source(ruta_script, echo = FALSE, local = FALSE)
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
  paso_05a = list("05a — UMAP",                        "Scripts/05a_reduccion_dimensional_umap.R"),
  paso_05b = list("05b — nMDS",                        "Scripts/05b_reduccion_dimensional_nMDS.R"),
  paso_05c = list("05c — t-SNE",                       "Scripts/05c_reduccion_dimensional_tSNE.R"),
  paso_05d = list("05d — PCA",                         "Scripts/05d_reduccion_dimensional_PCA.R"),
  paso_06a = list("06a — Clustering K-Means",          "Scripts/06a_clustering_kmeans.R"),
  paso_08  = list("08 — Etiquetas HGNC",               "Scripts/08_etiquetas.R"),
  paso_09  = list("09 — Enriquecimiento funcional",    "Scripts/09_enriquecimiento.R"),
  paso_10  = list("10 — Gráfico de burbujas",          "Scripts/10_grafico_enriquecimiento.R"),
  paso_11  = list("11 — Mapa de calor",                "Scripts/11_mapa_de_calor.R")
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
  res    <- ejecutar_paso(nombre, script)

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

# Guardar log en CSV dentro de Resultados/ para trazabilidad
ruta_log <- here("Resultados", paste0("pipeline_log_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
tryCatch(
  write.csv(log_pipeline, ruta_log, row.names = FALSE),
  error = function(e) cat("  Aviso: no se pudo guardar log CSV:", e$message, "\n")
)
cat("  Log guardado en:", ruta_log, "\n\n")
