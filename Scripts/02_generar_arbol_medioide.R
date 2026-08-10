# =============================================================================
# FILTRO DE ÁRBOLES Y CÁLCULO DEL ÁRBOL MEDIOIDE
# Lee el CORE_SET desde el Excel generado en el script 01,
# filtra los árboles correspondientes y calcula el medioide.
# =============================================================================

library(ape)
library(TreeDist)
library(openxlsx)

# --- Cargar configuración y funciones globales ---
source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# Definir ruta del Excel generado por script 01 (permite ejecución aislada)
ruta_especies_output <- file.path(DIR_RESULTS,
                                  paste0("especies_analisis_", NOMBRE_BDD, ".xlsx"))

if (!file.exists(ruta_especies_output)) {
  stop("No se encontró el Excel de análisis de especies: ", ruta_especies_output,
       "\nEjecuta primero el script 01_analisis_por_arbol.R")
}

# =============================================================================
# LEER CORE_SET DESDE EL EXCEL DEL SCRIPT 01
# Lee la primera fila de datos de "Resumen_Conteos" para obtener
# N_de_especies y N_de_arboles del conjunto más frecuente
# =============================================================================
cat("\n=== LEYENDO CORE_SET DESDE EXCEL ===\n")

resumen <- read.xlsx(ruta_especies_output, sheet = "Resumen_Conteos", startRow = 2)
# La hoja tiene título en fila 1, headers en fila 2, datos desde fila 3
# startRow = 2 hace que openxlsx trate la fila 2 como header y lea desde fila 3

primera_fila <- resumen[1, ]
CORE_SET      <- as.integer(primera_fila$N_de_especies)
N_ARBOLES_ESP <- as.integer(primera_fila$N_de_arboles)

cat("CORE_SET leído   :", CORE_SET, "especies\n")
cat("Árboles esperados:", N_ARBOLES_ESP, "\n")

# =============================================================================
# CARGAR ÁRBOLES DEL CORE SET (DESDE CACHÉ SI EXISTE)
# =============================================================================
ruta_cache_core_set <- file.path(DIR_CACHE, "conjunto_core_set.rds")

if (file.exists(ruta_cache_core_set)) {
  
  cat("\n=== CARGANDO CONJUNTO CORE SET DESDE CACHÉ ===\n")
  arboles <- readRDS(ruta_cache_core_set)
  cat("Árboles cargados desde caché:", length(arboles), "\n")
  
  if (length(arboles) != N_ARBOLES_ESP) {
    warning(
      "El conjunto en caché tiene ", length(arboles), " árboles, ",
      "pero el Excel espera ", N_ARBOLES_ESP, ". ",
      "Si cambiaron los datos de entrada, elimina el archivo de caché:\n  ",
      ruta_cache_core_set
    )
  }
  
} else {
  
  # ===========================================================================
  # LEER TODOS LOS ÁRBOLES Y CONTAR ESPECIES
  # ===========================================================================
  cat("\n=== CONTANDO ESPECIES POR ÁRBOL ===\n")
  
  zip_files <- list.files(
    file.path(DIR_INPUT, CARPETA_ARBOLES),
    pattern    = "\\.zip$",
    full.names = TRUE
  )
  cat("ZIPs encontrados:", length(zip_files), "\n")
  
  conteo <- data.frame(
    ruta       = zip_files,
    nombre     = basename(zip_files),
    n_especies = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(zip_files)) {
    tryCatch({
      archivos_internos <- unzip(zip_files[i], list = TRUE)$Name
      rootree_file <- archivos_internos[
        grepl(EXTENSION_ARBOLES, archivos_internos, fixed = TRUE)
      ][1]
      if (!is.na(rootree_file)) {
        arbol <- read.tree(unz(zip_files[i], rootree_file))
        conteo$n_especies[i] <- length(arbol$tip.label)
      }
    }, error = function(e) {
      message("  Error en: ", basename(zip_files[i]), " — ", e$message)
    })
  }
  
  cat("Árboles leídos correctamente:",
      sum(!is.na(conteo$n_especies)), "/", nrow(conteo), "\n")
  
  # ===========================================================================
  # FILTRAR ÁRBOLES CON EXACTAMENTE CORE_SET ESPECIES
  # ===========================================================================
  cat("\n=== FILTRANDO ÁRBOLES CON", CORE_SET, "ESPECIES ===\n")
  
  seleccionados <- conteo[
    !is.na(conteo$n_especies) & conteo$n_especies == CORE_SET,
  ]
  cat("Árboles encontrados con", CORE_SET, "especies:", nrow(seleccionados), "\n")
  
  if (nrow(seleccionados) != N_ARBOLES_ESP) {
    warning(
      "Se esperaban ", N_ARBOLES_ESP, " árboles según el Excel, ",
      "pero se encontraron ", nrow(seleccionados), "."
    )
  }
  
  if (nrow(seleccionados) == 0) {
    cat("\nDistribución de especies disponibles:\n")
    print(sort(table(conteo$n_especies), decreasing = TRUE))
    stop("No se encontraron árboles con ", CORE_SET, " especies. ",
         "Revisa el archivo Excel o el directorio de entrada.")
  }
  
  # ===========================================================================
  # CARGAR EN MEMORIA LOS ÁRBOLES SELECCIONADOS
  # ===========================================================================
  cat("\n=== CARGANDO ÁRBOLES SELECCIONADOS EN MEMORIA ===\n")
  
  arboles_raw <- leer_bosque_zip(
    directorio  = file.path(DIR_INPUT, CARPETA_ARBOLES),
    ext_interna = EXTENSION_ARBOLES,
    dir_cache = DIR_CACHE
  )
  
  # Conservar solo los árboles cuyo nombre de ZIP está en los seleccionados
  nombres_seleccionados <- tools::file_path_sans_ext(seleccionados$nombre)
  arboles <- arboles_raw[names(arboles_raw) %in% nombres_seleccionados]
  arboles <- arboles[!sapply(arboles, is.null)]
  class(arboles) <- "multiPhylo"
  
  cat("Árboles cargados para el análisis:", length(arboles), "\n")
  
  # --- Guardar en caché el conjunto filtrado con el core set de especies ---
  saveRDS(arboles, ruta_cache_core_set)
  cat("Conjunto core set guardado en caché:", ruta_cache_core_set, "\n")
}

# =============================================================================
# CALCULAR MATRIZ DE DISTANCIAS Y ENCONTRAR EL MEDIOIDE
# =============================================================================
cat("\n=== PASO 5: CALCULANDO MATRIZ DE DISTANCIAS ===\n")

matriz_distancias <- RobinsonFoulds(arboles, normalize = TRUE)
matriz_cuadrada   <- as.matrix(matriz_distancias)
suma_por_fila     <- rowSums(matriz_cuadrada)

indice_medioide   <- which.min(suma_por_fila)
nombre_medioide   <- names(arboles)[indice_medioide]
arbol_referencia  <- arboles[[nombre_medioide]]

cat("\n=== ÁRBOL MEDIOIDE ===\n")
cat("Índice         :", indice_medioide, "\n")
cat("Nombre         :", nombre_medioide, "\n")
cat("Suma distancias:", round(suma_por_fila[indice_medioide], 4), "\n")
cat("N° de tips     :", length(arbol_referencia$tip.label), "\n")

# =============================================================================
# RANKING COMPLETO Y EXPORTACIÓN A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RANKING ===\n")

nombres_internos <- sapply(seq_along(arboles), function(i) {
  nombre_interno <- arboles[[i]]$name
  if (is.null(nombre_interno) || nombre_interno == "") {
    return(names(arboles)[i])
  }
  return(nombre_interno)
})

ranking <- data.frame(
  posicion     = rank(suma_por_fila, ties.method = "first"),
  nombre_arbol = nombres_internos,
  suma_dist    = round(suma_por_fila, 4),
  dist_media   = round(suma_por_fila / (length(arboles) - 1), 4),
  stringsAsFactors = FALSE
)
ranking <- ranking[order(ranking$posicion), ]
rownames(ranking) <- NULL

wb <- createWorkbook()

titulo_dinamico <- paste0("RANKING MEDIOIDE — CORE SET: ", CORE_SET, " especies")

# Añadir hoja con anchos automáticos
agregar_hoja_formateada(wb, "Ranking_Medioide", titulo_dinamico, 
                        ranking, anchos_col = "auto")

# Guardar el libro
ruta_medioide_output <- file.path(DIR_RESULTS, paste0("ranking_medioide_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_medioide_output, overwrite = TRUE)

cat("Ranking guardado:", ruta_medioide_output, "\n")
cat("\n=== COMPLETADO ===\n")