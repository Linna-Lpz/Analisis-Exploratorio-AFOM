# ==============================================================================
# ARCHIVO: 00_funciones_globales.R
# leer_bosque_zip
# ==============================================================================

library(ape)
library(TreeDist)
library(openxlsx)

# ==============================================================================
# FUNCIÓN AUXILIAR: Caché genérico con .rds
# ==============================================================================
#' Ejecuta una expresión y cachea el resultado como .rds.
#' Si el caché ya existe, lo carga directamente.
#'
#' @param ruta_cache  Ruta completa del archivo .rds de caché
#' @param expresion   Bloque de código que produce el objeto (entre llaves)
#' @param forzar      Si TRUE, recalcula aunque exista caché (default FALSE)
#' @return            El objeto cacheado o recién calculado
usar_cache <- function(ruta_cache, expresion, forzar = FALSE) {
  
  if (!forzar && file.exists(ruta_cache)) {
    cat(sprintf("  [caché] Cargando desde: %s\n", basename(ruta_cache)))
    return(readRDS(ruta_cache))
  }
  
  cat(sprintf("  [caché] Calculando y guardando: %s\n", basename(ruta_cache)))
  
  # Crear directorio si no existe
  dir_cache <- dirname(ruta_cache)
  if (!dir.exists(dir_cache)) dir.create(dir_cache, recursive = TRUE)
  
  resultado <- expresion
  saveRDS(resultado, file = ruta_cache)
  
  return(resultado)
}

# ==============================================================================
# FUNCIÓN 1: Lectura y Estandarización de árboles — con caché opcional
# ==============================================================================
#' Leer árboles filogenéticos desde archivos ZIP, con caché .rds opcional
#' @param directorio    Ruta donde están los .zip
#' @param ext_interna   Extensión del archivo de árbol dentro del zip (ej. ".rootree", ".nwk")
#' @param dir_cache     Directorio donde guardar/leer el caché. NULL desactiva el caché.
#' @param forzar_recalc Si TRUE, ignora el caché existente y recalcula siempre.
leer_bosque_zip <- function(directorio,
                            ext_interna   = ".rootree",
                            dir_cache     = NULL,
                            forzar_recalc = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Verificar caché antes de hacer cualquier lectura
  # ---------------------------------------------------------------------------
  if (!is.null(dir_cache)) {
    if (!dir.exists(dir_cache)) dir.create(dir_cache, recursive = TRUE)
    
    # Nombre del caché derivado del directorio fuente (evita colisiones)
    nombre_cache <- paste0("bosque_", digest::digest(directorio), ".rds")
    ruta_cache   <- file.path(dir_cache, nombre_cache)
    
    if (file.exists(ruta_cache) && !forzar_recalc) {
      cat("Bosque encontrado en caché. Cargando desde:", ruta_cache, "\n")
      arboles_validos <- readRDS(ruta_cache)
      cat("Cargados", length(arboles_validos), "árboles desde caché.\n")
      return(arboles_validos)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Lectura real desde los ZIPs (solo si no hay caché válido)
  # ---------------------------------------------------------------------------
  zip_files <- list.files(path = directorio, pattern = "\\.zip$", full.names = TRUE)
  
  if (length(zip_files) == 0) stop("No se encontraron archivos .zip en ", directorio)
  
  arboles_raw <- lapply(zip_files, function(zfile) {
    tryCatch({
      archivos_internos <- unzip(zfile, list = TRUE)$Name
      
      archivo_objetivo <- archivos_internos[
        grepl(paste0("\\", ext_interna, "$"), archivos_internos, ignore.case = TRUE)
      ][1]
      
      if (is.na(archivo_objetivo)) {
        warning("Saltando: No hay archivo ", ext_interna, " dentro de ", basename(zfile))
        return(NULL)
      }
      
      read.tree(unz(zfile, archivo_objetivo))
      
    }, error = function(e) {
      message("Error crítico en ", basename(zfile), ": ", e$message)
      return(NULL)
    })
  })
  
  names(arboles_raw) <- sub("\\.zip$", "", basename(zip_files))
  
  arboles_validos <- arboles_raw[!sapply(arboles_raw, is.null)]
  class(arboles_validos) <- "multiPhylo"
  
  cat("Lectura completada:", length(arboles_validos), "árboles válidos extraídos.\n")
  
  # ---------------------------------------------------------------------------
  # Guardar en caché si corresponde
  # ---------------------------------------------------------------------------
  if (!is.null(dir_cache)) {
    saveRDS(arboles_validos, file = ruta_cache)
    cat("Bosque guardado en caché:", ruta_cache, "\n")
  }
  
  return(arboles_validos)
}


# ==============================================================================
# FUNCIÓN 2: Exportación y reporte en Excel
# ==============================================================================
#' @param wb
#' @param nombre_hoja
#' @param titulo_tabla
#' @param datos
#' @param anchos_col

agregar_hoja_formateada <- function(wb, nombre_hoja, titulo_tabla, datos, anchos_col = "auto") {
  
  # Definir estilos estandarizados internamente
  estilo_titulo <- createStyle(fontColour = "#FFFFFF", fgFill = "#2E4057",
                               textDecoration = "bold", halign = "center",
                               fontSize = 11)
  estilo_header <- createStyle(fontColour = "#FFFFFF", fgFill = "#4472C4",
                               textDecoration = "bold", halign = "center",
                               border = "Bottom", borderColour = "#FFFFFF")
  estilo_celdas <- createStyle(halign = "left", border = "TopBottomLeftRight",
                               borderColour = "#D9D9D9")
  
  # Inicializar la hoja
  addWorksheet(wb, nombre_hoja)
  n_cols <- ncol(datos)
  
  # Escribir y dar estilo al Título General (Fila 1)
  writeData(wb, nombre_hoja, titulo_tabla, startRow = 1, startCol = 1)
  addStyle(wb, nombre_hoja, estilo_titulo, rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  mergeCells(wb, nombre_hoja, cols = 1:n_cols, rows = 1)
  
  # Escribir Datos (Fila 2 en adelante)
  writeData(wb, nombre_hoja, datos, startRow = 2, startCol = 1, rowNames = FALSE)
  
  # Aplicar Estilo de Cabecera (Fila 2)
  addStyle(wb, nombre_hoja, estilo_header, rows = 2, cols = 1:n_cols, gridExpand = TRUE)
  
  # Aplicar Estilo de Celdas (Fila 3 en adelante)
  # Se añade un condicional por si el dataframe está vacío (ej. no hay duplicados)
  if (nrow(datos) > 0) {
    addStyle(wb, nombre_hoja, estilo_celdas, 
             rows = 3:(2 + nrow(datos)), cols = 1:n_cols, gridExpand = TRUE)
  }
  
  # Ajuste de anchos de columna
  if (length(anchos_col) == 1 && anchos_col == "auto") {
    setColWidths(wb, nombre_hoja, cols = 1:n_cols, widths = "auto")
  } else {
    setColWidths(wb, nombre_hoja, cols = 1:n_cols, widths = anchos_col)
  }
  
  return(wb)
}
message("Funciones globales cargadas correctamente.")