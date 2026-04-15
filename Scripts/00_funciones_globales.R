# ==============================================================================
# ARCHIVO: 00_funciones_globales.R
# leer_bosque_zip
# calcular_matriz_rf
# ==============================================================================

library(ape)
library(TreeDist)
library(openxlsx)

# Función 1: Lectura y Estandarización de árboles
# Recibe una ruta, devuelve un objeto multiphylo estandarizado
#' Leer árboles filogenéticos desde archivos ZIP
#' @param directorio Ruta donde están los .zip
#' @param ext_interna Extensión del archivo de árbol dentro del zip (ej. ".rootree", ".nwk")

leer_bosque_zip <- function(directorio, ext_interna = ".rootree") {
  
  # Listar todos los zips en el directorio indicado
  zip_files <- list.files(path = directorio, pattern = "\\.zip$", full.names = TRUE)
  
  if (length(zip_files) == 0) stop("No se encontraron archivos .zip en ", directorio)
  
  # Leer e iterar
  arboles_raw <- lapply(zip_files, function(zfile) {
    tryCatch({
      # Ver qué hay dentro del ZIP
      archivos_internos <- unzip(zfile, list = TRUE)$Name
      
      # Buscar dinámicamente el archivo que termine con la extensión deseada
      archivo_objetivo <- archivos_internos[grepl(paste0("\\", ext_interna, "$"), archivos_internos, ignore.case = TRUE)][1]
      
      if (is.na(archivo_objetivo)) {
        warning("Saltando: No hay archivo ", ext_interna, " dentro de ", basename(zfile))
        return(NULL)
      }
      
      # Leer el árbol
      read.tree(unz(zfile, archivo_objetivo))
      
    }, error = function(e) {
      message("Error crítico en ", basename(zfile), ": ", e$message)
      return(NULL)
    })
  })
  
  # Asignar nombres limpios (sin el .zip) a los elementos de la lista
  names(arboles_raw) <- sub("\\.zip$", "", basename(zip_files))
  
  # Limpiar los nulos (archivos que fallaron) y convertir a multiPhylo
  arboles_validos <- arboles_raw[!sapply(arboles_raw, is.null)]
  class(arboles_validos) <- "multiPhylo"
  
  cat("Lectura completada:", length(arboles_validos), "árboles válidos extraídos.\n")
  return(arboles_validos)
}

# Función 2: Exportación y reporte en Excel
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