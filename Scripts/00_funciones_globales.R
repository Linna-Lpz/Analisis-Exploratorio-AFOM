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

# ==============================================================================
# FUNCIÓN 3: Auxiliar para leer hoja de asignaciones desde Excel
# ==============================================================================
leer_asignaciones_xlsx <- function(ruta_xlsx, metodo) {
  
  if (!file.exists(ruta_xlsx)) {
    warning(sprintf("No se encontró el archivo de %s: %s", metodo, ruta_xlsx))
    return(NULL)
  }
  
  # Verificar que la hoja existe dentro del libro
  hojas_disponibles <- getSheetNames(ruta_xlsx)
  
  if (!"Asignaciones_K_Optimo" %in% hojas_disponibles) {
    warning(sprintf("Hoja 'Asignaciones_K_Optimo' no encontrada en %s.\nHojas disponibles: %s",
                    ruta_xlsx, paste(hojas_disponibles, collapse = ", ")))
    return(NULL)
  }
  
  df <- read.xlsx(ruta_xlsx,
                  sheet     = "Asignaciones_K_Optimo",
                  startRow  = 2,          # fila 1 = título de agregar_hoja_formateada
                  colNames  = TRUE)
  
  # Estandarizar: conservar solo Arbol + Cluster
  df$Arbol   <- as.character(df$Arbol)
  df$Cluster <- as.factor(df$Cluster)
  
  cat(sprintf("  [%s] %d árboles leídos desde %s\n", metodo, nrow(df), basename(ruta_xlsx)))
  
  return(df[, c("Arbol", "Cluster")])
}
  
  # ==============================================================================
  # FUNCIÓN 4: Leer y unir etiquetas de clustering al dataframe de coordenadas
  # ==============================================================================
  #' Lee las asignaciones de K-Means, PAM y CLARA desde sus Excel respectivos
  #' y las une al dataframe de coordenadas por la columna "Arbol".
  #'
  #' @param coords_df   Dataframe con coordenadas (debe tener columna "Arbol")
  #' @param dir_results Directorio donde están los Excel de resultados (DIR_RESULTS)
  #' @param nombre_kmeans Nombre del archivo Excel de K-Means (default: "kmeans_resultados.xlsx")
  #' @param nombre_pam    Nombre del archivo Excel de PAM   (default: "pam_resultados.xlsx")
  #' @param nombre_clara  Nombre del archivo Excel de CLARA (default: "clara_resultados.xlsx")
  #' @return Lista con dos elementos:
  #'         $coords_df : dataframe enriquecido con columnas Cluster_KMeans, Cluster_PAM, Cluster_CLARA
  #'         $k_optimos : named vector con el k usado por cada método (para subtítulos de gráficos)
  
  unir_etiquetas_clustering <- function(coords_df,
                                        dir_results,
                                        nombre_kmeans = "kmeans_resultados.xlsx",
                                        nombre_pam    = "pam_resultados.xlsx",
                                        nombre_clara  = "clara_resultados.xlsx") {
    
    cat("\n=== UNIENDO ETIQUETAS DE CLUSTERING ===\n")
    
    # Verificar que coords_df tiene la columna requerida
    if (!"Arbol" %in% colnames(coords_df)) {
      stop("coords_df debe contener una columna llamada 'Arbol'.")
    }
    
    # Leer los tres métodos
    kmeans_asig <- leer_asignaciones_xlsx(file.path(dir_results, nombre_kmeans), "K-Means")
    pam_asig    <- leer_asignaciones_xlsx(file.path(dir_results, nombre_pam),    "PAM")
    clara_asig  <- leer_asignaciones_xlsx(file.path(dir_results, nombre_clara),  "CLARA")
    
    # Unir cada método renombrando la columna Cluster antes del merge
    if (!is.null(kmeans_asig)) {
      kmeans_asig <- setNames(kmeans_asig, c("Arbol", "Cluster_KMeans"))
      coords_df   <- merge(coords_df, kmeans_asig, by = "Arbol", all.x = TRUE)
    }
    
    if (!is.null(pam_asig)) {
      pam_asig  <- setNames(pam_asig, c("Arbol", "Cluster_PAM"))
      coords_df <- merge(coords_df, pam_asig, by = "Arbol", all.x = TRUE)
    }
    
    if (!is.null(clara_asig)) {
      clara_asig <- setNames(clara_asig, c("Arbol", "Cluster_CLARA"))
      coords_df  <- merge(coords_df, clara_asig, by = "Arbol", all.x = TRUE)
    }
    
    # Vector de k óptimos para usar en subtítulos de gráficos
    k_optimos <- c(
      KMeans = if (!is.null(kmeans_asig)) length(unique(kmeans_asig$Cluster_KMeans)) else NA,
      PAM    = if (!is.null(pam_asig))    length(unique(pam_asig$Cluster_PAM))       else NA,
      CLARA  = if (!is.null(clara_asig))  length(unique(clara_asig$Cluster_CLARA))   else NA
    )
    
    cat("\nColumnas del dataframe enriquecido:\n")
    print(colnames(coords_df))
    cat("Primeras filas:\n")
    print(head(coords_df, 5))
    
    return(list(coords_df = coords_df, k_optimos = k_optimos))
  }
