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
leer_asignaciones_xlsx <- function(ruta_xlsx, metodo,
                                   nombre_hoja = "Asignaciones_K_Optimo") {
  
  if (!file.exists(ruta_xlsx)) {
    warning(sprintf("No se encontró el archivo de %s: %s", metodo, ruta_xlsx))
    return(NULL)
  }
  
  # Verificar que la hoja existe dentro del libro
  hojas_disponibles <- getSheetNames(ruta_xlsx)
  
  if (!nombre_hoja %in% hojas_disponibles) {
    warning(sprintf("Hoja '%s' no encontrada en %s.\nHojas disponibles: %s",
                    nombre_hoja, ruta_xlsx, paste(hojas_disponibles, collapse = ", ")))
    return(NULL)
  }
  
  df <- read.xlsx(ruta_xlsx,
                  sheet     = nombre_hoja,
                  startRow  = 2,          # fila 1 = título de agregar_hoja_formateada
                  colNames  = TRUE)
  
  # Estandarizar: conservar solo Arbol + Cluster
  df$Arbol   <- as.character(df$Arbol)
  df$Cluster <- as.factor(df$Cluster)
  
  cat(sprintf("  [%s] %d árboles leídos desde %s (hoja: %s)\n",
              metodo, nrow(df), basename(ruta_xlsx), nombre_hoja))
  
  return(df[, c("Arbol", "Cluster")])
}

# ==============================================================================
# FUNCIÓN 4: Leer y unir etiquetas de clustering al dataframe de coordenadas
# ==============================================================================
#' Lee las asignaciones de K-Means, PAM y CLARA desde sus Excel respectivos
#' y las une al dataframe de coordenadas por la columna "Arbol".
#' Además de las asignaciones del k óptimo, lee también hojas de k adicionales
#' (por defecto k=10 y k=15) cuando existen.
#'
#' @param coords_df   Dataframe con coordenadas (debe tener columna "Arbol")
#' @param dir_results Directorio donde están los Excel de resultados (DIR_RESULTS)
#' @param nombre_kmeans Nombre del archivo Excel de K-Means (default: "kmeans_resultados.xlsx")
#' @param nombre_pam    Nombre del archivo Excel de PAM   (default: "pam_resultados.xlsx")
#' @param nombre_clara  Nombre del archivo Excel de CLARA (default: "clara_resultados.xlsx")
#' @param k_extra       Vector de k adicionales a leer (default: c(10, 15))
#' @return Lista con tres elementos:
#'         $coords_df : dataframe enriquecido con columnas Cluster_KMeans, Cluster_PAM, Cluster_CLARA
#'                      y columnas adicionales como Cluster_KMeans_K10, Cluster_CLARA_K15, etc.
#'         $k_optimos : named vector con el k usado por cada método (para subtítulos de gráficos)
#'         $k_extra   : vector de k adicionales efectivamente leídos

unir_etiquetas_clustering <- function(coords_df,
                                      dir_results,
                                      nombre_kmeans = "kmeans_subdivision_iterativa.xlsx",
                                      nombre_pam    = "pam_resultados.xlsx",
                                      nombre_clara  = "clara_subdivision_iterativa.xlsx",
                                      k_extra       = c(10, 15)) {
  
  cat("\n=== UNIENDO ETIQUETAS DE CLUSTERING ===\n")
  
  if (!"Arbol" %in% colnames(coords_df))
    stop("coords_df debe contener una columna llamada 'Arbol'.")
  
  archivos_metodos <- list(
    KMeans = file.path(dir_results, nombre_kmeans),
    PAM    = file.path(dir_results, nombre_pam),
    CLARA  = file.path(dir_results, nombre_clara)
  )
  
  # --- K ÓPTIMO ---
  kmeans_asig <- leer_asignaciones_xlsx(archivos_metodos[["KMeans"]], "K-Means")
  pam_asig    <- leer_asignaciones_xlsx(archivos_metodos[["PAM"]],    "PAM")
  clara_asig  <- leer_asignaciones_xlsx(archivos_metodos[["CLARA"]],  "CLARA")
  
  if (!is.null(kmeans_asig)) {
    coords_df <- merge(coords_df, setNames(kmeans_asig, c("Arbol", "Cluster_KMeans")),
                       by = "Arbol", all.x = TRUE)
  }
  if (!is.null(pam_asig)) {
    coords_df <- merge(coords_df, setNames(pam_asig, c("Arbol", "Cluster_PAM")),
                       by = "Arbol", all.x = TRUE)
  }
  if (!is.null(clara_asig)) {
    coords_df <- merge(coords_df, setNames(clara_asig, c("Arbol", "Cluster_CLARA")),
                       by = "Arbol", all.x = TRUE)
  }
  
  # --- K EXTRA ---
  if (length(k_extra) > 0) {
    cat("\n--- Leyendo asignaciones para k extra:", paste(k_extra, collapse = ", "), "---\n")
    
    for (ke in k_extra) {
      nombre_hoja_extra <- paste0("Asignaciones_K_", ke)
      
      km_extra <- leer_asignaciones_xlsx(archivos_metodos[["KMeans"]], "K-Means",
                                         nombre_hoja = nombre_hoja_extra)
      if (!is.null(km_extra)) {
        coords_df <- merge(coords_df,
                           setNames(km_extra, c("Arbol", paste0("Cluster_KMeans_K", ke))),
                           by = "Arbol", all.x = TRUE)
      }
      
      pam_extra <- leer_asignaciones_xlsx(archivos_metodos[["PAM"]], "PAM",
                                          nombre_hoja = nombre_hoja_extra)
      if (!is.null(pam_extra)) {
        coords_df <- merge(coords_df,
                           setNames(pam_extra, c("Arbol", paste0("Cluster_PAM_K", ke))),
                           by = "Arbol", all.x = TRUE)
      }
      
      cl_extra <- leer_asignaciones_xlsx(archivos_metodos[["CLARA"]], "CLARA",
                                         nombre_hoja = nombre_hoja_extra)
      if (!is.null(cl_extra)) {
        coords_df <- merge(coords_df,
                           setNames(cl_extra, c("Arbol", paste0("Cluster_CLARA_K", ke))),
                           by = "Arbol", all.x = TRUE)
      }
    }
  }
  
  # --- K-MEANS SUBDIVISIÓN ITERATIVA --- (fuera del loop, se lee una sola vez)
  ruta_sub_iter_k <- file.path(dir_results, "kmeans_subdivision_iterativa.xlsx")
  if (file.exists(ruta_sub_iter_k)) {
    hojas_sub <- getSheetNames(ruta_sub_iter_k)
    if ("Asignaciones_K_Optimo" %in% hojas_sub) {
      sub_iter_df <- read.xlsx(ruta_sub_iter_k,
                               sheet    = "Asignaciones_K_Optimo",
                               startRow = 2,
                               colNames = TRUE)
      sub_iter_df$Arbol         <- as.character(sub_iter_df$Arbol)
      sub_iter_df$Cluster_Final <- as.factor(sub_iter_df$Cluster_Final)
      coords_df <- merge(coords_df,
                         setNames(sub_iter_df[, c("Arbol", "Cluster_Final")],
                                  c("Arbol", "Cluster_Kmeans_Iter")),
                         by = "Arbol", all.x = TRUE)
      cat(sprintf("  [KMEANS iterativa] %d árboles leídos\n", nrow(sub_iter_df)))
    }
  }
  
  # --- CLARA SUBDIVISIÓN ITERATIVA --- (fuera del loop, se lee una sola vez)
  ruta_sub_iter_c <- file.path(dir_results, "clara_subdivision_iterativa.xlsx")
  if (file.exists(ruta_sub_iter_c)) {
    hojas_sub <- getSheetNames(ruta_sub_iter_c)
    if ("Asignaciones_K_Optimo" %in% hojas_sub) {
      sub_iter_df <- read.xlsx(ruta_sub_iter_c,
                               sheet    = "Asignaciones_K_Optimo",
                               startRow = 2,
                               colNames = TRUE)
      sub_iter_df$Arbol         <- as.character(sub_iter_df$Arbol)
      sub_iter_df$Cluster_Final <- as.factor(sub_iter_df$Cluster_Final)
      coords_df <- merge(coords_df,
                         setNames(sub_iter_df[, c("Arbol", "Cluster_Final")],
                                  c("Arbol", "Cluster_CLARA_Iter")),
                         by = "Arbol", all.x = TRUE)
      cat(sprintf("  [CLARA iterativa] %d árboles leídos\n", nrow(sub_iter_df)))
    }
  }
  
  # --- MSTKNN SUBDIVISIÓN ITERATIVA --- (fuera del loop, se lee una sola vez)
  ruta_sub_iter_m <- file.path(dir_results, "mstknn_subdivision_iterativa.xlsx")
  if (file.exists(ruta_sub_iter_m)) {
    hojas_sub <- getSheetNames(ruta_sub_iter_m)
    if ("Asignaciones_Finales" %in% hojas_sub) {
      sub_iter_df <- read.xlsx(ruta_sub_iter_m,
                               sheet    = "Asignaciones_Finales",
                               startRow = 2,
                               colNames = TRUE)
      sub_iter_df$Arbol         <- as.character(sub_iter_df$Arbol)
      sub_iter_df$Cluster_Final <- as.factor(sub_iter_df$Cluster_Final)
      coords_df <- merge(coords_df,
                         setNames(sub_iter_df[, c("Arbol", "Cluster_Final")],
                                  c("Arbol", "Cluster_MSTKNN_Iter")),
                         by = "Arbol", all.x = TRUE)
      cat(sprintf("  [MSTKNN iterativa] %d árboles leídos\n", nrow(sub_iter_df)))
    }
  }
  
  # --- k_optimos ---
  k_optimos <- c(
    KMeans      = if ("Cluster_KMeans"      %in% colnames(coords_df))
      length(unique(na.omit(coords_df$Cluster_KMeans)))      else NA,
    PAM         = if ("Cluster_PAM"         %in% colnames(coords_df))
      length(unique(na.omit(coords_df$Cluster_PAM)))         else NA,
    CLARA       = if ("Cluster_CLARA"       %in% colnames(coords_df))
      length(unique(na.omit(coords_df$Cluster_CLARA)))       else NA,
    CLARA_Iter  = if ("Cluster_CLARA_Iter"  %in% colnames(coords_df))
      length(unique(na.omit(coords_df$Cluster_CLARA_Iter)))  else NA,
    KMeans_Iter = if ("Cluster_KMeans_Iter" %in% colnames(coords_df))
      length(unique(na.omit(coords_df$Cluster_KMeans_Iter))) else NA,
    MSTKNN_Iter = if ("Cluster_MSTKNN_Iter" %in% colnames(coords_df))
      length(unique(na.omit(coords_df$Cluster_MSTKNN_Iter))) else NA
  )
  
  cat("\nColumnas del dataframe enriquecido:\n")
  print(colnames(coords_df))
  cat("Primeras filas:\n")
  print(head(coords_df, 5))
  
  return(list(coords_df = coords_df, k_optimos = k_optimos, k_extra = k_extra))
}

# ==============================================================================
# FUNCIÓN 5: Obtener uso de memoria (cross-platform: Linux/Windows)
# ==============================================================================
#' Retorna métricas de memoria en MB.
#' @return Named list: ram_r_mb, ram_sistema_mb, ram_total_mb, pct_uso
obtener_memoria_mb <- function() {
  
  # RAM usada por el proceso R (siempre disponible)
  gc_info   <- gc(verbose = FALSE, reset = FALSE)
  ram_r_mb  <- sum(gc_info[, 2]) # La columna 2 siempre es 'used (Mb)'
  
  
  # RAM del sistema operativo
  ram_sistema_mb <- NA_real_
  ram_total_mb   <- NA_real_
  
  if (.Platform$OS.type == "unix" && file.exists("/proc/meminfo")) {
    # --- Linux (servidor) ---
    meminfo <- readLines("/proc/meminfo", n = 10)
    
    total_line <- grep("^MemTotal:", meminfo, value = TRUE)
    avail_line <- grep("^MemAvailable:", meminfo, value = TRUE)
    
    if (length(total_line) > 0) {
      ram_total_mb <- as.numeric(gsub("[^0-9]", "", total_line)) / 1024
    }
    if (length(avail_line) > 0) {
      ram_libre_mb   <- as.numeric(gsub("[^0-9]", "", avail_line)) / 1024
      ram_sistema_mb <- ram_total_mb - ram_libre_mb
    }
    
  } else if (.Platform$OS.type == "windows") {
    # --- Windows (local) ---
    tryCatch({
      wmic <- system2("wmic", c("OS", "get",
                                "FreePhysicalMemory,TotalVisibleMemorySize",
                                "/value"),
                      stdout = TRUE, stderr = TRUE)
      wmic <- wmic[nchar(trimws(wmic)) > 0]
      
      free_line  <- grep("^FreePhysicalMemory=", wmic, value = TRUE)
      total_line <- grep("^TotalVisibleMemorySize=", wmic, value = TRUE)
      
      if (length(total_line) > 0) {
        ram_total_mb <- as.numeric(gsub("[^0-9]", "", total_line)) / 1024
      }
      if (length(free_line) > 0) {
        ram_libre_mb   <- as.numeric(gsub("[^0-9]", "", free_line)) / 1024
        ram_sistema_mb <- ram_total_mb - ram_libre_mb
      }
    }, error = function(e) NULL)
  }
  
  pct_uso <- if (!is.na(ram_sistema_mb) && !is.na(ram_total_mb) && ram_total_mb > 0) {
    round(ram_sistema_mb / ram_total_mb * 100, 1)
  } else {
    NA_real_
  }
  
  list(
    ram_r_mb       = round(ram_r_mb, 1),
    ram_sistema_mb = round(ram_sistema_mb, 1),
    ram_total_mb   = round(ram_total_mb, 1),
    pct_uso        = pct_uso
  )
}

# ==============================================================================
# FUNCIÓN 6: Limpiar memoria con reporte
# ==============================================================================
#' Ejecuta gc() y reporta cuánta RAM se liberó.
#' @param silencioso Si TRUE, no imprime nada en consola
#' @return Lista con ram_antes y ram_despues (en MB)
limpiar_memoria <- function(silencioso = FALSE) {
  
  antes <- obtener_memoria_mb()
  gc(verbose = FALSE, reset = TRUE)
  gc(verbose = FALSE)  # doble gc para liberar objetos de finalizer
  despues <- obtener_memoria_mb()
  
  liberado <- round(antes$ram_r_mb - despues$ram_r_mb, 1)
  
  if (!silencioso) {
    cat(sprintf("  [gc] R: %.1f MB -> %.1f MB (liberados: %.1f MB)\n",
                antes$ram_r_mb, despues$ram_r_mb, liberado))
  }
  
  invisible(list(antes = antes, despues = despues, liberado_mb = liberado))
}

# ==============================================================================
# FUNCIÓN 7: Registrar memoria en log dataframe
# ==============================================================================
#' Agrega una fila al dataframe de log de memoria.
#' @param log_df    Dataframe de log existente (o NULL para inicializar)
#' @param paso      Nombre del paso (ej. "05a — UMAP")
#' @param momento   "ANTES" o "DESPUÉS"
#' @param imprimir  Si TRUE, imprime en consola
#' @return Dataframe actualizado con la nueva fila
registrar_memoria <- function(log_df = NULL, paso, momento, imprimir = TRUE) {
  
  mem <- obtener_memoria_mb()
  
  # Top 5 objetos más pesados en el entorno global
  objetos <- ls(envir = globalenv())
  if (length(objetos) > 0) {
    tamanios <- sapply(objetos, function(x) {
      tryCatch(object.size(get(x, envir = globalenv())),
               error = function(e) 0)
    })
    top5 <- head(sort(tamanios, decreasing = TRUE), 5)
    top5_txt <- paste(
      sprintf("%s(%.0fMB)", names(top5), top5 / 1024^2),
      collapse = ", "
    )
  } else {
    top5_txt <- ""
  }
  
  nueva_fila <- data.frame(
    Timestamp      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    Paso           = paso,
    Momento        = momento,
    RAM_R_MB       = mem$ram_r_mb,
    RAM_Sistema_MB = ifelse(is.na(mem$ram_sistema_mb), NA, mem$ram_sistema_mb),
    RAM_Total_MB   = ifelse(is.na(mem$ram_total_mb), NA, mem$ram_total_mb),
    Pct_Uso        = ifelse(is.na(mem$pct_uso), NA, mem$pct_uso),
    Objetos_Pesados = top5_txt,
    stringsAsFactors = FALSE
  )
  
  if (imprimir) {
    sistema_txt <- if (!is.na(mem$ram_sistema_mb)) {
      sprintf(" | Sistema: %.0f / %.0f MB (%.1f%%)",
              mem$ram_sistema_mb, mem$ram_total_mb, mem$pct_uso)
    } else {
      ""
    }
    cat(sprintf("  [RAM] %-8s | R: %7.1f MB%s\n",
                momento, mem$ram_r_mb, sistema_txt))
  }
  
  if (is.null(log_df)) {
    return(nueva_fila)
  }
  rbind(log_df, nueva_fila)
}