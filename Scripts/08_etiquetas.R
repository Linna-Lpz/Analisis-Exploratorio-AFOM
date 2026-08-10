# =============================================================================
# CONVERSIÓN NCBI GENE ID → HGNC APPROVED SYMBOL
# y exportación de hoja "Genes_Agrupados" al Excel de coordenadas PCA
# =============================================================================
library(httr)
library(jsonlite)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# PARÁMETROS
# =============================================================================
TAMANO_LOTE   <- 1000   # mygene.info acepta hasta 1000 IDs por request
ESPECIE       <- 9606   # 9606 = Homo sapiens (taxid NCBI)
# alternativa: 9443 = Primates (orden, no especie)
# Si quieres sin filtro de especie: ESPECIE <- NULL

# =============================================================================
# LEER ASIGNACIONES DEL EXCEL DE COORDENADAS PCA
#    (Se actualiza primero con las subdivisiones iterativas)
# =============================================================================
cat("=== CARGANDO Y ACTUALIZANDO ASIGNACIONES ===\n")

ruta_pca_excel <- file.path(DIR_RESULTS,
                            paste0("pca_coordenadas_", NOMBRE_BDD, ".xlsx"))

if (!file.exists(ruta_pca_excel)) {
  stop("No se encontró el Excel PCA: ", ruta_pca_excel,
       "\nEjecuta primero el script 05d_reduccion_dimensional_PCA.R")
}

# Leer coordenadas (que pueden no tener las columnas iterativas si 05d corrió antes)
wb_pca       <- loadWorkbook(ruta_pca_excel)
coords_full  <- read.xlsx(wb_pca,
                          sheet    = "Coordenadas_PCA",
                          startRow = 2,
                          colNames = TRUE)

# 1. Eliminar columnas de cluster antiguas para evitar duplicados (.x, .y)
cols_no_cluster <- grep("^Cluster_", colnames(coords_full), invert = TRUE, value = TRUE)
coords_base <- coords_full[, cols_no_cluster, drop = FALSE]

# 2. Refrescar todas las asignaciones desde los archivos de resultados
res_unidos <- unir_etiquetas_clustering(coords_base, DIR_RESULTS)
coords_full <- res_unidos$coords_df

# 3. Guardar las coordenadas actualizadas de vuelta en el Excel
removeWorksheet(wb_pca, "Coordenadas_PCA")
wb_pca <- agregar_hoja_formateada(
  wb           = wb_pca,
  nombre_hoja  = "Coordenadas_PCA",
  titulo_tabla = paste0("Coordenadas PCA — ", NOMBRE_BDD),
  datos        = coords_full,
  anchos_col   = "auto"
)
# Reordenar para que Coordenadas_PCA quede de primera
hojas_actuales <- names(wb_pca)
idx_coordenadas <- which(hojas_actuales == "Coordenadas_PCA")
worksheetOrder(wb_pca) <- c(idx_coordenadas, setdiff(seq_along(hojas_actuales), idx_coordenadas))
saveWorkbook(wb_pca, ruta_pca_excel, overwrite = TRUE)

cat("Excel de coordenadas PCA actualizado con asignaciones recientes.\n")

# Extraer columnas necesarias para HGNC
cols_cluster <- grep("^Cluster_", colnames(coords_full), value = TRUE)
cols_usar    <- c("Arbol", cols_cluster)

asig_df <- coords_full[, cols_usar, drop = FALSE]
asig_df$Arbol <- as.character(asig_df$Arbol)

cat(sprintf("Árboles cargados   : %d\n", nrow(asig_df)))
cat(sprintf("Columnas cluster   : %s\n", paste(cols_cluster, collapse = ", ")))

# =============================================================================
# EXTRAER NCBI GENE ID DESDE EL NOMBRE DEL ÁRBOL
# =============================================================================
cat("\n=== EXTRAYENDO GENE IDs ===\n")

# Patrón: número al inicio del nombre, antes del primer "_"
asig_df$Gene_ID <- sub("^(\\d+)_.*$", "\\1", asig_df$Arbol)

# Verificar extracción
n_validos <- sum(grepl("^\\d+$", asig_df$Gene_ID))
cat(sprintf("Gene IDs extraídos : %d de %d\n", n_validos, nrow(asig_df)))

ids_unicos <- unique(asig_df$Gene_ID[grepl("^\\d+$", asig_df$Gene_ID)])
cat(sprintf("Gene IDs únicos    : %d\n", length(ids_unicos)))

# =============================================================================
# 3. FUNCIÓN DE CONSULTA A MYGENE.INFO (en lotes)
# =============================================================================
consultar_mygene <- function(gene_ids, especie = 9606, tamano_lote = 1000) {
  
  n_total   <- length(gene_ids)
  n_lotes   <- ceiling(n_total / tamano_lote)
  resultado <- data.frame()
  
  cat(sprintf("Consultando mygene.info: %d IDs en %d lote(s)...\n",
              n_total, n_lotes))
  
  for (i in seq_len(n_lotes)) {
    
    idx_inicio <- (i - 1) * tamano_lote + 1
    idx_fin    <- min(i * tamano_lote, n_total)
    lote       <- gene_ids[idx_inicio:idx_fin]
    
    cat(sprintf("  Lote %d/%d (IDs %d a %d)...\n",
                i, n_lotes, idx_inicio, idx_fin))
    
    # Construir body del POST
    body_params <- list(
      ids    = paste(lote, collapse = ","),
      fields = "symbol,name,taxid,other_names",
      species = if (!is.null(especie)) as.character(especie) else "all"
    )
    
    respuesta <- tryCatch(
      POST("https://mygene.info/v3/gene",
           body    = body_params,
           encode  = "form",
           timeout(60)),
      error = function(e) {
        cat(sprintf("    ERROR en lote %d: %s\n", i, e$message))
        NULL
      }
    )
    
    if (is.null(respuesta) || status_code(respuesta) != 200) {
      cat(sprintf("    AVISO: lote %d falló (status %s). Se omite.\n",
                  i, status_code(respuesta)))
      next
    }
    
    contenido <- fromJSON(rawToChar(respuesta$content),
                          simplifyDataFrame = TRUE)
    
    # Extraer symbol para cada ID del lote
    lote_df <- data.frame(
      Gene_ID      = as.character(lote),
      HGNC_Symbol  = NA_character_,
      Gene_Name    = NA_character_,
      Tax_ID       = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (is.data.frame(contenido) && nrow(contenido) > 0) {
      
      for (j in seq_len(nrow(contenido))) {
        
        fila      <- contenido[j, ]
        query_id  <- as.character(fila$query)
        
        # Saltar si mygene no encontró el ID
        if (!is.null(fila$notfound) && isTRUE(fila$notfound)) next
        
        symbol <- if (!is.null(fila$symbol) && !is.na(fila$symbol))
          as.character(fila$symbol) else NA_character_
        name   <- if (!is.null(fila$name)   && !is.na(fila$name))
          as.character(fila$name)   else NA_character_
        taxid  <- if (!is.null(fila$taxid)  && !is.na(fila$taxid))
          as.character(fila$taxid)  else NA_character_
        
        idx <- which(lote_df$Gene_ID == query_id)
        if (length(idx) > 0) {
          lote_df$HGNC_Symbol[idx] <- symbol
          lote_df$Gene_Name[idx]   <- name
          lote_df$Tax_ID[idx]      <- taxid
        }
      }
    }
    
    resultado <- rbind(resultado, lote_df)
    Sys.sleep(0.5)  # pausa entre lotes para no saturar la API
  }
  
  return(resultado)
}

# =============================================================================
# 4. EJECUTAR CONSULTA — con caché para no repetir si ya existe
# =============================================================================
cat("\n=== CONSULTANDO MYGENE.INFO ===\n")

ruta_cache_hgnc <- file.path(DIR_CACHE, "hgnc_conversion.rds")

if (file.exists(ruta_cache_hgnc)) {
  cat("Conversión HGNC encontrada en caché. Cargando...\n")
  hgnc_df <- readRDS(ruta_cache_hgnc)
} else {
  cat("Consultando API (puede tardar varios minutos)...\n")
  tiempo_api <- proc.time()
  
  hgnc_df <- consultar_mygene(
    gene_ids    = ids_unicos,
    especie     = ESPECIE,
    tamano_lote = TAMANO_LOTE
  )
  
  tiempo_api <- as.numeric(proc.time() - tiempo_api)[3]
  saveRDS(hgnc_df, ruta_cache_hgnc)
  cat(sprintf("Consulta completada en %.1f segundos. Caché guardado.\n", tiempo_api))
}

# Diagnóstico de conversión
n_convertidos  <- sum(!is.na(hgnc_df$HGNC_Symbol))
n_sin_simbolo  <- sum(is.na(hgnc_df$HGNC_Symbol))

cat(sprintf("\nResultados conversión:\n"))
cat(sprintf("  Con símbolo HGNC    : %d (%.1f%%)\n",
            n_convertidos, 100 * n_convertidos / nrow(hgnc_df)))
cat(sprintf("  Sin símbolo (NA)    : %d (%.1f%%)\n",
            n_sin_simbolo, 100 * n_sin_simbolo / nrow(hgnc_df)))

# =============================================================================
# 5. UNIR CONVERSIÓN CON ASIGNACIONES
# =============================================================================
cat("\n=== UNIENDO DATOS ===\n")

genes_df <- merge(asig_df, hgnc_df,
                  by = "Gene_ID", all.x = TRUE)

# Reordenar columnas: Arbol | Gene_ID | HGNC_Symbol | Gene_Name | Tax_ID | Clusters
cols_orden <- c("Arbol", "Gene_ID", "HGNC_Symbol", "Gene_Name", "Tax_ID", cols_cluster)
genes_df   <- genes_df[, cols_orden]
genes_df   <- genes_df[order(genes_df$Arbol), ]
rownames(genes_df) <- NULL

cat(sprintf("Tabla final: %d filas x %d columnas\n",
            nrow(genes_df), ncol(genes_df)))
print(head(genes_df[, 1:5], 5))

# =============================================================================
# 6. AGREGAR HOJA AL EXCEL PCA EXISTENTE
# =============================================================================
cat("\n=== AGREGANDO HOJA AL EXCEL ===\n")

wb <- loadWorkbook(ruta_pca_excel)

# Eliminar hoja si ya existe (para poder sobreescribir)
if ("Genes_Agrupados" %in% getSheetNames(ruta_pca_excel)) {
  removeWorksheet(wb, "Genes_Agrupados")
  cat("Hoja 'Genes_Agrupados' existente reemplazada.\n")
}

wb <- agregar_hoja_formateada(
  wb           = wb,
  nombre_hoja  = "Genes_Agrupados",
  titulo_tabla = paste0("Genes Agrupados — NCBI Gene ID → HGNC Approved Symbol (", NOMBRE_BDD, ")"),
  datos        = genes_df,
  anchos_col   = "auto"
)

saveWorkbook(wb, ruta_pca_excel, overwrite = TRUE)
cat("Hoja 'Genes_Agrupados' agregada a:", ruta_pca_excel, "\n")

# =============================================================================
# 7. EXPORTAR TAMBIÉN LISTA PLANA DE SÍMBOLOS HGNC POR CLUSTER
#    Lista lista para pegar en http://bioinformatics.sdstate.edu/go/
# =============================================================================
cat("\n=== EXPORTANDO LISTAS POR CLUSTER ===\n")

algoritmo <- if(exists("ALGORITMO_DOWNSTREAM")) ALGORITMO_DOWNSTREAM else "AUTO"

col_lista <- if (algoritmo == "MST-kNN" && "Cluster_MSTKNN_Iter" %in% cols_cluster) {
  "Cluster_MSTKNN_Iter"
} else if (algoritmo == "K-Means" && "Cluster_KMeans_Iter" %in% cols_cluster) {
  "Cluster_KMeans_Iter"
} else if (algoritmo == "CLARA" && "Cluster_CLARA_Iter" %in% cols_cluster) {
  "Cluster_CLARA_Iter"
} else if (algoritmo == "PAM" && "Cluster_PAM" %in% cols_cluster) {
  "Cluster_PAM"
} else {
  # Fallback si está en AUTO o no se encuentra la columna iterada
  if ("Cluster_KMeans_Iter" %in% cols_cluster) "Cluster_KMeans_Iter" else
    if ("Cluster_MSTKNN_Iter" %in% cols_cluster) "Cluster_MSTKNN_Iter" else
      cols_cluster[1]
}

cat(sprintf("Algoritmo objetivo: %s\n", algoritmo))
cat(sprintf("Columna de clustering para listas: %s\n", col_lista))

clusters_unicos <- sort(unique(na.omit(genes_df[[col_lista]])))

# Un Excel separado con una hoja por cluster — solo símbolos HGNC válidos
wb_listas <- createWorkbook()

resumen_listas <- data.frame()

for (cid in clusters_unicos) {
  
  simbolos_cluster <- genes_df$HGNC_Symbol[
    !is.na(genes_df[[col_lista]]) & genes_df[[col_lista]] == cid &
      !is.na(genes_df$HGNC_Symbol)
  ]
  
  if (length(simbolos_cluster) == 0) next
  
  lista_df <- data.frame(HGNC_Symbol = simbolos_cluster,
                         stringsAsFactors = FALSE)
  
  nombre_hoja <- paste0("Cluster_", cid)
  # openxlsx limita nombres de hoja a 31 caracteres
  if (nchar(nombre_hoja) > 31) nombre_hoja <- substr(nombre_hoja, 1, 31)
  
  wb_listas <- agregar_hoja_formateada(
    wb           = wb_listas,
    nombre_hoja  = nombre_hoja,
    titulo_tabla = paste0("Genes Cluster ", cid,
                          " (n=", length(simbolos_cluster), ") — HGNC Symbols"),
    datos        = lista_df,
    anchos_col   = 20
  )
  
  resumen_listas <- rbind(resumen_listas, data.frame(
    Cluster          = cid,
    N_Arboles        = sum(!is.na(genes_df[[col_lista]]) &
                             genes_df[[col_lista]] == cid),
    N_Con_HGNC       = length(simbolos_cluster),
    Pct_Con_HGNC     = round(100 * length(simbolos_cluster) /
                               sum(!is.na(genes_df[[col_lista]]) &
                                     genes_df[[col_lista]] == cid), 1)
  ))
}

# Hoja de resumen al inicio
wb_listas <- agregar_hoja_formateada(
  wb           = wb_listas,
  nombre_hoja  = "Resumen",
  titulo_tabla = paste0("Resumen de genes HGNC por cluster — ", col_lista),
  datos        = resumen_listas,
  anchos_col   = "auto"
)

ruta_listas <- file.path(DIR_RESULTS,
                         paste0("hgnc_listas_", gsub("-","", tolower(algoritmo)), "_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb_listas, ruta_listas, overwrite = TRUE)

cat("Listas por cluster guardadas en:", ruta_listas, "\n")
cat(sprintf("\nResumen:\n"))
print(resumen_listas)
cat("\n=== COMPLETADO ===\n")