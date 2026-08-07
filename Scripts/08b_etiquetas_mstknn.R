# =============================================================================
# 08b_etiquetas_mstknn.R
# ASIGNACIÓN DE SÍMBOLOS HGNC A LOS CLUSTERS MST-kNN ITERATIVOS
# Reutiliza el caché HGNC existente — no requiere nueva consulta a la API
# =============================================================================
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# 1. CARGAR ASIGNACIONES MST-kNN ITERATIVAS
# =============================================================================
cat("=== CARGANDO ASIGNACIONES MST-kNN ITERATIVAS ===\n")

ruta_cache_asig <- file.path(DIR_CACHE, "mstknn_iter_asignaciones.rds")
ruta_excel_iter <- file.path(DIR_RESULTS, "mstknn_subdivision_iterativa.xlsx")

if (file.exists(ruta_cache_asig)) {
  cat("Cargando desde caché .rds...\n")
  asig_df <- readRDS(ruta_cache_asig)
} else if (file.exists(ruta_excel_iter)) {
  cat("Cargando desde Excel...\n")
  asig_df <- read.xlsx(ruta_excel_iter, sheet = "Asignaciones_Finales",
                       startRow = 2, colNames = TRUE)
} else {
  stop("No se encontraron asignaciones MST-kNN iterativas. ",
       "Ejecuta primero 07_iterar_mstknn.R")
}

asig_df$Arbol         <- as.character(asig_df$Arbol)
asig_df$Cluster_Final <- as.integer(asig_df$Cluster_Final)

n_clusters <- length(unique(asig_df$Cluster_Final))
cat(sprintf("Árboles cargados: %d | Clusters: %d\n", nrow(asig_df), n_clusters))

# =============================================================================
# 2. EXTRAER NCBI GENE ID DESDE EL NOMBRE DEL ÁRBOL
# =============================================================================
cat("\n=== EXTRAYENDO GENE IDs ===\n")

asig_df$Gene_ID <- sub("^(\\d+)_.*$", "\\1", asig_df$Arbol)

n_validos <- sum(grepl("^\\d+$", asig_df$Gene_ID))
cat(sprintf("Gene IDs numéricos extraídos: %d de %d\n", n_validos, nrow(asig_df)))

# =============================================================================
# 3. CARGAR CONVERSIÓN HGNC DESDE CACHÉ (no volver a consultar la API)
# =============================================================================
cat("\n=== CARGANDO CONVERSIÓN HGNC DESDE CACHÉ ===\n")

# Intentar ambos cachés disponibles
ruta_hgnc1 <- file.path(DIR_CACHE, "hgnc_conversion.rds")
ruta_hgnc2 <- file.path(DIR_CACHE, "hgnc_conversion2.rds")

if (file.exists(ruta_hgnc1)) {
  hgnc_df <- readRDS(ruta_hgnc1)
  cat("Caché HGNC cargado:", ruta_hgnc1, "\n")
} else if (file.exists(ruta_hgnc2)) {
  hgnc_df <- readRDS(ruta_hgnc2)
  cat("Caché HGNC cargado:", ruta_hgnc2, "\n")
} else {
  stop("No se encontró caché HGNC. Ejecuta primero 08_etiquetas.R")
}

hgnc_df$Gene_ID <- as.character(hgnc_df$Gene_ID)
n_conv <- sum(!is.na(hgnc_df$HGNC_Symbol))
cat(sprintf("Conversiones HGNC disponibles: %d (%.1f%%)\n",
            n_conv, 100 * n_conv / nrow(hgnc_df)))

# =============================================================================
# 4. UNIR ASIGNACIONES CON SÍMBOLOS HGNC
# =============================================================================
cat("\n=== UNIENDO DATOS ===\n")

genes_df <- merge(asig_df, hgnc_df[, c("Gene_ID", "HGNC_Symbol", "Gene_Name")],
                  by = "Gene_ID", all.x = TRUE)

# Reordenar columnas
genes_df <- genes_df[, c("Arbol", "Gene_ID", "HGNC_Symbol", "Gene_Name",
                          "Cluster_Final")]
genes_df <- genes_df[order(genes_df$Cluster_Final, genes_df$Arbol), ]
rownames(genes_df) <- NULL

n_con_hgnc <- sum(!is.na(genes_df$HGNC_Symbol))
cat(sprintf("Genes con símbolo HGNC: %d de %d (%.1f%%)\n",
            n_con_hgnc, nrow(genes_df),
            100 * n_con_hgnc / nrow(genes_df)))

# =============================================================================
# 5. EXPORTAR LISTAS HGNC POR CLUSTER
# =============================================================================
cat("\n=== EXPORTANDO LISTAS POR CLUSTER ===\n")

clusters_unicos <- sort(unique(genes_df$Cluster_Final))
wb_listas       <- createWorkbook()
resumen_listas  <- data.frame()

for (cid in clusters_unicos) {

  mask_cid <- !is.na(genes_df$Cluster_Final) & genes_df$Cluster_Final == cid
  simbolos_cluster <- genes_df$HGNC_Symbol[mask_cid & !is.na(genes_df$HGNC_Symbol)]

  n_arboles <- sum(mask_cid)
  n_hgnc    <- length(simbolos_cluster)

  resumen_listas <- rbind(resumen_listas, data.frame(
    Cluster          = cid,
    N_Arboles        = n_arboles,
    N_Con_HGNC       = n_hgnc,
    Pct_HGNC         = round(100 * n_hgnc / max(n_arboles, 1), 1),
    Apto_Enriquecimiento = n_hgnc >= 15,
    stringsAsFactors = FALSE
  ))

  if (n_hgnc == 0) next

  lista_df    <- data.frame(HGNC_Symbol = simbolos_cluster, stringsAsFactors = FALSE)
  nombre_hoja <- paste0("C", cid)
  if (nchar(nombre_hoja) > 31) nombre_hoja <- substr(nombre_hoja, 1, 31)

  wb_listas <- agregar_hoja_formateada(
    wb           = wb_listas,
    nombre_hoja  = nombre_hoja,
    titulo_tabla = paste0("Cluster ", cid, " (n=", n_hgnc, " genes HGNC)"),
    datos        = lista_df,
    anchos_col   = 20
  )
}

# Hoja de resumen al inicio
wb_listas <- agregar_hoja_formateada(
  wb           = wb_listas,
  nombre_hoja  = "Resumen",
  titulo_tabla = paste0("Resumen de genes HGNC por cluster MST-kNN iterativo (",
                        n_clusters, " clusters)"),
  datos        = resumen_listas,
  anchos_col   = "auto"
)

ruta_listas <- file.path(DIR_RESULTS,
                         paste0("hgnc_listas_mstknn_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb_listas, ruta_listas, overwrite = TRUE)
cat("Listas guardadas en:", ruta_listas, "\n")

# También guardar tabla completa genes_df como referencia
ruta_genes <- file.path(DIR_CACHE, "mstknn_iter_genes_hgnc.rds")
saveRDS(genes_df, ruta_genes)
cat("Tabla genes+HGNC guardada en:", ruta_genes, "\n")

cat(sprintf("\nClusters aptos (n >= 15): %d de %d\n",
            sum(resumen_listas$Apto_Enriquecimiento), n_clusters))
cat("\n=== COMPLETADO: 08b_etiquetas_mstknn.R ===\n")
