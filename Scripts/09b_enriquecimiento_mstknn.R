# =============================================================================
# 09b_enriquecimiento_mstknn.R
# ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL — gprofiler2 para clusters MST-kNN
# =============================================================================
library(gprofiler2)
library(openxlsx)
library(ggplot2)
library(patchwork)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# PARÁMETROS
# =============================================================================
ORGANISMO         <- "hsapiens"
FUENTES           <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")
P_VALOR_UMBRAL    <- 0.05
MIN_GENES_TERMINO <- 3
MAX_TERMINOS_PLOT <- 15
MIN_GENES_CLUSTER <- 15  # < 15 se omiten por falta de potencia estadística

# =============================================================================
# 1. CARGAR TABLA DE GENES (MST-kNN)
# =============================================================================
cat("=== CARGANDO GENES MST-kNN ===\n")

ruta_genes <- file.path(DIR_CACHE, "mstknn_iter_genes_hgnc.rds")
if (!file.exists(ruta_genes)) {
  stop("No se encontró ", ruta_genes, ". Ejecuta 08b_etiquetas_mstknn.R primero.")
}

genes_df <- readRDS(ruta_genes)

# Solo clusters con n >= MIN_GENES_CLUSTER
tamanos <- table(genes_df$Cluster_Final[!is.na(genes_df$HGNC_Symbol)])
clusters_validos <- as.integer(names(tamanos[tamanos >= MIN_GENES_CLUSTER]))
clusters_unicos <- sort(unique(na.omit(genes_df$Cluster_Final)))

cat(sprintf("Clusters totales        : %d\n", length(clusters_unicos)))
cat(sprintf("Clusters aptos (n >= %d): %d\n", MIN_GENES_CLUSTER, length(clusters_validos)))
cat(sprintf("Fondo (background)      : %d genes HGNC únicos\n",
            length(unique(na.omit(genes_df$HGNC_Symbol)))))

# =============================================================================
# 2. FUNCIÓN DE ENRIQUECIMIENTO
# =============================================================================
enriquecer_cluster <- function(simbolos, cluster_id, organismo,
                               fuentes, p_umbral, min_genes, bg_universe) {

  simbolos_validos <- simbolos[!is.na(simbolos) & simbolos != ""]

  cat(sprintf("    Cluster %s: %d genes validos → consultando...\n",
              cluster_id, length(simbolos_validos)))

  resultado <- tryCatch(
    gost(query              = simbolos_validos,
         organism           = organismo,
         sources            = fuentes,
         correction_method  = "fdr",
         user_threshold     = p_umbral,
         significant        = TRUE,
         measure_underrepresentation = FALSE,
         evcodes            = TRUE,
         custom_bg          = bg_universe),
    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
  )

  if (is.null(resultado) || is.null(resultado$result) || nrow(resultado$result) == 0) {
    return(NULL)
  }

  res_df <- resultado$result
  res_df <- res_df[res_df$intersection_size >= min_genes, ]
  if (nrow(res_df) == 0) return(NULL)

  res_df$Cluster <- cluster_id
  res_df$N_Genes_Cluster <- length(simbolos_validos)

  cols_exportar <- c("Cluster", "N_Genes_Cluster", "source", "term_id",
                     "term_name", "p_value", "intersection_size",
                     "term_size", "query_size", "intersection")
  cols_exportar <- cols_exportar[cols_exportar %in% colnames(res_df)]
  res_df <- res_df[, cols_exportar]

  colnames(res_df)[colnames(res_df) == "source"]            <- "Fuente"
  colnames(res_df)[colnames(res_df) == "term_id"]           <- "Term_ID"
  colnames(res_df)[colnames(res_df) == "term_name"]         <- "Termino"
  colnames(res_df)[colnames(res_df) == "p_value"]           <- "FDR"
  colnames(res_df)[colnames(res_df) == "intersection_size"] <- "N_Genes_en_Termino"
  colnames(res_df)[colnames(res_df) == "term_size"]         <- "Tamano_Termino"
  colnames(res_df)[colnames(res_df) == "query_size"]        <- "Tamano_Query"
  colnames(res_df)[colnames(res_df) == "intersection"]      <- "Genes_en_Termino"

  res_df <- res_df[order(res_df$FDR), ]
  rownames(res_df) <- NULL
  return(res_df)
}

# =============================================================================
# 3. EJECUTAR ENRIQUECIMIENTO (solo para clusters válidos)
# =============================================================================
cat("\n=== EJECUTANDO ENRIQUECIMIENTO FUNCIONAL ===\n")

universo_hgnc <- unique(na.omit(genes_df$HGNC_Symbol))
universo_hgnc <- universo_hgnc[universo_hgnc != ""]

dir_cache_enrich <- file.path(DIR_CACHE, "enrichment_mstknn")
if (!dir.exists(dir_cache_enrich)) dir.create(dir_cache_enrich, recursive = TRUE)

resultados_todos <- list()
resumen_enrich   <- data.frame()
tiempo_total     <- proc.time()

for (cid in clusters_unicos) {

  simbolos_cid <- genes_df$HGNC_Symbol[!is.na(genes_df$Cluster_Final) &
                                         genes_df$Cluster_Final == cid]

  if (!(cid %in% clusters_validos)) {
    # No analizar
    res_cid <- NULL
  } else {
    ruta_cache_cid <- file.path(dir_cache_enrich, paste0("enrich_cluster_", cid, ".rds"))

    if (file.exists(ruta_cache_cid)) {
      res_cid <- readRDS(ruta_cache_cid)
      cat(sprintf("  Cluster %s: cargado desde caché.\n", cid))
    } else {
      res_cid <- enriquecer_cluster(simbolos_cid, cid, ORGANISMO, FUENTES,
                                    P_VALOR_UMBRAL, MIN_GENES_TERMINO, universo_hgnc)
      saveRDS(res_cid, ruta_cache_cid)
      Sys.sleep(0.3)
    }
  }

  resultados_todos[[as.character(cid)]] <- res_cid

  resumen_enrich <- rbind(resumen_enrich, data.frame(
    Cluster        = cid,
    N_Genes_Input  = length(simbolos_cid),
    N_Genes_HGNC   = sum(!is.na(simbolos_cid)),
    Apto           = cid %in% clusters_validos,
    N_Terminos_Sig = if (!is.null(res_cid)) nrow(res_cid) else 0,
    N_GO_BP        = if (!is.null(res_cid)) sum(res_cid$Fuente == "GO:BP") else 0,
    N_GO_MF        = if (!is.null(res_cid)) sum(res_cid$Fuente == "GO:MF") else 0,
    N_GO_CC        = if (!is.null(res_cid)) sum(res_cid$Fuente == "GO:CC") else 0,
    N_KEGG         = if (!is.null(res_cid)) sum(res_cid$Fuente == "KEGG")  else 0,
    N_Reactome     = if (!is.null(res_cid)) sum(res_cid$Fuente == "REAC")  else 0,
    stringsAsFactors = FALSE
  ))
}

tiempo_total <- as.numeric(proc.time() - tiempo_total)[3]
resultados_df <- do.call(rbind, Filter(Negate(is.null), resultados_todos))
rownames(resultados_df) <- NULL

cat(sprintf("Clusters con términos significativos: %d\n", sum(resumen_enrich$N_Terminos_Sig > 0)))

# =============================================================================
# 4. GRÁFICOS
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS ===\n")
dir_graficos <- file.path(DIR_RESULTS, "enrichment_plots_mstknn")
if (!dir.exists(dir_graficos)) dir.create(dir_graficos, recursive = TRUE)

paleta_fuentes <- c("GO:BP" = "#4DAF4A", "GO:MF" = "#377EB8", "GO:CC" = "#984EA3",
                    "KEGG"  = "#FF7F00", "REAC"  = "#E41A1C")

clusters_con_resultados <- names(Filter(Negate(is.null), resultados_todos))

for (cid in clusters_con_resultados) {
  res_plot <- resultados_todos[[cid]]
  if (is.null(res_plot) || nrow(res_plot) == 0) next

  res_plot <- head(res_plot[order(res_plot$FDR), ], MAX_TERMINOS_PLOT)
  res_plot$Termino_corto <- ifelse(nchar(res_plot$Termino) > 50,
                                   paste0(substr(res_plot$Termino, 1, 47), "..."),
                                   res_plot$Termino)
  duplicados <- duplicated(res_plot$Termino_corto) | duplicated(res_plot$Termino_corto, fromLast = TRUE)
  res_plot$Termino_corto[duplicados] <- paste0(res_plot$Termino_corto[duplicados], " [", res_plot$Term_ID[duplicados], "]")

  res_plot$Termino_corto <- factor(res_plot$Termino_corto, levels = rev(unique(res_plot$Termino_corto)))
  res_plot$log10_FDR <- -log10(res_plot$FDR)

  p_enrich <- ggplot(res_plot, aes(x = log10_FDR, y = Termino_corto, fill = Fuente, size = N_Genes_en_Termino)) +
    geom_point(shape = 21, alpha = 0.85) +
    scale_fill_manual(values = paleta_fuentes, name = "Fuente") +
    scale_size_continuous(name = "Genes en\ntérmino", range = c(3, 10)) +
    geom_vline(xintercept = -log10(P_VALOR_UMBRAL), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    labs(title = paste0("Enriquecimiento Funcional — Cluster ", cid, " (MST-kNN)"),
         x = expression(-log[10](FDR)), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(dir_graficos, paste0("enrichment_mstknn_cluster_", cid, ".png")),
         plot = p_enrich, width = 10, height = max(5, nrow(res_plot) * 0.4 + 2), dpi = 300)
}

# =============================================================================
# 5. EXPORTAR EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")
wb <- createWorkbook()
wb <- agregar_hoja_formateada(wb, "Resumen", "Resumen Enriquecimiento MST-kNN", resumen_enrich, "auto")
if (nrow(resultados_df) > 0) {
  wb <- agregar_hoja_formateada(wb, "Todos_los_Terminos", "Términos MST-kNN", resultados_df, "auto")
}
for (cid in clusters_con_resultados) {
  res_cid <- resultados_todos[[cid]]
  if (is.null(res_cid) || nrow(res_cid) == 0) next
  nombre_hoja <- paste0("C", cid)
  wb <- agregar_hoja_formateada(wb, nombre_hoja, paste0("Cluster ", cid), res_cid, "auto")
}

ruta_excel <- file.path(DIR_RESULTS, paste0("enrichment_funcional_mstknn_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Excel guardado:", ruta_excel, "\n")
cat("\n=== COMPLETADO: 09b_enriquecimiento_mstknn.R ===\n")
