# =============================================================================
# ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL — gprofiler2
# GO terms, KEGG, Reactome por cluster
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
ORGANISMO         <- "hsapiens"     # Homo sapiens
FUENTES           <- c("GO:BP",     # Gene Ontology — Biological Process
                       "GO:MF",    # Gene Ontology — Molecular Function
                       "GO:CC",    # Gene Ontology — Cellular Component
                       "KEGG",     # KEGG Pathways
                       "REAC")     # Reactome Pathways
P_VALOR_UMBRAL    <- 0.05           # umbral de significancia (FDR corregido)
MIN_GENES_TERMINO <- 3              # mínimo de genes del cluster en el término
MAX_TERMINOS_PLOT <- 15             # top N términos a graficar por cluster

# Columna de clustering a analizar
# Prioriza CLARA iterativa, luego CLARA óptimo
COL_CLUSTER <- "Cluster_KMeans_Iter" # otra op Cluster_Kmeans_Iter

# =============================================================================
# 1. CARGAR TABLA DE GENES AGRUPADOS
# =============================================================================
cat("=== CARGANDO GENES AGRUPADOS ===\n")

ruta_pca_excel <- file.path(DIR_RESULTS,
                            paste0("pca_coordenadas_", NOMBRE_BDD, ".xlsx"))

if (!file.exists(ruta_pca_excel)) {
  stop("No se encontró el Excel PCA: ", ruta_pca_excel)
}

hojas <- getSheetNames(ruta_pca_excel)
if (!"Genes_Agrupados" %in% hojas) {
  stop("No se encontró la hoja 'Genes_Agrupados'. ",
       "Ejecuta primero el script de conversión HGNC.")
}

genes_df <- read.xlsx(ruta_pca_excel,
                      sheet    = "Genes_Agrupados",
                      startRow = 2,
                      colNames = TRUE)

genes_df$Arbol       <- as.character(genes_df$Arbol)
genes_df$HGNC_Symbol <- as.character(genes_df$HGNC_Symbol)

# Verificar columna de clustering
if (!COL_CLUSTER %in% colnames(genes_df)) {
  cols_disponibles <- grep("^Cluster_", colnames(genes_df), value = TRUE)
  cat(sprintf("AVISO: '%s' no encontrada. Columnas disponibles: %s\n",
              COL_CLUSTER, paste(cols_disponibles, collapse = ", ")))
  COL_CLUSTER <- cols_disponibles[1]
  cat(sprintf("Usando: %s\n", COL_CLUSTER))
}

clusters_unicos <- sort(unique(na.omit(genes_df[[COL_CLUSTER]])))
cat(sprintf("Clusters a analizar : %d\n", length(clusters_unicos)))
cat(sprintf("Genes totales       : %d\n", nrow(genes_df)))
cat(sprintf("Con símbolo HGNC    : %d\n", sum(!is.na(genes_df$HGNC_Symbol))))

# =============================================================================
# 2. FUNCIÓN DE ENRIQUECIMIENTO POR CLUSTER
# =============================================================================
enriquecer_cluster <- function(simbolos, cluster_id, organismo,
                               fuentes, p_umbral, min_genes, bg_universe) {
  
  # Filtrar NAs
  simbolos_validos <- simbolos[!is.na(simbolos) & simbolos != ""]
  
  if (length(simbolos_validos) < 3) {
    cat(sprintf("    Cluster %s: menos de 3 genes válidos. Saltando.\n", cluster_id))
    return(NULL)
  }
  
  cat(sprintf("    Cluster %s: %d genes → consultando gprofiler2...\n",
              cluster_id, length(simbolos_validos)))
  
  resultado <- tryCatch(
    gost(query              = simbolos_validos,
         organism           = organismo,
         sources            = fuentes,
         correction_method  = "fdr",
         user_threshold     = p_umbral,
         significant        = TRUE,
         measure_underrepresentation = FALSE,
         evcodes            = TRUE,    # incluir lista de genes por término
         custom_bg          = bg_universe),
    error = function(e) {
      cat(sprintf("    ERROR en cluster %s: %s\n", cluster_id, e$message))
      NULL
    }
  )
  
  if (is.null(resultado) || is.null(resultado$result) ||
      nrow(resultado$result) == 0) {
    cat(sprintf("    Cluster %s: sin términos significativos.\n", cluster_id))
    return(NULL)
  }
  
  res_df <- resultado$result
  
  # Filtrar por mínimo de genes del cluster en el término
  res_df <- res_df[res_df$intersection_size >= min_genes, ]
  
  if (nrow(res_df) == 0) return(NULL)
  
  # Columnas relevantes + identificador de cluster
  res_df$Cluster         <- cluster_id
  res_df$N_Genes_Cluster <- length(simbolos_validos)
  
  # Seleccionar y renombrar columnas para claridad
  cols_exportar <- c("Cluster", "N_Genes_Cluster", "source", "term_id",
                     "term_name", "p_value",
                     "intersection_size", "term_size", "query_size",
                     "effective_domain_size", "intersection")
  
  cols_exportar <- cols_exportar[cols_exportar %in% colnames(res_df)]
  res_df        <- res_df[, cols_exportar]
  
  colnames(res_df)[colnames(res_df) == "source"]           <- "Fuente"
  colnames(res_df)[colnames(res_df) == "term_id"]          <- "Term_ID"
  colnames(res_df)[colnames(res_df) == "term_name"]        <- "Termino"
  colnames(res_df)[colnames(res_df) == "p_value"]          <- "FDR"
  colnames(res_df)[colnames(res_df) == "intersection_size"]<- "N_Genes_en_Termino"
  colnames(res_df)[colnames(res_df) == "term_size"]        <- "Tamano_Termino"
  colnames(res_df)[colnames(res_df) == "query_size"]       <- "Tamano_Query"
  colnames(res_df)[colnames(res_df) == "intersection"]     <- "Genes_en_Termino"
  
  res_df <- res_df[order(res_df$FDR), ]
  rownames(res_df) <- NULL
  
  return(res_df)
}

# =============================================================================
# 3. EJECUTAR ENRIQUECIMIENTO — con caché por cluster
# =============================================================================
cat("\n=== EJECUTANDO ENRIQUECIMIENTO FUNCIONAL ===\n")
cat(sprintf("Fuentes    : %s\n", paste(FUENTES, collapse = ", ")))
cat(sprintf("FDR umbral : %.2f\n", P_VALOR_UMBRAL))
cat(sprintf("Min genes  : %d\n", MIN_GENES_TERMINO))

universo_hgnc <- unique(na.omit(genes_df$HGNC_Symbol))
universo_hgnc <- universo_hgnc[universo_hgnc != ""]
cat(sprintf("Fondo (bg) : %d genes únicos\n\n", length(universo_hgnc)))

dir_cache_enrich <- file.path(DIR_CACHE, "enrichment")
if (!dir.exists(dir_cache_enrich)) dir.create(dir_cache_enrich, recursive = TRUE)

resultados_todos <- list()
resumen_enrich   <- data.frame()
tiempo_total     <- proc.time()

for (cid in clusters_unicos) {
  
  ruta_cache_cid <- file.path(dir_cache_enrich,
                              paste0("enrich_cluster_", cid, ".rds"))
  
  simbolos_cid <- genes_df$HGNC_Symbol[
    !is.na(genes_df[[COL_CLUSTER]]) & genes_df[[COL_CLUSTER]] == cid
  ]
  
  if (file.exists(ruta_cache_cid)) {
    cat(sprintf("  Cluster %s: cargando desde caché.\n", cid))
    res_cid <- readRDS(ruta_cache_cid)
  } else {
    res_cid <- enriquecer_cluster(
      simbolos    = simbolos_cid,
      cluster_id  = cid,
      organismo   = ORGANISMO,
      fuentes     = FUENTES,
      p_umbral    = P_VALOR_UMBRAL,
      min_genes   = MIN_GENES_TERMINO,
      bg_universe = universo_hgnc
    )
    saveRDS(res_cid, ruta_cache_cid)
    Sys.sleep(0.3)  # pausa para no saturar la API
  }
  
  resultados_todos[[as.character(cid)]] <- res_cid
  
  # Fila de resumen
  resumen_enrich <- rbind(resumen_enrich, data.frame(
    Cluster           = cid,
    N_Genes_Input     = length(simbolos_cid),
    N_Genes_HGNC      = sum(!is.na(simbolos_cid)),
    N_Terminos_Sig    = if (!is.null(res_cid)) nrow(res_cid) else 0,
    N_GO_BP           = if (!is.null(res_cid)) sum(res_cid$Fuente == "GO:BP") else 0,
    N_GO_MF           = if (!is.null(res_cid)) sum(res_cid$Fuente == "GO:MF") else 0,
    N_GO_CC           = if (!is.null(res_cid)) sum(res_cid$Fuente == "GO:CC") else 0,
    N_KEGG            = if (!is.null(res_cid)) sum(res_cid$Fuente == "KEGG")  else 0,
    N_Reactome        = if (!is.null(res_cid)) sum(res_cid$Fuente == "REAC")  else 0,
    Top_Termino       = if (!is.null(res_cid) && nrow(res_cid) > 0)
      res_cid$Termino[1] else "—",
    stringsAsFactors = FALSE
  ))
}

tiempo_total <- as.numeric(proc.time() - tiempo_total)[3]
cat(sprintf("\nEnriquecimiento completado en %.1f segundos.\n", tiempo_total))

# Tabla combinada de todos los clusters
resultados_df <- do.call(rbind, Filter(Negate(is.null), resultados_todos))
rownames(resultados_df) <- NULL

cat(sprintf("Clusters con términos significativos: %d de %d\n",
            sum(resumen_enrich$N_Terminos_Sig > 0), length(clusters_unicos)))
cat(sprintf("Total términos significativos        : %d\n", nrow(resultados_df)))

# =============================================================================
# 4. GRÁFICOS — top términos por cluster (solo clusters con resultados)
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS ===\n")

dir_graficos <- file.path(DIR_RESULTS, "enrichment_plots")
if (!dir.exists(dir_graficos)) dir.create(dir_graficos, recursive = TRUE)

# Paleta por fuente
paleta_fuentes <- c(
  "GO:BP" = "#4DAF4A",
  "GO:MF" = "#377EB8",
  "GO:CC" = "#984EA3",
  "KEGG"  = "#FF7F00",
  "REAC"  = "#E41A1C"
)

clusters_con_resultados <- names(Filter(Negate(is.null), resultados_todos))

for (cid in clusters_con_resultados) {
  
  res_plot <- resultados_todos[[cid]]
  if (is.null(res_plot) || nrow(res_plot) == 0) next
  
  # Top N términos por FDR
  res_plot <- res_plot[order(res_plot$FDR), ]
  res_plot <- head(res_plot, MAX_TERMINOS_PLOT)
  
  # Truncar nombres largos
  res_plot$Termino_corto <- ifelse(
    nchar(res_plot$Termino) > 50,
    paste0(substr(res_plot$Termino, 1, 47), "..."),
    res_plot$Termino
  )
  
  # Desambiguar duplicados añadiendo sufijo con Term_ID
  duplicados <- duplicated(res_plot$Termino_corto) |
    duplicated(res_plot$Termino_corto, fromLast = TRUE)
  
  res_plot$Termino_corto[duplicados] <- paste0(
    res_plot$Termino_corto[duplicados], " [", res_plot$Term_ID[duplicados], "]"
  )
  
  # Truncar de nuevo si el sufijo alargó demasiado
  res_plot$Termino_corto <- ifelse(
    nchar(res_plot$Termino_corto) > 65,
    paste0(substr(res_plot$Termino_corto, 1, 55), "...[",
           res_plot$Term_ID, "]"),
    res_plot$Termino_corto
  )
  
  res_plot$Termino_corto <- factor(res_plot$Termino_corto,
                                   levels = rev(unique(res_plot$Termino_corto)))
  res_plot$log10_FDR     <- -log10(res_plot$FDR)
  
  p_enrich <- ggplot(res_plot,
                     aes(x = log10_FDR, y = Termino_corto,
                         fill = Fuente, size = N_Genes_en_Termino)) +
    geom_point(shape = 21, alpha = 0.85) +
    scale_fill_manual(values = paleta_fuentes, name = "Fuente") +
    scale_size_continuous(name = "Genes en\ntérmino", range = c(3, 10)) +
    geom_vline(xintercept = -log10(P_VALOR_UMBRAL),
               linetype = "dashed", color = "gray50", linewidth = 0.5) +
    labs(
      title    = paste0("Enriquecimiento Funcional — Cluster ", cid),
      subtitle = sprintf("Top %d términos  |  FDR < %.2f  |  %d genes input",
                         nrow(res_plot), P_VALOR_UMBRAL,
                         resumen_enrich$N_Genes_HGNC[resumen_enrich$Cluster == cid]),
      x        = expression(-log[10](FDR)),
      y        = NULL,
      caption  = paste0("Fuentes: ", paste(FUENTES, collapse = ", "))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "gray40", size = 9),
      plot.caption  = element_text(color = "gray60", size = 8),
      legend.position = "right",
      axis.text.y   = element_text(size = 8)
    )
  
  ruta_plot <- file.path(dir_graficos,
                         paste0("enrichment_cluster_", cid, "_",
                                NOMBRE_BDD, ".png"))
  ggsave(ruta_plot, plot = p_enrich,
         width = 10, height = max(5, nrow(res_plot) * 0.4 + 2),
         dpi = 300)
}

cat(sprintf("Gráficos guardados en: %s\n", dir_graficos))

# =============================================================================
# 5. EXPORTAR EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

wb <- createWorkbook()

# Hoja 1: Resumen por cluster
wb <- agregar_hoja_formateada(wb, "Resumen",
                              paste0("Resumen de Enriquecimiento por Cluster — ",
                                     NOMBRE_BDD),
                              resumen_enrich,
                              anchos_col = "auto")

# Hoja 2: Todos los términos combinados
wb <- agregar_hoja_formateada(wb, "Todos_los_Terminos",
                              paste0("Términos Significativos — Todos los Clusters",
                                     " (FDR < ", P_VALOR_UMBRAL, ")"),
                              resultados_df,
                              anchos_col = "auto")

# Hojas individuales por cluster (solo los que tienen resultados)
for (cid in clusters_con_resultados) {
  res_cid <- resultados_todos[[cid]]
  if (is.null(res_cid) || nrow(res_cid) == 0) next
  
  nombre_hoja <- paste0("Cluster_", cid)
  if (nchar(nombre_hoja) > 31) nombre_hoja <- substr(nombre_hoja, 1, 31)
  
  wb <- agregar_hoja_formateada(
    wb           = wb,
    nombre_hoja  = nombre_hoja,
    titulo_tabla = paste0("Enriquecimiento Cluster ", cid,
                          " — ", nrow(res_cid), " términos significativos"),
    datos        = res_cid,
    anchos_col   = "auto"
  )
}

# Hoja de parámetros
parametros_df <- data.frame(
  Parametro = c("organismo", "fuentes", "fdr_umbral", "min_genes_termino",
                "col_cluster", "n_clusters_analizados",
                "n_clusters_con_resultados", "tiempo_total_s"),
  Valor     = c(ORGANISMO, paste(FUENTES, collapse = ", "),
                P_VALOR_UMBRAL, MIN_GENES_TERMINO,
                COL_CLUSTER, length(clusters_unicos),
                length(clusters_con_resultados),
                round(tiempo_total, 1))
)

wb <- agregar_hoja_formateada(wb, "Parametros",
                              "Parámetros del Análisis",
                              parametros_df,
                              anchos_col = c(30, 40))

ruta_excel <- file.path(DIR_RESULTS,
                        paste0("enrichment_funcional_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Excel guardado:", ruta_excel, "\n")
cat("\n=== COMPLETADO ===\n")