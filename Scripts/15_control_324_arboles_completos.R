# =============================================================================
# 15_control_324_arboles_completos.R
#
# PROPOSITO (Observacion 1.2 del informe de tesis):
#   Los 324 arboles con las 190 especies son un grupo de control exacto:
#   no requieren ningun injerto. Se agrupa este subconjunto limpio,
#   se corre enriquecimiento funcional sobre los clusteres resultantes
#   y se reporta si la coherencia biologica aparece sin injerto.
#
# LOGICA DE INTERPRETACION:
#   Resultado POSITIVO (terminos GO significativos) -> confirmacion definitiva
#   Resultado NEGATIVO                              -> ambiguo (menor potencia
#                                                     con 324 genes vs 15868)
#
# FUENTES (todo ya calculado):
#   - ranking_medioide_OrthoMaM.xlsx -> IDs de los 324 arboles completos
#   - cache/bosque_*.rds             -> arboles para recalcular RF 324x324
#   - cache/hgnc_conversion.rds      -> mapeo NCBI -> HGNC (reutilizado)
#
# SALIDAS:
#   - Resultados/control_324_clustering_OrthoMaM.xlsx
#   - Resultados/control_324_enrichment_OrthoMaM.xlsx
#   - Resultados/control_324_burbujas_OrthoMaM.png
#   - Resultados/control_324_resumen_OrthoMaM.csv
# =============================================================================

library(here)
library(ape)
library(TreeDist)
library(cluster)
library(gprofiler2)
library(openxlsx)
library(ggplot2)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

cat("=============================================================\n")
cat(" CONTROL: 324 arboles completos sin injerto\n")
cat(" Observacion 1.2 - Retroalimentacion informe de tesis\n")
cat("=============================================================\n\n")


# =============================================================================
# 1. IDs DE LOS 324 ARBOLES COMPLETOS
# =============================================================================
cat("--- [1] Cargando IDs de los 324 arboles completos ---\n")

ruta_medioide <- file.path(DIR_RESULTS,
                           paste0("ranking_medioide_", NOMBRE_BDD, ".xlsx"))
if (!file.exists(ruta_medioide))
  stop("Falta: ", ruta_medioide, "\nEjecuta 02_generar_arbol_medioide.R")

ranking_df  <- read.xlsx(ruta_medioide, sheet = 1, startRow = 2)
ids_324     <- as.character(ranking_df$nombre_arbol)
cat(sprintf("  %d arboles completos identificados\n", length(ids_324)))
cat(sprintf("  Ejemplo ID: %s\n", ids_324[1]))


# =============================================================================
# 2. CARGAR BOSQUE Y EXTRAER SOLO LOS 324 ARBOLES
#    Usamos el bosque original (sin injertar) para calcular RF limpio
# =============================================================================
cat("\n--- [2] Cargando bosque original y extrayendo subconjunto ---\n")

archivos_bosque <- list.files(DIR_CACHE,
                              pattern = "^bosque_.*\\.rds$",
                              full.names = TRUE)
if (length(archivos_bosque) == 0)
  stop("No hay bosque en cache. Ejecuta el pipeline primero.")

bosque_rds <- archivos_bosque[which.min(file.size(archivos_bosque))]
cat(sprintf("  Cargando: %s (%.1f MB)\n",
            basename(bosque_rds), file.size(bosque_rds) / 1e6))

bosque_total <- readRDS(bosque_rds)
cat(sprintf("  Total arboles en bosque: %d\n", length(bosque_total)))

# Filtrar solo los 324 con todas las especies
bosque_324 <- bosque_total[names(bosque_total) %in% ids_324]
bosque_324 <- bosque_324[!sapply(bosque_324, is.null)]
class(bosque_324) <- "multiPhylo"
rm(bosque_total); gc(verbose = FALSE)

cat(sprintf("  Arboles seleccionados: %d\n", length(bosque_324)))

# Verificar que todos tienen 190 hojas
n_hojas_ctrl <- sapply(bosque_324, function(t) length(t$tip.label))
cat(sprintf("  Hojas por arbol: min=%d max=%d (deben ser todos 190)\n",
            min(n_hojas_ctrl), max(n_hojas_ctrl)))


# =============================================================================
# 3. CALCULAR MATRIZ RF 324x324
# =============================================================================
cat("\n--- [3] Calculando matriz RF 324x324 ---\n")

ruta_cache_324 <- file.path(DIR_CACHE, "matriz_rf_324.rds")

if (file.exists(ruta_cache_324)) {
  cat("  Matriz 324x324 encontrada en cache. Cargando...\n")
  mat_324 <- readRDS(ruta_cache_324)
} else {
  cat("  Calculando RF normalizado (puede tomar 1-2 min)...\n")
  t_mat <- system.time({
    dist_obj <- RobinsonFoulds(bosque_324, normalize = TRUE)
    mat_324  <- as.matrix(dist_obj)
  })
  saveRDS(mat_324, ruta_cache_324)
  cat(sprintf("  Calculado en %.1f s. Cache guardado.\n", t_mat["elapsed"]))
}

rm(bosque_324); gc(verbose = FALSE)
cat(sprintf("  Matriz: %d x %d\n", nrow(mat_324), ncol(mat_324)))
cat(sprintf("  Rango valores: [%.4f, %.4f]\n",
            min(mat_324[upper.tri(mat_324)]),
            max(mat_324[upper.tri(mat_324)])))


# =============================================================================
# 4. CLUSTERING CLARA — seleccion de k por Silhouette
# =============================================================================
cat("\n--- [4] Clustering CLARA sobre submatriz 324x324 ---\n")

K_RANGO   <- 2:10
SEED      <- 42
MUESTRAS  <- 50

set.seed(SEED)
resultados_k <- list()
silhouettes_k <- numeric(length(K_RANGO))
names(silhouettes_k) <- as.character(K_RANGO)

for (k in K_RANGO) {
  res_k <- tryCatch(
    clara(mat_324, k = k, metric = "euclidean",
          samples  = MUESTRAS,
          sampsize = min(nrow(mat_324), 40 + 2 * k),
          keep.data = FALSE, rngR = TRUE),
    error = function(e) { cat("  ERROR k=", k, ":", e$message, "\n"); NULL }
  )
  if (!is.null(res_k)) {
    resultados_k[[as.character(k)]] <- res_k
    silhouettes_k[as.character(k)]  <- res_k$silinfo$avg.width
    cat(sprintf("  k=%2d  Silhouette=%.4f\n", k, res_k$silinfo$avg.width))
  }
}

k_optimo <- as.integer(names(which.max(silhouettes_k)))
cat(sprintf("\n  K optimo: %d (Silhouette = %.4f)\n",
            k_optimo, silhouettes_k[as.character(k_optimo)]))

clustering_final <- setNames(
  as.integer(resultados_k[[as.character(k_optimo)]]$clustering),
  rownames(mat_324)
)

# Tabla de asignaciones
asig_df <- data.frame(
  Arbol   = names(clustering_final),
  Cluster = clustering_final,
  stringsAsFactors = FALSE
)
asig_df <- asig_df[order(asig_df$Cluster), ]
rownames(asig_df) <- NULL

tamanos <- as.data.frame(table(Cluster = asig_df$Cluster))
tamanos$Cluster <- as.integer(as.character(tamanos$Cluster))
colnames(tamanos)[2] <- "Tamano"
cat("  Tamanos de clusteres:\n")
print(tamanos)

# Silhouette por k — tabla
sil_df <- data.frame(
  k          = K_RANGO,
  Silhouette = as.numeric(silhouettes_k),
  stringsAsFactors = FALSE
)

# Exportar clustering
wb_cl <- createWorkbook()
wb_cl <- agregar_hoja_formateada(wb_cl, "Asignaciones",
  paste0("Asignaciones Finales — 324 arboles completos (k=", k_optimo, ")"),
  asig_df, anchos_col = "auto")
wb_cl <- agregar_hoja_formateada(wb_cl, "Tamanos",
  "Tamanos de Clusteres",
  tamanos, anchos_col = "auto")
wb_cl <- agregar_hoja_formateada(wb_cl, "Silhouette_por_k",
  "Silhouette promedio por k",
  sil_df, anchos_col = "auto")

ruta_cl <- file.path(DIR_RESULTS,
  paste0("control_324_clustering_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb_cl, ruta_cl, overwrite = TRUE)
cat("  Clustering guardado:", ruta_cl, "\n")

rm(resultados_k, mat_324); gc(verbose = FALSE)


# =============================================================================
# 5. CONVERSION NCBI -> HGNC (usando cache existente)
# =============================================================================
cat("\n--- [5] Mapeando Gene IDs a HGNC ---\n")

ruta_cache_hgnc <- file.path(DIR_CACHE, "hgnc_conversion.rds")
if (!file.exists(ruta_cache_hgnc))
  stop("No se encontro el cache HGNC: ", ruta_cache_hgnc,
       "\nEjecuta primero 08_etiquetas.R")

hgnc_df <- readRDS(ruta_cache_hgnc)
cat(sprintf("  HGNC cache cargado: %d entradas\n", nrow(hgnc_df)))

# Extraer Gene_ID del nombre del arbol (numero antes del primer _)
asig_df$Gene_ID <- sub("^(\\d+)_.*$", "\\1", asig_df$Arbol)

# Unir
genes_ctrl <- merge(asig_df, hgnc_df[, c("Gene_ID", "HGNC_Symbol")],
                    by = "Gene_ID", all.x = TRUE)
genes_ctrl <- genes_ctrl[order(genes_ctrl$Cluster), ]
rownames(genes_ctrl) <- NULL

n_hgnc <- sum(!is.na(genes_ctrl$HGNC_Symbol))
cat(sprintf("  Con simbolo HGNC: %d de %d (%.1f%%)\n",
            n_hgnc, nrow(genes_ctrl),
            100 * n_hgnc / nrow(genes_ctrl)))


# =============================================================================
# 6. ENRIQUECIMIENTO FUNCIONAL POR CLUSTER
# =============================================================================
cat("\n--- [6] Enriquecimiento funcional (gprofiler2) ---\n")

ORGANISMO      <- "hsapiens"
FUENTES        <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")
P_UMBRAL       <- 0.05
MIN_GENES      <- 3

clusters_unicos <- sort(unique(asig_df$Cluster))
cat(sprintf("  Clusters a analizar: %d\n", length(clusters_unicos)))

# Cache de enriquecimiento especifico para el control 324
dir_cache_ctrl <- file.path(DIR_CACHE, "enrichment_324")
if (!dir.exists(dir_cache_ctrl))
  dir.create(dir_cache_ctrl, recursive = TRUE)

resultados_enrich <- list()
resumen_enrich    <- data.frame()

for (cid in clusters_unicos) {
  simbolos_cid <- genes_ctrl$HGNC_Symbol[
    genes_ctrl$Cluster == cid & !is.na(genes_ctrl$HGNC_Symbol)
  ]
  simbolos_cid <- simbolos_cid[simbolos_cid != ""]

  ruta_cache_cid <- file.path(dir_cache_ctrl,
    paste0("ctrl324_enrich_cluster_", cid, ".rds"))

  if (file.exists(ruta_cache_cid)) {
    cat(sprintf("  Cluster %d: cargando desde cache.\n", cid))
    res_cid <- readRDS(ruta_cache_cid)
  } else {
    cat(sprintf("  Cluster %d: %d genes HGNC -> gprofiler2...\n",
                cid, length(simbolos_cid)))
    if (length(simbolos_cid) < MIN_GENES) {
      cat(sprintf("    Menos de %d genes. Saltando.\n", MIN_GENES))
      res_cid <- NULL
    } else {
      res_gost <- tryCatch(
        gost(query             = simbolos_cid,
             organism          = ORGANISMO,
             sources           = FUENTES,
             correction_method = "fdr",
             user_threshold    = P_UMBRAL,
             significant       = TRUE,
             evcodes           = FALSE),
        error = function(e) {
          cat(sprintf("    ERROR: %s\n", e$message)); NULL
        }
      )
      if (!is.null(res_gost) && !is.null(res_gost$result) &&
          nrow(res_gost$result) > 0) {
        res_cid <- res_gost$result
        res_cid <- res_cid[res_cid$intersection_size >= MIN_GENES, ]
        res_cid$Cluster <- cid
        res_cid <- res_cid[order(res_cid$p_value), ]
        cat(sprintf("    %d terminos significativos\n", nrow(res_cid)))
      } else {
        cat("    Sin terminos significativos.\n")
        res_cid <- NULL
      }
    }
    saveRDS(res_cid, ruta_cache_cid)
    Sys.sleep(0.4)
  }

  resultados_enrich[[as.character(cid)]] <- res_cid

  resumen_enrich <- rbind(resumen_enrich, data.frame(
    Cluster        = cid,
    N_Genes_Input  = nrow(genes_ctrl[genes_ctrl$Cluster == cid, ]),
    N_Genes_HGNC   = length(simbolos_cid),
    N_Terminos_Sig = if (!is.null(res_cid)) nrow(res_cid) else 0L,
    N_GO_BP        = if (!is.null(res_cid)) sum(res_cid$source == "GO:BP") else 0L,
    N_KEGG         = if (!is.null(res_cid)) sum(res_cid$source == "KEGG")  else 0L,
    Top_Termino    = if (!is.null(res_cid) && nrow(res_cid) > 0)
                       res_cid$term_name[1] else "—",
    stringsAsFactors = FALSE
  ))
}

n_clusters_sig <- sum(resumen_enrich$N_Terminos_Sig > 0)
n_terminos_tot <- sum(resumen_enrich$N_Terminos_Sig)
cat(sprintf("\n  Clusters con terminos GO significativos: %d de %d\n",
            n_clusters_sig, length(clusters_unicos)))
cat(sprintf("  Total terminos significativos: %d\n", n_terminos_tot))


# =============================================================================
# 7. GRAFICO DE BURBUJAS
# =============================================================================
cat("\n--- [7] Generando grafico de burbujas ---\n")

resultados_df <- do.call(rbind,
  Filter(Negate(is.null), resultados_enrich))
rownames(resultados_df) <- NULL

ruta_png <- file.path(DIR_RESULTS,
  paste0("control_324_burbujas_", NOMBRE_BDD, ".png"))

if (!is.null(resultados_df) && nrow(resultados_df) > 0) {

  # Filtrar a GO:BP y top 20 terminos por frecuencia en clusteres
  gobp <- resultados_df[resultados_df$source == "GO:BP", ]

  if (nrow(gobp) > 0) {
    freq_terminos <- sort(table(gobp$term_name), decreasing = TRUE)
    top_terminos  <- names(freq_terminos)[seq_len(min(20, length(freq_terminos)))]
    gobp_plot     <- gobp[gobp$term_name %in% top_terminos, ]

    # Truncar nombres
    gobp_plot$term_short <- ifelse(
      nchar(gobp_plot$term_name) > 45,
      paste0(substr(gobp_plot$term_name, 1, 42), "..."),
      gobp_plot$term_name
    )
    gobp_plot$log10_fdr  <- -log10(gobp_plot$p_value)
    gobp_plot$Cluster    <- factor(gobp_plot$Cluster)

    p_burbuja <- ggplot(gobp_plot,
                        aes(x = Cluster,
                            y = reorder(term_short, log10_fdr),
                            size = intersection_size,
                            color = log10_fdr)) +
      geom_point(alpha = 0.85) +
      scale_color_gradient(low = "#74C4DE", high = "#C0392B",
                           name = expression(-log[10](FDR))) +
      scale_size_continuous(name = "Genes en\ntermino", range = c(2, 9)) +
      labs(
        title    = "Control: enriquecimiento GO:BP en los 324 arboles sin injerto",
        subtitle = sprintf(
          "k = %d clusteres | %d terminos GO:BP significativos | n = 324 arboles (190 sp. completas)",
          k_optimo, nrow(gobp_plot)
        ),
        x        = "Cluster",
        y        = NULL,
        caption  = paste0(
          "FDR < 0.05 (gprofiler2) | Obs. 1.2: grupo de control libre de injerto\n",
          "Resultado positivo confirma coherencia funcional independiente del injerto."
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 9),
        plot.caption  = element_text(color = "gray60", size = 8),
        axis.text.y   = element_text(size = 8),
        legend.position = "right"
      )

    ggsave(ruta_png, plot = p_burbuja,
           width = 10,
           height = max(6, length(unique(gobp_plot$term_short)) * 0.38 + 3),
           dpi = 300)
    cat("  Grafico guardado:", ruta_png, "\n")
  } else {
    cat("  Sin terminos GO:BP para graficar.\n")
  }
} else {
  cat("  Sin terminos significativos en ningun cluster.\n")
}


# =============================================================================
# 8. EXPORTAR EXCEL DE ENRIQUECIMIENTO
# =============================================================================
cat("\n--- [8] Exportando Excel de enriquecimiento ---\n")

wb_en <- createWorkbook()
wb_en <- agregar_hoja_formateada(wb_en, "Resumen",
  paste0("Resumen enriquecimiento — 324 arboles completos (k=", k_optimo, ")"),
  resumen_enrich, anchos_col = "auto")

if (!is.null(resultados_df) && nrow(resultados_df) > 0) {
  cols_keep <- intersect(
    c("Cluster", "source", "term_id", "term_name", "p_value",
      "intersection_size", "term_size", "query_size"),
    colnames(resultados_df)
  )
  wb_en <- agregar_hoja_formateada(wb_en, "Todos_Terminos",
    "Todos los terminos significativos",
    resultados_df[, cols_keep], anchos_col = "auto")

  for (cid in names(resultados_enrich)) {
    res_cid <- resultados_enrich[[cid]]
    if (is.null(res_cid) || nrow(res_cid) == 0) next
    nombre_h <- paste0("Cluster_", cid)
    wb_en <- agregar_hoja_formateada(wb_en, nombre_h,
      paste0("Cluster ", cid, " — ", nrow(res_cid), " terminos"),
      res_cid[, cols_keep], anchos_col = "auto")
  }
}

ruta_en <- file.path(DIR_RESULTS,
  paste0("control_324_enrichment_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb_en, ruta_en, overwrite = TRUE)
cat("  Excel guardado:", ruta_en, "\n")


# =============================================================================
# 9. CSV RESUMEN + TEXTO LaTeX
# =============================================================================
cat("\n--- [9] Resumen y texto LaTeX ---\n")

resumen_csv <- data.frame(
  Metrica = c("n_arboles_control", "n_arboles_con_hgnc",
              "k_optimo", "silhouette_k_optimo",
              "n_clusters_con_go_sig", "n_clusters_total",
              "n_terminos_gobp_sig", "n_terminos_total_sig"),
  Valor   = c(nrow(asig_df),
              n_hgnc,
              k_optimo,
              round(silhouettes_k[as.character(k_optimo)], 4),
              n_clusters_sig,
              length(clusters_unicos),
              if (!is.null(resultados_df) && nrow(resultados_df) > 0)
                sum(resultados_df$source == "GO:BP") else 0L,
              n_terminos_tot),
  stringsAsFactors = FALSE
)

ruta_csv <- file.path(DIR_RESULTS,
  paste0("control_324_resumen_", NOMBRE_BDD, ".csv"))
write.csv(resumen_csv, ruta_csv, row.names = FALSE)
cat("  CSV guardado:", ruta_csv, "\n")

cat("\n=============================================================\n")
cat(" RESULTADO FINAL\n")
cat("=============================================================\n")

cat(sprintf("  Arboles analizados   : %d (100%% completos, sin injerto)\n",
            nrow(asig_df)))
cat(sprintf("  k optimo             : %d (Silhouette = %.4f)\n",
            k_optimo, silhouettes_k[as.character(k_optimo)]))
cat(sprintf("  Clusters con GO sig  : %d / %d\n",
            n_clusters_sig, length(clusters_unicos)))
cat(sprintf("  Terminos sig totales : %d\n", n_terminos_tot))

cat("\n=============================================================\n")
cat(" TEXTO PARA SECCION 3.5 (LaTeX)\n")
cat("=============================================================\n\n")

if (n_clusters_sig > 0) {
  interpretacion <- paste0(
    "Este resultado es \\textbf{positivo}: ",
    n_clusters_sig, " de ", length(clusters_unicos),
    " cl\\'{u}steres presentaron al menos un t\\'{e}rmino GO ",
    "significativo (FDR $< 0{,}05$), acumulando ", n_terminos_tot,
    " t\\'{e}rminos en total. Dado que este subconjunto no fue ",
    "sometido a ning\\'{u}n procedimiento de injerto filogen\\'{e}tico, ",
    "la aparici\\'{o}n de coherencia funcional confirma que la se\\~{n}al ",
    "biol\\'{o}gica detectada en el an\\'{a}lisis principal es intr\\'{i}nseca ",
    "a la variabilidad topol\\'{o}gica de los genes ort\\'{o}logos y no un ",
    "artefacto de la homogeneizaci\\'{o}n por injerto."
  )
} else {
  interpretacion <- paste0(
    "Este resultado es \\textbf{ambiguo}: ninguno de los ",
    length(clusters_unicos),
    " cl\\'{u}steres alcanz\\'{o} significancia estad\\'{i}stica en ",
    "el enriquecimiento funcional (FDR $< 0{,}05$). Sin embargo, con ",
    "\\'{u}nicamente ", nrow(asig_df), " genes (frente a los ",
    "$15.868$ del an\\'{a}lisis principal), la potencia estad\\'{i}stica ",
    "del enriquecimiento es inherentemente menor. Un resultado negativo ",
    "en este contexto no refuta la coherencia funcional observada en el ",
    "corpus completo; simplemente indica que el subconjunto de control ",
    "no tiene tama\\~{n}o suficiente para detectarla de forma robusta."
  )
}

cat(paste0(
  "\\subsection{Grupo de control: \\'{a}rboles sin injerto (324 \\'{a}rboles completos)}\n\n",
  "Como experimento de control definitivo para la Observaci\\'{o}n~1.2, se ",
  "agreg\\'{o} el subconjunto de los \\textbf{324 \\'{a}rboles filogen\\'{e}ticos} ",
  "que ya conten\\'{i}an las \\textbf{190 especies} de forma completa (sin ",
  "requerir ning\\'{u}n injerto). Sobre estos \\'{a}rboles se calcul\\'{o} una ",
  "nueva matriz de distancias de Robinson-Foulds $324 \\\\times 324$ exclusivamente ",
  "con topolog\\'{i}as originales, se aplic\\'{o} CLARA con selecci\\'{o}n autom\\'{a}tica ",
  "de $k$ por el coeficiente Silhouette ($k_{\\\\text{optimo}} = ", k_optimo, "$, ",
  "Silhouette $= ", round(silhouettes_k[as.character(k_optimo)], 4), "$) y se ",
  "ejecut\\'{o} el an\\'{a}lisis de enriquecimiento funcional mediante \\texttt{gprofiler2} ",
  "(FDR $< 0{,}05$, GO:BP, KEGG y Reactome).\n\n",
  interpretacion, "\n\n",
  "\\begin{figure}[H]\n",
  "    \\centering\n",
  "    \\includegraphics[width=0.85\\textwidth]{images/control_324_burbujas_", NOMBRE_BDD, ".png}\n",
  "    \\caption{Enriquecimiento funcional (GO:BP) sobre los ",
  nrow(asig_df), " \\'{a}rboles filogen\\'{e}ticos completos sin injerto ",
  "($k = ", k_optimo, "$ cl\\'{u}steres). El tama\\~{n}o de cada burbuja es ",
  "proporcional al n\\'{u}mero de genes del cl\\'{u}ster presentes en el t\\'{e}rmino ",
  "y el color representa la significancia ($-\\\\log_{10}$(FDR)).}\n",
  "    \\label{fig:control_324_burbujas}\n",
  "\\end{figure}\n"
))

cat("\n=== COMPLETADO ===\n")