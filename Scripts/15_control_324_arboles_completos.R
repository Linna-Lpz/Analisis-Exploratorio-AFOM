# =============================================================================
# 15_control_324_arboles_completos.R
# =============================================================================

library(here)
library(ape)
library(TreeDist)
library(cluster)
library(gprofiler2)
library(openxlsx)
library(ggplot2)
library(mstknnclust)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# Metodo de clustering a utilizar: "clara", "kmeans" o "mst-knn"
METODO_CLUSTERING <- "mst-knn"

# =============================================================================
# 1. CARGAR CONJUNTO CORE SET (ARBOLES COMPLETOS) DESDE CACHE
#    Generado en 02_filtro_arboles_medioide.R
# =============================================================================
cat("--- [1] Cargando conjunto core set (arboles completos) desde cache ---\n")

ruta_cache_core_set <- file.path(DIR_CACHE, "conjunto_core_set.rds")
if (!file.exists(ruta_cache_core_set))
  stop("Falta: ", ruta_cache_core_set,
       "\nEjecuta primero 02_filtro_arboles_medioide.R")

bosque_324 <- readRDS(ruta_cache_core_set)
class(bosque_324) <- "multiPhylo"

cat(sprintf("  Arboles cargados: %d\n", length(bosque_324)))
cat(sprintf("  Ejemplo ID: %s\n", names(bosque_324)[1]))

# Guardamos los nombres reales para poder reasignarlos a la matriz RF
# mas adelante, por si RobinsonFoulds() no los conserva como dimnames
nombres_ids <- names(bosque_324)

# Verificar que todos los arboles tienen el mismo numero de hojas (core set)
n_hojas_ctrl <- sapply(bosque_324, function(t) length(t$tip.label))
cat(sprintf("  Hojas por arbol: min=%d max=%d (deben ser todas iguales)\n",
            min(n_hojas_ctrl), max(n_hojas_ctrl)))

if (min(n_hojas_ctrl) != max(n_hojas_ctrl)) {
  warning("Los arboles del conjunto core set no tienen todos el mismo numero ",
          "de hojas. Revisa el cache: ", ruta_cache_core_set)
}


# =============================================================================
# 2. CALCULAR MATRIZ RF 324x324
# =============================================================================
cat("\n--- [2] Calculando matriz RF 324x324 ---\n")

ruta_cache_324 <- file.path(DIR_CACHE, "matriz_rf_324.rds")

if (file.exists(ruta_cache_324)) {
  cat("  Matriz 324x324 encontrada en cache. Cargando...\n")
  mat_324 <- readRDS(ruta_cache_324)
} else {
  cat("  Calculando Distancia (puede tomar 1-2 min)...\n")
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

# --- Asegurar que la matriz tenga los nombres reales de los arboles ---
# RobinsonFoulds() devuelve un objeto "dist"; si este no trae la etiqueta
# "Labels", as.matrix() NO deja rownames en NULL: les asigna nombres
# genericos "1","2",...,"n". Por eso no basta con chequear is.null();
# forzamos siempre los nombres reales capturados en nombres_ids, que estan
# en el mismo orden en que se paso bosque_324 a RobinsonFoulds().
if (!identical(rownames(mat_324), nombres_ids)) {
  cat("  AVISO: mat_324 no tenia los nombres reales de los arboles (traia nombres genericos). Reasignando desde nombres_ids.\n")
  rownames(mat_324) <- nombres_ids
  colnames(mat_324) <- nombres_ids
}


# =============================================================================
# 3. CLUSTERING (clara, kmeans o mst-knn)
# =============================================================================
cat(sprintf("\n--- [3] Clustering %s sobre submatriz 324x324 ---\n", toupper(METODO_CLUSTERING)))

K_RANGO   <- 2:10
SEED      <- 42
MUESTRAS  <- 50

set.seed(SEED)
resultados_k <- list()
silhouettes_k <- numeric()

if (METODO_CLUSTERING %in% c("clara", "kmeans")) {
  # Calcular embedding MDS para tener un espacio Euclidiano válido
  N_DIMS <- min(20, nrow(mat_324) - 1)
  cat(sprintf("  Calculando embedding MDS (k=%d) para %s...\n", N_DIMS, toupper(METODO_CLUSTERING)))
  mds_coords <- cmdscale(as.dist(mat_324), k = N_DIMS)
  
  silhouettes_k <- numeric(length(K_RANGO))
  names(silhouettes_k) <- as.character(K_RANGO)
  
  for (k in K_RANGO) {
    res_k <- tryCatch({
      if (METODO_CLUSTERING == "clara") {
        clara(x = mds_coords, k = k, metric = "euclidean",
              samples  = MUESTRAS,
              sampsize = min(nrow(mds_coords), 40 + 2 * k),
              keep.data = FALSE, rngR = TRUE)
      } else {
        kmeans(x = mds_coords, centers = k, nstart = 25)
      }
    }, error = function(e) { cat("  ERROR k=", k, ":", e$message, "\n"); NULL })
    
    if (!is.null(res_k)) {
      resultados_k[[as.character(k)]] <- res_k
      
      # Calcular silhouette sobre la matriz de distancias ORIGINAL para ambos métodos
      if (METODO_CLUSTERING == "kmeans") {
        asignaciones_temp <- res_k$cluster
      } else {
        asignaciones_temp <- res_k$clustering
      }
      
      sil <- cluster::silhouette(asignaciones_temp, as.dist(mat_324))
      sil_avg <- mean(sil[, 3])
      
      silhouettes_k[as.character(k)] <- sil_avg
      cat(sprintf("  k=%2d  Silhouette=%.4f\n", k, sil_avg))
    }
  }
  
  k_optimo <- as.integer(names(which.max(silhouettes_k)))
  cat(sprintf("\n  K optimo: %d (Silhouette = %.4f)\n",
              k_optimo, silhouettes_k[as.character(k_optimo)]))
  
  if (METODO_CLUSTERING == "clara") {
    clustering_final <- setNames(
      as.integer(resultados_k[[as.character(k_optimo)]]$clustering),
      rownames(mat_324)
    )
  } else {
    clustering_final <- setNames(
      as.integer(resultados_k[[as.character(k_optimo)]]$cluster),
      rownames(mat_324)
    )
  }
  
  sil_df <- data.frame(
    k          = K_RANGO,
    Silhouette = as.numeric(silhouettes_k),
    stringsAsFactors = FALSE
  )
  
} else if (METODO_CLUSTERING == "mst-knn") {
  mst_res <- mstknnclust::mst.knn(mat_324)
  k_optimo <- mst_res$cnumber
  
  clustering_final <- setNames(
    as.integer(mst_res$cluster),
    names(mst_res$cluster)
  )
  
  # Calcular silhouette para el mst-knn optimo
  sil <- cluster::silhouette(clustering_final, as.dist(mat_324[names(clustering_final), names(clustering_final)]))
  silhouettes_k <- setNames(mean(sil[, 3]), as.character(k_optimo))
  
  cat(sprintf("\n  K optimo: %d (Silhouette = %.4f)\n",
              k_optimo, silhouettes_k[as.character(k_optimo)]))
  
  sil_df <- data.frame(
    k          = k_optimo,
    Silhouette = as.numeric(silhouettes_k),
    stringsAsFactors = FALSE
  )
}

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
  k          = sil_df$k,
  Silhouette = as.numeric(sil_df$Silhouette),
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
# 4. CONVERSION NCBI -> HGNC (usando cache existente)
# =============================================================================
cat("\n--- [4] Mapeando Gene IDs a HGNC ---\n")

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
# 5. ENRIQUECIMIENTO FUNCIONAL POR CLUSTER
# =============================================================================
cat("\n--- [5] Enriquecimiento funcional (gprofiler2) ---\n")

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
      universo_hgnc <- unique(na.omit(hgnc_df$HGNC_Symbol))
      universo_hgnc <- universo_hgnc[universo_hgnc != ""]
      
      res_gost <- tryCatch(
        gost(query             = simbolos_cid,
             organism          = ORGANISMO,
             sources           = FUENTES,
             correction_method = "fdr",
             user_threshold    = P_UMBRAL,
             significant       = TRUE,
             custom_bg         = universo_hgnc,
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
# 6. GRAFICO DE BURBUJAS
# =============================================================================
cat("\n--- [6] Generando grafico de burbujas ---\n")

resultados_df <- do.call(rbind,
                         Filter(Negate(is.null), resultados_enrich))
rownames(resultados_df) <- NULL

ruta_png <- file.path(DIR_RESULTS,
                      paste0("control_324_burbujas_", NOMBRE_BDD, ".png"))

if (!is.null(resultados_df) && nrow(resultados_df) > 0) {
  
  # Filtrar a GO:BP
  gobp <- resultados_df[resultados_df$source == "GO:BP", ]
  
  if (nrow(gobp) > 0) {
    # Top 20 términos por clúster según p_value
    gobp_plot <- do.call(rbind, lapply(split(gobp, gobp$Cluster), function(df_cluster) {
      head(df_cluster[order(df_cluster$p_value), ], 20)
    }))
    rownames(gobp_plot) <- NULL
    
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
          "k = %d clusteres | %d terminos GO:BP en total (%d burbujas ilustradas) | n = 324 arboles (190 sp. completas)",
          k_optimo, nrow(gobp), nrow(gobp_plot)
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
# 7. EXPORTAR EXCEL DE ENRIQUECIMIENTO
# =============================================================================
cat("\n--- [7] Exportando Excel de enriquecimiento ---\n")

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
# 8. CSV RESUMEN
# =============================================================================
cat("\n--- [8] Resumen ---\n")

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

cat("=============================================================")
cat(" RESULTADO FINAL")
cat("=============================================================")

cat(sprintf("  Arboles analizados   : %d (100%% completos, sin injerto)\n",
            nrow(asig_df)))
cat(sprintf("  k optimo             : %d (Silhouette = %.4f)\n",
            k_optimo, silhouettes_k[as.character(k_optimo)]))
cat(sprintf("  Clusters con GO sig  : %d / %d\n",
            n_clusters_sig, length(clusters_unicos)))
cat(sprintf("  Terminos sig totales : %d\n", n_terminos_tot))

cat("\n=== COMPLETADO ===\n")