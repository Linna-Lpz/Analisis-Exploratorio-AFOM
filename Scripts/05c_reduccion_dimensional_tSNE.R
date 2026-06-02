# =============================================================================
# 05c_reduccion_dimensional_tSNE.R
# t-SNE sobre la matriz RF normalizada — con caché .rds
# Salidas: coordenadas en Excel, gráfico ggplot2, caché .rds
# =============================================================================

library(Rtsne)
library(ggplot2)
library(patchwork)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# CARGAR MATRIZ RF DESDE CACHÉ (generada en script 04)
# =============================================================================
cat("=== CARGANDO MATRIZ RF ===\n")

ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché: ", ruta_cache_matriz,
       "\nEjecuta primero el script 04_calcular_matriz.R")
}

matriz_cuadrada <- readRDS(ruta_cache_matriz)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# Verificar duplicados exactos — Rtsne falla si hay filas idénticas
n_duplicados <- sum(duplicated(matriz_cuadrada))
if (n_duplicados > 0) {
  cat(sprintf(
    "[aviso] Se detectaron %d filas duplicadas en la matriz.\n", n_duplicados
  ))
  cat("  Rtsne puede fallar con duplicados. Considera check_duplicates = FALSE.\n")
}

# =============================================================================
# CALCULAR t-SNE — con caché condicional
# =============================================================================
cat("\n=== CALCULANDO t-SNE ===\n")

ruta_cache_tsne <- file.path(DIR_CACHE, "tsne_coords.rds")

# --- Hiperparámetros ---
# Perplexity: controla el balance entre estructura local y global.
# Regla práctica: entre 5 y 50. Con n grande se puede subir hasta 100.
# Debe cumplir: perplexity < (n - 1) / 3
PERPLEXITY   <- min(50, floor((nrow(matriz_cuadrada) - 1) / 3) - 1)
N_COMPONENTS <- 2       # dimensiones de salida
MAX_ITER     <- 1000    # iteraciones de optimización
THETA        <- 0.0     # aproximación Barnes-Hut: 0 = exacto (lento), 1 = rápido (menos preciso)
SEED         <- 42
N_THREADS    <- parallel::detectCores() - 1  # usar todos los núcleos menos uno

cat(sprintf(
  "Hiperparámetros: perplexity = %d | max_iter = %d | theta = %.1f | threads = %d\n",
  PERPLEXITY, MAX_ITER, THETA, N_THREADS
))

if (file.exists(ruta_cache_tsne)) {
  cat("t-SNE encontrado en caché. Cargando...\n")
  tiempo_tsne <- system.time({
    tsne_resultado <- readRDS(ruta_cache_tsne)
  })
  cat("Cargado desde caché.\n")
  
} else {
  cat("Calculando t-SNE (puede tardar varios minutos con",
      nrow(matriz_cuadrada), "árboles)...\n\n")
  
  tiempo_tsne <- system.time({
    set.seed(SEED)
    tsne_resultado <- Rtsne(
      X                = matriz_cuadrada,
      dims             = N_COMPONENTS,
      perplexity       = PERPLEXITY,
      max_iter         = MAX_ITER,
      theta            = THETA,
      is_distance      = TRUE,      # matriz de distancias precomputada
      pca              = FALSE,     # omitir PCA: ya tenemos distancias
      normalize        = FALSE,     # no renormalizar: RF ya está normalizada [0,1]
      check_duplicates = n_duplicados == 0,  # desactivar si hay duplicados
      verbose          = TRUE,
      num_threads      = N_THREADS
    )
    saveRDS(tsne_resultado, file = ruta_cache_tsne)
  })
  cat("\nt-SNE calculado y guardado en caché.\n")
}

cat(sprintf("Tiempo t-SNE : %.1f segundos\n", tiempo_tsne["elapsed"]))
cat(sprintf("Costo final  : %.6f\n", tail(tsne_resultado$itercosts, 1)))

# =============================================================================
# CONSTRUIR DATAFRAME DE COORDENADAS
# =============================================================================
coords_df <- data.frame(
  Arbol  = rownames(matriz_cuadrada),
  tSNE_1 = tsne_resultado$Y[, 1],
  tSNE_2 = tsne_resultado$Y[, 2],
  stringsAsFactors = FALSE
)

cat("\nPrimeras filas de coordenadas:\n")
print(head(coords_df, 5))

# =============================================================================
# LEER Y UNIR ETIQUETAS DE CLUSTERING (k óptimo + k=10, k=15)
# =============================================================================
resultado    <- unir_etiquetas_clustering(coords_df, DIR_RESULTS)
coords_df    <- resultado$coords_df
k_optimos    <- resultado$k_optimos
k_extra      <- resultado$k_extra

k_km  <- ifelse(is.na(k_optimos["KMeans"]), "?", k_optimos["KMeans"])
k_pam <- ifelse(is.na(k_optimos["PAM"]),    "?", k_optimos["PAM"])
k_cl  <- ifelse(is.na(k_optimos["CLARA"]),  "?", k_optimos["CLARA"])

# =============================================================================
# FUNCIÓN PARA GRAFICAR CLUSTERING EN tSNE
# =============================================================================
graficar_clustering_tsne <- function(df, col_cluster, titulo, subtitulo = "") {
  if (!col_cluster %in% colnames(df)) {
    warning(sprintf("Columna '%s' no encontrada. Saltando gráfico.", col_cluster))
    return(NULL)
  }
  n_clusters <- length(unique(na.omit(df[[col_cluster]])))
  paleta <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
              "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
              "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A")
  
  ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = .data[[col_cluster]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = paleta[seq_len(n_clusters)],
                       name = "Cluster", na.value = "grey80") +
    labs(title = titulo, subtitle = subtitulo,
         x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(color = "gray40", size = 9),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# =============================================================================
# GRÁFICO PRINCIPAL — distribución de árboles (sin clustering)
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS ===\n")

costo_final <- round(tail(tsne_resultado$itercosts, 1), 4)

p_tsne <- ggplot(coords_df, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(
    color = "forestgreen",
    alpha = 0.6,
    size  = 1.5
  ) +
  labs(
    title    = paste0("t-SNE — Matriz RF Normalizada (", NOMBRE_BDD, ")"),
    subtitle = paste0(
      "perplexity = ", PERPLEXITY,
      "  |  max_iter = ", MAX_ITER,
      "  |  theta = ", THETA,
      "  |  n = ", nrow(coords_df), " árboles"
    ),
    x       = "t-SNE 1",
    y       = "t-SNE 2",
    caption = paste0(
      "Distancia Robinson-Foulds normalizada | Rtsne (Barnes-Hut)",
      "  |  Costo final KL = ", costo_final
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8)
  )

ruta_grafico <- file.path(DIR_RESULTS, paste0("tsne_", NOMBRE_BDD, ".png"))
ggsave(ruta_grafico, plot = p_tsne, width = 10, height = 7, dpi = 300)
cat("Gráfico t-SNE guardado:", ruta_grafico, "\n")

# =============================================================================
# GRÁFICOS CLUSTERING — K ÓPTIMO
# =============================================================================
cat("\n=== GRÁFICOS CLUSTERING tSNE (K ÓPTIMO) ===\n")

p_km <- graficar_clustering_tsne(coords_df, "Cluster_KMeans", "K-Means",
                                  sprintf("k = %s  |  n = %d árboles", k_km, nrow(coords_df)))
p_pa <- graficar_clustering_tsne(coords_df, "Cluster_PAM", "PAM (K-Medoids)",
                                  sprintf("k = %s  |  n = %d árboles", k_pam, nrow(coords_df)))
p_cl <- graficar_clustering_tsne(coords_df, "Cluster_CLARA", "CLARA",
                                  sprintf("k = %s  |  n = %d árboles", k_cl, nrow(coords_df)))

plots_disp <- Filter(Negate(is.null), list(p_km, p_pa, p_cl))
if (length(plots_disp) > 0) {
  panel <- wrap_plots(plots_disp, ncol = min(length(plots_disp), 3)) +
    plot_annotation(
      title   = paste0("Comparación de Clustering — t-SNE (", NOMBRE_BDD, ")"),
      caption = paste0("Distancia RF normalizada  |  Rtsne  |  perplexity=", PERPLEXITY),
      theme   = theme(
        plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.caption = element_text(color = "gray50", size = 8, hjust = 0.5)
      )
    )
  ruta_panel <- file.path(DIR_RESULTS, paste0("tsne_comparativo_clustering_", NOMBRE_BDD, ".png"))
  ggsave(ruta_panel, plot = panel, width = 7 * length(plots_disp), height = 6, dpi = 300)
  cat("Panel comparativo guardado:", ruta_panel, "\n")
}

# =============================================================================
# GRÁFICOS CLUSTERING — K EXTRA (k=10, k=15, etc.)
# =============================================================================
for (ke in k_extra) {
  cat(sprintf("\n=== GRÁFICOS tSNE (K = %d) ===\n", ke))
  
  col_km  <- paste0("Cluster_KMeans_K", ke)
  col_pam <- paste0("Cluster_PAM_K",    ke)
  col_cl  <- paste0("Cluster_CLARA_K",  ke)
  
  pe_km <- graficar_clustering_tsne(coords_df, col_km, "K-Means",
                                     sprintf("k = %d  |  n = %d árboles", ke, nrow(coords_df)))
  pe_pa <- graficar_clustering_tsne(coords_df, col_pam, "PAM (K-Medoids)",
                                     sprintf("k = %d  |  n = %d árboles", ke, nrow(coords_df)))
  pe_cl <- graficar_clustering_tsne(coords_df, col_cl, "CLARA",
                                     sprintf("k = %d  |  n = %d árboles", ke, nrow(coords_df)))
  
  plots_extra <- Filter(Negate(is.null), list(pe_km, pe_pa, pe_cl))
  if (length(plots_extra) > 0) {
    panel_extra <- wrap_plots(plots_extra, ncol = min(length(plots_extra), 3)) +
      plot_annotation(
        title   = paste0("Comparación de Clustering (k=", ke, ") — t-SNE (", NOMBRE_BDD, ")"),
        caption = paste0("Distancia RF normalizada  |  Rtsne  |  perplexity=", PERPLEXITY),
        theme   = theme(
          plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.caption = element_text(color = "gray50", size = 8, hjust = 0.5)
        )
      )
    ruta_extra <- file.path(DIR_RESULTS, paste0("tsne_comparativo_K", ke, "_", NOMBRE_BDD, ".png"))
    ggsave(ruta_extra, plot = panel_extra, width = 7 * length(plots_extra), height = 6, dpi = 300)
    cat("Panel comparativo tSNE k=", ke, " guardado:", ruta_extra, "\n")
  }
}

# =============================================================================
# GRÁFICO DE CONVERGENCIA — evolución del costo KL por iteración
# Permite verificar si el algoritmo convergió antes de max_iter
# =============================================================================
convergencia_df <- data.frame(
  iteracion = seq(50, MAX_ITER, by = 50)[seq_along(tsne_resultado$itercosts)],
  costo_KL  = tsne_resultado$itercosts
)

p_conv <- ggplot(convergencia_df, aes(x = iteracion, y = costo_KL)) +
  geom_line(color = "forestgreen", linewidth = 0.8) +
  geom_point(size = 1, color = "forestgreen") +
  labs(
    title    = paste0("Convergencia t-SNE (", NOMBRE_BDD, ")"),
    subtitle = paste0(
      "Costo final = ", costo_final,
      "  |  Si la curva se aplana antes de ", MAX_ITER,
      " iter → convergió correctamente"
    ),
    x       = "Iteración",
    y       = "Divergencia KL",
    caption = "Un descenso suave que se estabiliza indica buena convergencia"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8)
  )

ruta_conv <- file.path(DIR_RESULTS, paste0("tsne_convergencia_", NOMBRE_BDD, ".png"))
ggsave(ruta_conv, plot = p_conv, width = 8, height = 5, dpi = 300)
cat("Gráfico de convergencia guardado:", ruta_conv, "\n")

# =============================================================================
# EXPORTAR COORDENADAS Y DIAGNÓSTICO A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

parametros_df <- data.frame(
  Parametro = c("metodo", "libreria", "perplexity", "max_iter", "theta",
                "is_distance", "pca", "normalize", "seed", "n_threads",
                "n_arboles", "costo_KL_final", "tiempo_calculo_s"),
  Valor     = c("t-SNE (Barnes-Hut)", "Rtsne", PERPLEXITY, MAX_ITER, THETA,
                "TRUE", "FALSE", "FALSE", SEED, N_THREADS,
                nrow(coords_df), costo_final, round(tiempo_tsne["elapsed"], 1)),
  stringsAsFactors = FALSE
)

convergencia_export_df <- convergencia_df
colnames(convergencia_export_df) <- c("Iteracion", "Costo_KL")

wb <- createWorkbook()

agregar_hoja_formateada(wb, "Coordenadas_tSNE",
                        paste0("Coordenadas t-SNE + Etiquetas de Clustering — ", NOMBRE_BDD),
                        coords_df,
                        anchos_col = "auto")

agregar_hoja_formateada(wb, "Diagnostico",
                        "Diagnóstico y Parámetros t-SNE",
                        parametros_df,
                        anchos_col = c(25, 30))

agregar_hoja_formateada(wb, "Convergencia_KL",
                        "Evolución del Costo KL por Iteración",
                        convergencia_export_df,
                        anchos_col = c(15, 18))

ruta_excel <- file.path(DIR_RESULTS, paste0("tsne_coordenadas_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Coordenadas guardadas:", ruta_excel, "\n")

cat("\n=== COMPLETADO ===\n")