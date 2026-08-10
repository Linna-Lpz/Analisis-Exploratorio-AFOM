# =============================================================================
# 05_reduccion_dimensional.R
# VISUALIZACIÓN UMAP + ETIQUETAS DE CLUSTERING (K-Means, PAM, CLARA)
# =============================================================================
library(uwot)
library(ggplot2)
library(patchwork)   # para paneles comparativos
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# CARGAR MATRIZ RF
# =============================================================================
cat("=== CARGANDO MATRIZ RF ===\n")

ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

# Si no existe caché .rds, leer desde CSV
if (file.exists(ruta_cache_matriz)) {
  matriz_cuadrada <- readRDS(ruta_cache_matriz)
  cat("Matriz cargada desde caché RDS.\n")
} else {
  cat("Caché RDS no encontrado. Leyendo desde CSV...\n")
  matriz_cuadrada <- as.matrix(
    read.table(file.path(DIR_RESULTS, "matriz_rf_conjunto.csv"),
               sep = ";", header = TRUE, row.names = 1, check.names = FALSE)
  )
  saveRDS(matriz_cuadrada, ruta_cache_matriz)
  cat("Matriz guardada en caché RDS para uso futuro.\n")
}

cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# CALCULAR UMAP (con caché condicional)
# =============================================================================
cat("\n=== CALCULANDO UMAP ===\n")

ruta_cache_umap <- file.path(DIR_CACHE, "umap_coords.rds")

N_NEIGHBORS  <- 15
MIN_DIST     <- 0.1
N_COMPONENTS <- 2
SEED         <- 42

if (file.exists(ruta_cache_umap)) {
  cat("UMAP encontrado en caché. Cargando...\n")
  umap_coords <- readRDS(ruta_cache_umap)
  tiempo_umap <- NA
} else {
  cat("Calculando UMAP (puede tardar varios minutos)...\n")
  tiempo_inicio <- proc.time()
  
  set.seed(SEED)
  umap_coords <- umap(
    X           = as.dist(matriz_cuadrada),
    n_neighbors = N_NEIGHBORS,
    min_dist    = MIN_DIST,
    n_components = N_COMPONENTS,
    n_threads   = parallel::detectCores() - 1,
    verbose     = TRUE
  )
  
  tiempo_umap <- as.numeric(proc.time() - tiempo_inicio)[3]
  saveRDS(umap_coords, ruta_cache_umap)
  cat(sprintf("UMAP calculado en %.1f segundos y guardado en caché.\n", tiempo_umap))
}

cat("Dimensiones UMAP:", nrow(umap_coords), "x", ncol(umap_coords), "\n")

# =============================================================================
# CONSTRUIR DATAFRAME BASE DE COORDENADAS
# =============================================================================
coords_df <- data.frame(
  Arbol  = rownames(matriz_cuadrada),
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  stringsAsFactors = FALSE
)

# =============================================================================
# LEER Y UNIR ETIQUETAS DE CLUSTERING
# =============================================================================
resultado    <- unir_etiquetas_clustering(coords_df, DIR_RESULTS)
coords_df    <- resultado$coords_df
k_optimos    <- resultado$k_optimos
k_extra      <- resultado$k_extra

# Extraer k por método para subtítulos
k_km  <- ifelse(is.na(k_optimos["KMeans"]), "?", k_optimos["KMeans"])
k_pam <- ifelse(is.na(k_optimos["PAM"]),    "?", k_optimos["PAM"])
k_cl  <- ifelse(is.na(k_optimos["CLARA"]),  "?", k_optimos["CLARA"])

# =============================================================================
# FUNCIÓN PARA GRAFICAR UN MÉTODO
# =============================================================================
graficar_clustering <- function(df, col_cluster, titulo, subtitulo = "") {
  
  # Verificar que la columna existe
  if (!col_cluster %in% colnames(df)) {
    warning(sprintf("Columna '%s' no encontrada. Saltando gráfico.", col_cluster))
    return(NULL)
  }
  
  n_clusters <- length(unique(na.omit(df[[col_cluster]])))
  
  # Paleta de colores discreta — hasta 20 clusters
  paleta <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
              "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
              "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A")
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[col_cluster]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(
      values = paleta[seq_len(n_clusters)],
      name   = "Cluster",
      na.value = "grey80"
    ) +
    labs(
      title    = titulo,
      subtitle = subtitulo,
      x        = "UMAP 1",
      y        = "UMAP 2"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(color = "gray40", size = 9),
      legend.position = "bottom",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# =============================================================================
# GENERAR GRÁFICOS INDIVIDUALES
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS (K ÓPTIMO) ===\n")

p_kmeans <- graficar_clustering(
  df          = coords_df,
  col_cluster = "Cluster_KMeans",
  titulo      = "K-Means",
  subtitulo   = sprintf("k = %s  |  n = %d árboles", k_km, nrow(coords_df))
)

p_clara <- graficar_clustering(
  df          = coords_df,
  col_cluster = "Cluster_CLARA",
  titulo      = "CLARA",
  subtitulo   = sprintf("k = %s  |  n = %d árboles", k_cl, nrow(coords_df))
)

graficos_individuales <- list(
  list(plot = p_kmeans, nombre = "umap_kmeans"),
  list(plot = p_pam,    nombre = "umap_pam"),
  list(plot = p_clara,  nombre = "umap_clara")
)

for (g in graficos_individuales) {
  if (!is.null(g$plot)) {
    ruta_g <- file.path(DIR_RESULTS, paste0(g$nombre, "_", NOMBRE_BDD, ".png"))
    ggsave(ruta_g, plot = g$plot, width = 8, height = 6, dpi = 300)
    cat("Gráfico guardado:", ruta_g, "\n")
  }
}

# =============================================================================
# GRÁFICO UMAP — MST-KNN SUBDIVISIÓN ITERATIVA
# =============================================================================
if ("Cluster_MSTKNN_Iter" %in% colnames(coords_df)) {
  
  n_iter <- length(unique(na.omit(coords_df$Cluster_MSTKNN_Iter)))
  cat(sprintf("\n=== GENERANDO GRÁFICO MST-KNN ITERATIVA (k=%d) ===\n", n_iter))
  
  # Paleta ampliada — hcl.colors para k grande
  paleta_iter <- if (n_iter <= 20) {
    c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
      "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
      "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
      "#B3B3B3","#1B9E77","#D95F02","#7570B3","#E7298A")[seq_len(n_iter)]
  } else {
    hcl.colors(n_iter, palette = "Dynamic")
  }
  
  p_iter <- ggplot(coords_df,
                   aes(x = UMAP_1, y = UMAP_2, color = Cluster_MSTKNN_Iter)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = paleta_iter,
                       name   = "Cluster",
                       na.value = "grey80") +
    labs(
      title    = "MST-KNN — Subdivisión Iterativa",
      subtitle = sprintf("k = %d  |  n = %d árboles", n_iter, nrow(coords_df)),
      x        = "UMAP 1",
      y        = "UMAP 2"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(color = "gray40", size = 9),
      # ocultar leyenda si hay muchos clusters — satura el gráfico
      legend.position = if (n_iter > 15) "none" else "bottom"
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol = 5   # leyenda en 5 columnas si se muestra
    ))
  
  ruta_iter <- file.path(DIR_RESULTS,
                         paste0("umap_mstknn_iterativa_", NOMBRE_BDD, ".png"))
  ggsave(ruta_iter, plot = p_iter, width = 9, height = 7, dpi = 300)
  cat("Gráfico MST-KNN iterativa guardado:", ruta_iter, "\n")
  
} else {
  cat("Columna 'Cluster_MSTKNN_Iter' no encontrada. Saltando gráfico iterativo.\n")
}

# =============================================================================
# EXPORTAR DATAFRAME COMPLETO A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO COORDENADAS + ETIQUETAS ===\n")

hiperparametros_df <- data.frame(
  Parametro = c("n_neighbors", "min_dist", "n_components", "seed",
                "n_arboles",   "tiempo_calculo_s"),
  Valor     = c(N_NEIGHBORS, MIN_DIST, N_COMPONENTS, SEED,
                nrow(coords_df),
                ifelse(is.na(tiempo_umap), "desde caché", round(tiempo_umap, 1))),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()

wb <- agregar_hoja_formateada(wb, "Coordenadas_y_Clusters",
                              paste0("Coordenadas UMAP + Etiquetas de Clustering — ", NOMBRE_BDD),
                              coords_df,
                              anchos_col = "auto")

wb <- agregar_hoja_formateada(wb, "Hiperparametros_UMAP",
                              "Hiperparámetros UMAP",
                              hiperparametros_df,
                              anchos_col = c(25, 20))

ruta_excel <- file.path(DIR_RESULTS, paste0("umap_clustering_completo_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)

cat("Excel guardado:", ruta_excel, "\n")
cat("\n=== COMPLETADO ===\n")