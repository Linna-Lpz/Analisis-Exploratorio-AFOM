# =============================================================================
# 05d_reduccion_dimensional_PCA.R
# PCA sobre la matriz RF normalizada — con caché .rds
# Salidas: coordenadas en Excel, gráfico ggplot2, scree plot, caché .rds
# =============================================================================

library(ggplot2)
library(openxlsx)
library(patchwork)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# CARGAR MATRIZ RF DESDE CACHÉ
# =============================================================================
cat("=== CARGANDO MATRIZ RF ===\n")

ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché: ", ruta_cache_matriz,
       "\nEjecuta primero el script 04_calcular_matriz.R")
}

matriz_cuadrada <- readRDS(ruta_cache_matriz)
cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# CALCULAR PCA — con caché condicional
# =============================================================================
cat("\n=== CALCULANDO PCA ===\n")

ruta_cache_pca <- file.path(DIR_CACHE, "pca_resultado.rds")

# Número de componentes a retener para exportar coordenadas
# (guarda 10 para tener flexibilidad, se grafica solo PC1 y PC2)
N_COMPONENTS_GUARDAR <- 10

if (file.exists(ruta_cache_pca)) {
  cat("PCA encontrado en caché. Cargando...\n")
  tiempo_pca <- system.time({
    pca_resultado <- readRDS(ruta_cache_pca)
  })
  cat("Cargado desde caché.\n")
  
} else {
  cat("Calculando PCA...\n")
  tiempo_pca <- system.time({
    
    # PCA sobre la matriz de distancias RF
    # center = TRUE : resta la media de cada columna (estándar en PCA)
    # scale  = FALSE: no escala por SD — las columnas ya están en la misma
    #                 unidad [0,1] al ser distancias RF normalizadas
    pca_resultado <- prcomp(
      x      = matriz_cuadrada,
      center = TRUE,
      scale. = FALSE
    )
    saveRDS(pca_resultado, file = ruta_cache_pca)
  })
  cat("PCA calculado y guardado en caché.\n")
}

cat(sprintf("Tiempo PCA: %.1f segundos\n", tiempo_pca["elapsed"]))

# =============================================================================
# PROPORCIÓN VARIANZA EXPLICADA (PVE)
# =============================================================================
cat("\n=== VARIANZA EXPLICADA ===\n")

varianza_prop    <- summary(pca_resultado)$importance["Proportion of Variance", ]
varianza_acum    <- summary(pca_resultado)$importance["Cumulative Proportion", ]
desv_std         <- summary(pca_resultado)$importance["Standard deviation", ]

var_pc1 <- round(varianza_prop["PC1"] * 100, 2)
var_pc2 <- round(varianza_prop["PC2"] * 100, 2)
var_total_2d <- round(varianza_acum["PC2"] * 100, 2)

# Número de componentes necesarios para alcanzar umbrales de varianza
n_comp_70 <- which(varianza_acum >= 0.70)[1]
n_comp_80 <- which(varianza_acum >= 0.80)[1]
n_comp_90 <- which(varianza_acum >= 0.90)[1]

cat(sprintf("PC1: %.2f%%  |  PC2: %.2f%%  |  Total 2D: %.2f%%\n",
            var_pc1, var_pc2, var_total_2d))
cat(sprintf("Componentes para 70%% varianza: %d\n", n_comp_70))
cat(sprintf("Componentes para 80%% varianza: %d\n", n_comp_80))
cat(sprintf("Componentes para 90%% varianza: %d\n", n_comp_90))

# =============================================================================
# CONSTRUIR DATAFRAME DE COORDENADAS (primeras 10 PCs)
# =============================================================================
n_comp_real <- min(N_COMPONENTS_GUARDAR, ncol(pca_resultado$x))

coords_matrix <- pca_resultado$x[, 1:n_comp_real]
colnames(coords_matrix) <- paste0("PC", 1:n_comp_real)

coords_df <- data.frame(
  Arbol = rownames(matriz_cuadrada),
  as.data.frame(coords_matrix),
  stringsAsFactors = FALSE
)

cat("\nPrimeras filas de coordenadas (PC1 y PC2):\n")
print(head(coords_df[, 1:3], 5))

# =============================================================================
# LEER Y UNIR ETIQUETAS DE CLUSTERING
# =============================================================================
resultado    <- unir_etiquetas_clustering(coords_df, DIR_RESULTS)
coords_df    <- resultado$coords_df
k_optimos    <- resultado$k_optimos

# Extraer k por método para subtítulos
k_km  <- ifelse(is.na(k_optimos["KMeans"]), "?", k_optimos["KMeans"])
k_pam <- ifelse(is.na(k_optimos["PAM"]),    "?", k_optimos["PAM"])
k_cl  <- ifelse(is.na(k_optimos["CLARA"]),  "?", k_optimos["CLARA"])


# =============================================================================
# GRÁFICO PRINCIPAL — PC1 vs PC2
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS ===\n")

graficar_clustering_pca <- function(df, col_cluster, titulo, subtitulo = "") {
  if (!col_cluster %in% colnames(df)) {
    warning(sprintf("Columna '%s' no encontrada. Saltando gráfico.", col_cluster))
    return(NULL)
  }
  n_clusters <- length(unique(na.omit(df[[col_cluster]])))
  paleta <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
              "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
  
  ggplot(df, aes(x = PC1, y = PC2, color = .data[[col_cluster]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = paleta[seq_len(n_clusters)],
                       name = "Cluster", na.value = "grey80") +
    labs(
      title    = titulo,
      subtitle = subtitulo,
      x        = paste0("PC1 (", var_pc1, "%)"),
      y        = paste0("PC2 (", var_pc2, "%)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(color = "gray40", size = 9),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

p_kmeans <- graficar_clustering_pca(coords_df, "Cluster_KMeans", "K-Means",
                                    sprintf("k = %s  |  n = %d árboles  |  PC1+PC2 = %.1f%%", k_km,  nrow(coords_df), var_total_2d))
p_pam    <- graficar_clustering_pca(coords_df, "Cluster_PAM",    "PAM (K-Medoids)",
                                    sprintf("k = %s  |  n = %d árboles  |  PC1+PC2 = %.1f%%", k_pam, nrow(coords_df), var_total_2d))
p_clara  <- graficar_clustering_pca(coords_df, "Cluster_CLARA",  "CLARA",
                                    sprintf("k = %s  |  n = %d árboles  |  PC1+PC2 = %.1f%%", k_cl,  nrow(coords_df), var_total_2d))

plots_disponibles <- Filter(Negate(is.null), list(p_kmeans, p_pam, p_clara))

panel_comparativo <- wrap_plots(plots_disponibles, ncol = min(length(plots_disponibles), 3)) +
  plot_annotation(
    title   = paste0("Comparación de Algoritmos de Clustering — PCA (", NOMBRE_BDD, ")"),
    caption = paste0("Distancia Robinson-Foulds normalizada  |  prcomp(center=TRUE, scale=FALSE)"),
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.caption = element_text(color = "gray50", size = 8, hjust = 0.5)
    )
  )

ruta_grafico <- file.path(DIR_RESULTS, paste0("pca_comparativo_clustering_", NOMBRE_BDD, ".png"))
ggsave(ruta_grafico,
       plot   = panel_comparativo,
       width  = 7 * length(plots_disponibles),
       height = 6,
       dpi    = 300)
cat("Panel comparativo PCA guardado:", ruta_grafico, "\n")

# =============================================================================
# SCREE PLOT — varianza explicada por componente
# Herramienta clave para decidir cuántos componentes usar en k-means
# =============================================================================
n_scree <- min(30, length(varianza_prop))  # mostrar hasta 30 PCs

scree_df <- data.frame(
  componente     = paste0("PC", 1:n_scree),
  varianza_pct   = as.numeric(varianza_prop[1:n_scree]) * 100,
  varianza_acum  = as.numeric(varianza_acum[1:n_scree]) * 100,
  stringsAsFactors = FALSE
)
scree_df$componente <- factor(scree_df$componente,
                              levels = scree_df$componente)

p_scree <- ggplot(scree_df, aes(x = componente)) +
  
  # Barras: varianza individual por componente
  geom_col(aes(y = varianza_pct),
           fill = "mediumpurple", alpha = 0.7) +
  
  # Línea: varianza acumulada
  geom_line(aes(y = varianza_acum, group = 1),
            color = "firebrick", linewidth = 0.8) +
  geom_point(aes(y = varianza_acum),
             color = "firebrick", size = 2) +
  
  # Líneas de referencia en 70%, 80%, 90%
  geom_hline(yintercept = c(70, 80, 90),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  annotate("text", x = n_scree, y = c(71, 81, 91),
           label = c("70%", "80%", "90%"),
           hjust = 1, size = 3, color = "gray40") +
  
  scale_y_continuous(
    name     = "Varianza individual (%)",
    limits   = c(0, 100),
    sec.axis = sec_axis(~ ., name = "Varianza acumulada (%)")
  ) +
  labs(
    title    = paste0("Scree Plot PCA — ", NOMBRE_BDD),
    subtitle = paste0(
      "Componentes para 70%: ", n_comp_70,
      "  |  80%: ", n_comp_80,
      "  |  90%: ", n_comp_90
    ),
    x       = "Componente principal",
    caption = "Barras: varianza individual  |  Línea roja: varianza acumulada"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 8)
  )

ruta_scree <- file.path(DIR_RESULTS, paste0("pca_scree_", NOMBRE_BDD, ".png"))
ggsave(ruta_scree, plot = p_scree, width = 10, height = 6, dpi = 300)
cat("Scree plot guardado:", ruta_scree, "\n")

# =============================================================================
# EXPORTAR COORDENADAS Y DIAGNÓSTICO A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

# Tabla de varianza por componente (todas las PCs)
n_total_comp <- length(varianza_prop)
varianza_df <- data.frame(
  Componente        = paste0("PC", 1:n_total_comp),
  Desv_std          = round(as.numeric(desv_std), 6),
  Varianza_pct      = round(as.numeric(varianza_prop) * 100, 4),
  Varianza_acum_pct = round(as.numeric(varianza_acum) * 100, 4),
  stringsAsFactors  = FALSE
)

# Tabla de parámetros
parametros_df <- data.frame(
  Parametro = c("metodo", "funcion_R", "center", "scale",
                "n_componentes_guardados", "var_PC1_pct", "var_PC2_pct",
                "var_total_2D_pct", "n_comp_para_70pct",
                "n_comp_para_80pct", "n_comp_para_90pct",
                "n_arboles", "tiempo_calculo_s"),
  Valor     = c("PCA", "prcomp", "TRUE", "FALSE",
                n_comp_real, var_pc1, var_pc2,
                var_total_2d, n_comp_70,
                n_comp_80, n_comp_90,
                nrow(coords_df), round(tiempo_pca["elapsed"], 1)),
  stringsAsFactors = FALSE
)

# Anchos dinámicos para coordenadas según n_comp_real
anchos_coords <- c(40, rep(16, n_comp_real))

wb <- createWorkbook()

agregar_hoja_formateada(wb, "Coordenadas_PCA",
                        paste0("Coordenadas PCA — ", NOMBRE_BDD),
                        coords_df,
                        anchos_col = anchos_coords)

agregar_hoja_formateada(wb, "Varianza_Explicada",
                        "Varianza Explicada por Componente",
                        varianza_df,
                        anchos_col = c(15, 18, 18, 22))

agregar_hoja_formateada(wb, "Diagnostico",
                        "Diagnóstico y Parámetros PCA",
                        parametros_df,
                        anchos_col = c(30, 25))

ruta_excel <- file.path(DIR_RESULTS, paste0("pca_coordenadas_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Coordenadas guardadas:", ruta_excel, "\n")

cat("\n=== COMPLETADO ===\n")