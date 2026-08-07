# =============================================================================
# MAPA DE CALOR: CLUSTERS × FUNCIONES BIOLÓGICAS
# Visualización tipo heatmap con dendrogramas jerárquicos
# =============================================================================
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# 1. PARÁMETROS
# =============================================================================
FUENTE_ANALIZAR <- "GO:BP"   # Procesos Biológicos
MIN_CLUSTERS    <- 2         # Término debe aparecer en >= N clusters distintos
TOP_N_TERMINOS  <- 25        # Máximo de términos a mostrar en el mapa

# =============================================================================
# 2. CARGAR DATOS DE ENRIQUECIMIENTO
# =============================================================================
cat("=== CARGANDO DATOS DE ENRIQUECIMIENTO ===\n")

ruta_excel <- file.path(DIR_RESULTS,
                        paste0("enrichment_funcional_", NOMBRE_BDD, ".xlsx"))

if (!file.exists(ruta_excel)) {
  stop("No se encontró el Excel de enriquecimiento: ", ruta_excel,
       "\nEjecuta primero el script 09_enriquecimiento.R")
}

df_todos <- read.xlsx(ruta_excel, sheet = "Todos_los_Terminos", startRow = 2)

cat(sprintf("Registros totales: %d\n", nrow(df_todos)))
cat(sprintf("Fuentes disponibles: %s\n",
            paste(unique(df_todos$Fuente), collapse = ", ")))

# =============================================================================
# 3. FILTRAR Y SELECCIONAR TÉRMINOS COMPARTIDOS
# =============================================================================
cat("\n=== FILTRANDO TÉRMINOS ===\n")

df_filtrado <- df_todos %>%
  filter(Fuente == FUENTE_ANALIZAR) %>%
  mutate(log10_FDR = -log10(as.numeric(FDR)))

# Identificar términos compartidos entre clusters
terminos_compartidos <- df_filtrado %>%
  group_by(Termino) %>%
  summarise(
    N_Clusters = n_distinct(Cluster),
    FDR_Medio  = mean(log10_FDR, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  filter(N_Clusters >= MIN_CLUSTERS) %>%
  arrange(desc(N_Clusters), desc(FDR_Medio)) %>%
  slice_head(n = TOP_N_TERMINOS)

cat(sprintf("Términos compartidos (>= %d clusters): %d\n",
            MIN_CLUSTERS, nrow(terminos_compartidos)))
cat(sprintf("Términos seleccionados (top %d): %d\n",
            TOP_N_TERMINOS, nrow(terminos_compartidos)))

if (nrow(terminos_compartidos) == 0) {
  cat(sprintf("\nNo hay términos funcionales compartidos entre %d o más clusters para la ontología %s.\n", MIN_CLUSTERS, FUENTE_ANALIZAR))
  cat("El análisis de enriquecimiento es tan específico que cada clúster tiene su propia huella biológica única.\n")
  cat("No se puede generar el mapa de calor cruzado, abortando ejecución con éxito.\n")
  quit(save = "no", status = 0)
}

# Filtrar datos a solo esos términos
df_heatmap <- df_filtrado %>%
  filter(Termino %in% terminos_compartidos$Termino)

# =============================================================================
# 4. PIVOTAR A MATRIZ ANCHA (Clusters × Términos)
# =============================================================================
cat("\n=== CONSTRUYENDO MATRIZ ===\n")

# Truncar nombres largos para legibilidad
df_heatmap$Termino_Corto <- str_trunc(df_heatmap$Termino,
                                       width = 55, side = "right")

# Si hay duplicados en Termino_Corto para el mismo cluster, conservar el
# más significativo (menor FDR → mayor -log10)
df_heatmap <- df_heatmap %>%
  group_by(Cluster, Termino_Corto) %>%
  slice_max(log10_FDR, n = 1, with_ties = FALSE) %>%
  ungroup()

# Pivotar: filas = Cluster, columnas = Termino_Corto, valores = -log10(FDR)
matriz_ancha <- df_heatmap %>%
  select(Cluster, Termino_Corto, log10_FDR) %>%
  pivot_wider(
    names_from  = Termino_Corto,
    values_from = log10_FDR,
    values_fill = NA_real_
  ) %>%
  arrange(as.numeric(Cluster))

# Convertir a matriz numérica con nombres de fila
clusters_ids   <- matriz_ancha$Cluster
matriz_numerica <- as.matrix(matriz_ancha[, -1])
rownames(matriz_numerica) <- paste0("Cluster ", clusters_ids)

cat(sprintf("Matriz resultante: %d clusters × %d términos\n",
            nrow(matriz_numerica), ncol(matriz_numerica)))
cat(sprintf("Celdas con dato: %d  |  Celdas NA: %d\n",
            sum(!is.na(matriz_numerica)), sum(is.na(matriz_numerica))))

# Guardar matriz con NA para visualización
matriz_visual <- matriz_numerica

# Reemplazar NA por 0: sin enriquecimiento = -log10(1) = 0 (sin significancia)
# Necesario porque hclust no puede calcular distancias con NA
matriz_numerica[is.na(matriz_numerica)] <- 0

# =============================================================================
# 5. GENERAR MAPA DE CALOR
# =============================================================================
cat("\n=== GENERANDO MAPA DE CALOR ===\n")

# Paleta de colores: plasma (coherente con gráfico de burbujas)
n_colores  <- 100
paleta_col <- viridisLite::plasma(n_colores, direction = 1)

# Dimensiones adaptativas según tamaño de la matriz
alto_fig  <- max(10, nrow(matriz_numerica) * 0.35 + 6)
ancho_fig <- max(16, ncol(matriz_numerica) * 0.60 + 8)

# Tamaño de fuente adaptativo
fontsize_fila <- if (nrow(matriz_numerica) > 80) 5 else
                 if (nrow(matriz_numerica) > 40) 7 else 9
fontsize_col  <- if (ncol(matriz_numerica) > 20) 8 else 10

# Ruta de salida
ruta_heatmap <- file.path(DIR_RESULTS,
                          paste0("Mapa_Calor_Funcional_", NOMBRE_BDD, ".png"))

# Dendrogramas calculados sobre matriz_numerica (sin NA)
dend_rows <- hclust(dist(t(matriz_numerica), method = "euclidean"), method = "ward.D2")  # términos en filas
dend_cols <- hclust(dist(matriz_numerica,    method = "euclidean"), method = "ward.D2")  # clusters en cols


# Generar y guardar heatmap
pheatmap(
  mat               = t(matriz_visual),
  
  margins           = c(18, 35),
  
  # Dendrogramas — clustering jerárquico en ambos ejes
  cluster_rows      = dend_rows,
  cluster_cols      = dend_cols,
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  
  # Escala de color
  color             = paleta_col,
  na_col            = "grey90",        # celdas sin dato en gris claro
  
  # Etiquetas
  fontsize_row      = fontsize_fila,
  fontsize_col      = fontsize_col,
  angle_col         = 45,              # rotar nombres de columna
  
  # Leyenda
  legend             = TRUE,
  legend_breaks      = pretty(range(matriz_numerica, na.rm = TRUE), n = 5),
  legend_labels      = pretty(range(matriz_numerica, na.rm = TRUE), n = 5),
  
  # Bordes
  border_color      = "grey80",
  
  # Título
  main              = paste0("Convergencia Funcional entre Clusters — ",
                             FUENTE_ANALIZAR, " (", NOMBRE_BDD, ")"),
  
  # Guardar directamente a archivo
  filename          = ruta_heatmap,
  width             = ancho_fig,
  height            = alto_fig,
  dpi               = 300
)

cat("Mapa de calor guardado:", ruta_heatmap, "\n")

# =============================================================================
# 6. RESUMEN
# =============================================================================
cat("\n=== RESUMEN ===\n")
cat(sprintf("Fuente analizada    : %s\n", FUENTE_ANALIZAR))
cat(sprintf("Min. clusters/término: %d\n", MIN_CLUSTERS))
cat(sprintf("Términos graficados : %d\n", ncol(matriz_numerica)))
cat(sprintf("Clusters graficados : %d\n", nrow(matriz_numerica)))
cat(sprintf("Rango -log10(FDR)   : [%.2f, %.2f]\n",
            min(matriz_numerica, na.rm = TRUE),
            max(matriz_numerica, na.rm = TRUE)))
cat(sprintf("Archivo             : %s\n", ruta_heatmap))
cat("\n=== COMPLETADO ===\n")
