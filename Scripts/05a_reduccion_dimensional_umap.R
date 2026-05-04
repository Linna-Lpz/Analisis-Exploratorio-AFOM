# =============================================================================
# 05_reduccion_dimensional.R
# UMAP sobre la matriz RF normalizada — con caché .rds
# Salidas: coordenadas en Excel, gráfico ggplot2, caché .rds
# =============================================================================

library(uwot)
library(ggplot2)
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

# =============================================================================
# CALCULAR UMAP — con caché condicional
# =============================================================================
cat("\n=== CALCULANDO UMAP ===\n")

ruta_cache_umap <- file.path(DIR_CACHE, "umap2_coords.rds")

# Hiperparámetros UMAP — ajustables según el tamaño del dataset
N_NEIGHBORS  <- 15    # Vecinos locales: más alto = estructura más global
MIN_DIST     <- 0.1   # Compactación de clusters: más bajo = clusters más densos
N_COMPONENTS <- 2     # Dimensiones de salida (2 para visualización y k-means)
SEED         <- 42

if (file.exists(ruta_cache_umap)) {
  cat("UMAP encontrado en caché. Cargando...\n")
  tiempo_umap <- system.time({
    umap_coords <- readRDS(ruta_cache_umap)
  })
} else {
  cat("Calculando UMAP (puede tardar varios minutos)...\n")
  tiempo_umap <- system.time({
    set.seed(SEED)
    umap_coords <- umap(
      X            = as.dist(matriz_cuadrada),   # <- conversión clave
      n_neighbors  = N_NEIGHBORS,
      min_dist     = MIN_DIST,
      n_components = N_COMPONENTS,
      n_threads   = parallel::detectCores() - 1,  # paralelizar
      verbose      = TRUE
    )
    saveRDS(umap_coords, file = ruta_cache_umap)
  })
  cat("UMAP calculado y guardado en caché.\n")
}

cat("Tiempo UMAP:", round(tiempo_umap["elapsed"], 1), "segundos\n")
cat("Dimensiones resultado:", nrow(umap_coords), "x", ncol(umap_coords), "\n")

# =============================================================================
# CONSTRUIR DATAFRAME DE COORDENADAS
# =============================================================================
coords_df <- data.frame(
  nombre_arbol = rownames(matriz_cuadrada),
  UMAP_1       = umap_coords[, 1],
  UMAP_2       = umap_coords[, 2],
  stringsAsFactors = FALSE
)

cat("\nPrimeras filas de coordenadas:\n")
print(head(coords_df, 5))

# =============================================================================
# GRÁFICO ESTÁTICO CON ggplot2
# =============================================================================
cat("\n=== GENERANDO GRÁFICO ===\n")

p <- ggplot(coords_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "steelblue",
    alpha = 0.6,
    size  = 1.5
  ) +
  labs(
    title    = paste0("UMAP — Matriz RF Normalizada (", NOMBRE_BDD, ")"),
    subtitle = paste0(
      "n_neighbors = ", N_NEIGHBORS,
      "  |  min_dist = ", MIN_DIST,
      "  |  n = ", nrow(coords_df), " árboles"
    ),
    x = "UMAP 1",
    y = "UMAP 2",
    caption = "Distancia Robinson-Foulds normalizada"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8)
  )

ruta_grafico <- file.path(DIR_RESULTS, paste0("umap_", NOMBRE_BDD, ".png"))
ggsave(ruta_grafico,
       plot   = p,
       width  = 10,
       height = 7,
       dpi    = 300)

cat("Gráfico guardado:", ruta_grafico, "\n")

# =============================================================================
# EXPORTAR COORDENADAS A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO COORDENADAS ===\n")

# Resumen de hiperparámetros para documentar en el Excel
hiperparametros_df <- data.frame(
  Parametro = c("n_neighbors", "min_dist", "n_components", "input", "seed",
                "n_arboles", "tiempo_calculo_s"),
  Valor     = c(N_NEIGHBORS, MIN_DIST, N_COMPONENTS, "dist (as.dist)", SEED,
                nrow(coords_df), round(tiempo_umap["elapsed"], 1)),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()

agregar_hoja_formateada(wb, "Coordenadas_UMAP",
                        paste0("Coordenadas UMAP — ", NOMBRE_BDD),
                        coords_df,
                        anchos_col = c(40, 18, 18))

agregar_hoja_formateada(wb, "Hiperparametros",
                        "Hiperparámetros UMAP",
                        hiperparametros_df,
                        anchos_col = c(25, 20))

ruta_excel <- file.path(DIR_RESULTS, paste0("umap_coordenadas_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)

cat("Coordenadas guardadas:", ruta_excel, "\n")

cat("\n=== COMPLETADO ===\n")