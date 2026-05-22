# =============================================================================
# 05b_reduccion_dimensional_nMDS.R
# nMDS (no métrico) sobre la matriz RF normalizada — con caché .rds
# Salidas: coordenadas en Excel, gráfico ggplot2, Shepard plot, caché .rds
# =============================================================================

library(vegan)
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
dist_rf         <- as.dist(matriz_cuadrada)

cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# CALCULAR nMDS — con caché condicional
# =============================================================================
cat("\n=== CALCULANDO nMDS ===\n")

ruta_cache_nmds <- file.path(DIR_CACHE, "nmds_coords.rds")

# --- Hiperparámetros ---
N_COMPONENTS <- 3      # dimensiones de salida
TRYMAX       <- 20     # máximo de arranques aleatorios
MAXIT        <- 500    # máximo de iteraciones por arranque
SEED         <- 42

# Referencia de stress (Kruskal 1964)
interpretar_stress <- function(s) {
  dplyr::case_when(
    s < 0.05 ~ "Excelente",
    s < 0.10 ~ "Bueno",
    s < 0.20 ~ "Aceptable",
    TRUE     ~ "Deficiente — considerar más dimensiones"
  )
}

if (file.exists(ruta_cache_nmds)) {
  cat("nMDS encontrado en caché. Cargando...\n")
  tiempo_nmds <- system.time({
    nmds_resultado <- readRDS(ruta_cache_nmds)
  })
  cat("Cargado desde caché.\n")
  
} else {
  cat("Calculando nMDS (puede tardar varios minutos con", nrow(matriz_cuadrada), "árboles)...\n")
  cat("trymax =", TRYMAX, "| maxit =", MAXIT, "| k =", N_COMPONENTS, "\n\n")
  
  tiempo_nmds <- system.time({
    set.seed(SEED)
    nmds_resultado <- metaMDS(
      comm          = dist_rf,   # matriz de distancias precomputada
      k             = N_COMPONENTS,
      trymax        = TRYMAX,
      maxit         = MAXIT,
      autotransform = FALSE,     # desactivar: no son datos de comunidad
      wascores      = FALSE,     # desactivar: no hay scores de especies
      noshare       = FALSE,     # desactivar: stepacross no aplica con dist precomputada
      trace         = 1          # mostrar progreso en consola (0 = silencioso)
    )
    saveRDS(nmds_resultado, file = ruta_cache_nmds)
  })
  cat("\nnMDS calculado y guardado en caché.\n")
}

# =============================================================================
# EVALUAR CALIDAD — STRESS
# =============================================================================
cat("\n=== EVALUACIÓN DE CALIDAD ===\n")

stress_val    <- nmds_resultado$stress
interpretacion <- interpretar_stress(stress_val)

cat(sprintf("Stress final  : %.4f\n", stress_val))
cat(sprintf("Interpretación: %s\n",   interpretacion))
cat(sprintf("Convergencia  : %s\n",
            ifelse(nmds_resultado$converged, "SÍ", "NO — aumentar trymax o maxit")))
cat(sprintf("Tiempo        : %.1f segundos\n", tiempo_nmds["elapsed"]))

# =============================================================================
# CONSTRUIR DATAFRAME DE COORDENADAS
# =============================================================================
coords_nmds <- scores(nmds_resultado, display = "sites")

coords_df <- data.frame(
  nombre_arbol = rownames(matriz_cuadrada),
  NMDS_1       = coords_nmds[, 1],
  NMDS_2       = coords_nmds[, 2],
  stringsAsFactors = FALSE
)

cat("\nPrimeras filas de coordenadas:\n")
print(head(coords_df, 5))

# =============================================================================
# GRÁFICO PRINCIPAL — distribución de árboles
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS ===\n")

p_nmds <- ggplot(coords_df, aes(x = NMDS_1, y = NMDS_2)) +
  geom_point(
    color = "steelblue",
    alpha = 0.6,
    size  = 1.5
  ) +
  labs(
    title    = paste0("nMDS — Matriz RF Normalizada (", NOMBRE_BDD, ")"),
    subtitle = paste0(
      "Stress = ", round(stress_val, 4),
      "  (", interpretacion, ")",
      "  |  n = ", nrow(coords_df), " árboles"
    ),
    x       = "nMDS 1",
    y       = "nMDS 2",
    caption = "Distancia Robinson-Foulds normalizada | metaMDS (vegan)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8)
  )

ruta_grafico <- file.path(DIR_RESULTS, paste0("nmds_", NOMBRE_BDD, ".png"))
ggsave(ruta_grafico, plot = p_nmds, width = 10, height = 7, dpi = 300)
cat("Gráfico nMDS guardado:", ruta_grafico, "\n")

# =============================================================================
# SHEPARD PLOT — diagnóstico de ajuste
# Compara distancias originales RF vs distancias en espacio nMDS
# Un buen ajuste = puntos cercanos a la línea escalonada
# =============================================================================
shepard_data <- data.frame(
  dist_original = as.vector(dist_rf),
  dist_nmds     = as.vector(dist(scores(nmds_resultado)))
)

p_shepard <- ggplot(shepard_data, aes(x = dist_original, y = dist_nmds)) +
  geom_point(alpha = 0.1, size = 0.8, color = "gray40") +
  geom_smooth(method = "loess", se = FALSE,
              color = "firebrick", linewidth = 0.8) +
  labs(
    title    = paste0("Shepard Plot — nMDS (", NOMBRE_BDD, ")"),
    subtitle = paste0("Stress = ", round(stress_val, 4),
                      "  |  Ajuste ideal = relación monótona"),
    x       = "Distancia RF original",
    y       = "Distancia en espacio nMDS",
    caption = "Línea roja: ajuste LOESS — dispersión baja indica buen ajuste"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8)
  )

ruta_shepard <- file.path(DIR_RESULTS, paste0("nmds_shepard_", NOMBRE_BDD, ".png"))
ggsave(ruta_shepard, plot = p_shepard, width = 8, height = 6, dpi = 300)
cat("Shepard plot guardado:", ruta_shepard, "\n")

# =============================================================================
# EXPORTAR COORDENADAS Y DIAGNÓSTICO A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

parametros_df <- data.frame(
  Parametro = c("metodo", "engine", "k", "trymax", "maxit", "seed",
                "stress", "interpretacion_stress", "convergencia",
                "n_arboles", "tiempo_calculo_s"),
  Valor     = c("nMDS no métrico", "monoMDS (vegan)", N_COMPONENTS,
                TRYMAX, MAXIT, SEED,
                round(stress_val, 4), interpretacion,
                ifelse(nmds_resultado$converged, "Sí", "No"),
                nrow(coords_df), round(tiempo_nmds["elapsed"], 1)),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()

agregar_hoja_formateada(wb, "Coordenadas_nMDS",
                        paste0("Coordenadas nMDS — ", NOMBRE_BDD),
                        coords_df,
                        anchos_col = c(40, 18, 18))

agregar_hoja_formateada(wb, "Diagnostico",
                        "Diagnóstico y Parámetros nMDS",
                        parametros_df,
                        anchos_col = c(30, 35))

ruta_excel <- file.path(DIR_RESULTS, paste0("nmds_coordenadas_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Coordenadas guardadas:", ruta_excel, "\n")

cat("\n=== COMPLETADO ===\n")