# =============================================================================
# SCRIPT DE VALIDACIÓN: CONTRASTE ENTRE SESGO DE INJERTO Y ENRIQUECIMIENTO FUNCIONAL
# =============================================================================

library(dplyr)
library(openxlsx)
library(ggplot2)
library(here)
library(rlang)

# 1. Cargar configuraciones y rutas globales
source(here::here("Scripts", "config.R"))

cat("\n=== INICIANDO VALIDACIÓN DE SESGO DE INJERTO VS ENRIQUECIMIENTO ===\n")

# Definir el algoritmo a evaluar (puede ajustarse según necesidad)
algoritmo <- "mstknn"

# 2. Cargar datos
ruta_conteos <- file.path(DIR_RESULTS, paste0("especies_analisis_", NOMBRE_BDD, ".xlsx"))
ruta_clusters <- file.path(DIR_RESULTS, paste0("pca_coordenadas_", NOMBRE_BDD, ".xlsx"))
ruta_enriquecimiento <- file.path(DIR_RESULTS, paste0("enrichment_funcional_", gsub("-","",tolower(algoritmo)), "_", NOMBRE_BDD, ".xlsx"))

# a) Cobertura original por árbol (conteo de hojas/especies antes del injerto)
archivos_bosque <- list.files(DIR_CACHE, pattern = "^bosque_.*\\.rds$", full.names = TRUE)
if (length(archivos_bosque) == 0) stop("No hay bosque en cache.")
bosque_rds <- archivos_bosque[which.min(file.size(archivos_bosque))]
bosque <- readRDS(bosque_rds)
n_hojas_vec <- sapply(bosque, function(t) length(t$tip.label))

df_cobertura <- data.frame(
  Gen = names(n_hojas_vec),
  Cobertura = as.numeric(n_hojas_vec),
  stringsAsFactors = FALSE
)
rm(bosque); gc(verbose = FALSE)

# b) Asignaciones gen -> clúster (leemos desde Genes_Agrupados)
df_asignaciones <- read.xlsx(ruta_clusters, sheet = "Genes_Agrupados", startRow = 2)

# Asegurar que exista la columna 'Cluster' para el join
if (!"Cluster" %in% colnames(df_asignaciones)) {
  if ("Cluster_MSTKNN_Iter" %in% colnames(df_asignaciones)) {
    df_asignaciones <- df_asignaciones %>% rename(Cluster = Cluster_MSTKNN_Iter)
  } else {
    # Tomar la última columna que empiece con Cluster por defecto
    col_cluster <- tail(grep("Cluster", colnames(df_asignaciones), value = TRUE), 1)
    df_asignaciones <- df_asignaciones %>% rename(Cluster = !!sym(col_cluster))
  }
}

# c) Resumen de enriquecimiento por clúster para saber cuáles tienen señal
df_enriquece_resumen <- read.xlsx(ruta_enriquecimiento, sheet = "Resumen", startRow = 2)

# 3. Homogeneizar nombres de columnas dinámicamente para las asignaciones
nombres_asig <- colnames(df_asignaciones)
col_id_asig_gen <- ifelse("Gen" %in% nombres_asig, "Gen", ifelse("Arbol" %in% nombres_asig, "Arbol", nombres_asig[1]))

df_asignaciones <- df_asignaciones %>% rename(Gen = !!sym(col_id_asig_gen))

# 4. Integrar la información a nivel de árbol/gen
analisis_integrado <- df_asignaciones %>%
  left_join(df_cobertura, by = "Gen") %>%
  left_join(df_enriquece_resumen, by = "Cluster") %>%
  mutate(
    Tiene_Enriquecimiento = ifelse(!is.na(N_Terminos_Sig) & N_Terminos_Sig > 0, 
                                   "Con Significado Biológico (Verde)", 
                                   "Sin Enriquecimiento / Ruido (Rojo)")
  ) %>%
  filter(!is.na(Cobertura))


cat("Total de genes procesados en la validación:", nrow(analisis_integrado), "\n")

# 5. Generar Estadísticas Descriptivas del Sesgo
cat("\n--- ESTADÍSTICAS DE VALIDACIÓN ---\n")
print(
  analisis_integrado %>%
    group_by(Tiene_Enriquecimiento) %>%
    summarise(
      Conteo_Genes = n(),
      Media_Cobertura = mean(Cobertura, na.rm = TRUE),
      Mediana_Cobertura = median(Cobertura, na.rm = TRUE)
    )
)

# 6. Calcular test de Wilcoxon sobre la cobertura real (solo si hay al menos 2 grupos)
if (length(unique(analisis_integrado$Tiene_Enriquecimiento)) >= 2) {
  wt <- wilcox.test(Cobertura ~ Tiene_Enriquecimiento, data = analisis_integrado)
  p_val_text <- sprintf("Prueba de Wilcoxon, p-value = %.4e", wt$p.value)
} else {
  p_val_text <- "Prueba de Wilcoxon: No aplicable (solo 1 grupo)"
}
cat("\nResultado Wilcoxon:", p_val_text, "\n")


# 7. Generar Gráfico de Validación (Boxplot)
library(stringr)  # para envolver texto largo

subtitulo_wrap <- str_wrap(
  "Distribución de la cobertura original de los árboles según significancia funcional de su clúster",
  width = 60
)

y_max <- max(analisis_integrado$Cobertura, na.rm = TRUE)

p_val <- ggplot(analisis_integrado, aes(x = Tiene_Enriquecimiento, y = Cobertura, fill = Tiene_Enriquecimiento)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black", outlier.shape = 16) +
  annotate("text", x = 1.5, y = y_max * 1.08, label = p_val_text, size = 4.5) +
  coord_cartesian(ylim = c(0, y_max * 1.18), clip = "off") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Validación de Resiliencia del Pipeline frente al Sesgo de Injerto",
    subtitle = subtitulo_wrap,
    x = "Estado del Clúster",
    y = "Número de Especies (Cobertura Original)",
    fill = "Estado del Clúster"
  ) +
  scale_fill_manual(values = c("Con Significado Biológico (Verde)" = "#2ecc71", 
                               "Sin Enriquecimiento / Ruido (Rojo)" = "#e74c3c")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    plot.margin = margin(t = 15, r = 20, b = 10, l = 10)
  )

# 8. Guardar el gráfico (más ancho para que quepa el título/subtítulo completos)
ruta_plot <- file.path(DIR_RESULTS, "validacion_sesgo_injerto_boxplot.png")
ggsave(ruta_plot, plot = p_val, width = 10, height = 6.5, dpi = 300)
cat("\nGráfico de validación guardado exitosamente en:", ruta_plot, "\n")