# =============================================================================
# SCRIPT DE VALIDACIÓN: CONTRASTE ENTRE SESGO DE INJERTO Y ENRIQUECIMIENTO FUNCIONAL
# =============================================================================

library(dplyr)
library(openxlsx)
library(ggplot2)
library(here)

# 1. Cargar configuraciones y rutas globales
source(here::here("Scripts", "config.R"))

cat("\n=== INICIANDO VALIDACIÓN DE SESGO DE INJERTO VS ENRIQUECIMIENTO ===\n")

# 2. Cargar el conteo de especies por árbol original desde la fuente o archivos exportados
# Nota: Asumiendo que puedes leer el resumen de conteos de especies por árbol
ruta_conteos <- file.path(DIR_RESULTS, "ESPECIES", "especies_analisis_OrthoMaM.xlsx")

# O alternativamente, mapear desde la base de datos de coordenadas o genes agrupados
# Vamos a cargar los datos de genes por clúster y el resumen de enriquecimiento
ruta_clusters <- file.path(DIR_RESULTS, "hgnc_listas_por_cluster_OrthoMaM.xlsx")
ruta_enriquecimiento <- file.path(DIR_RESULTS, "enrichment_funcional_OrthoMaM.xlsx")

# Leer resumen de clusters y resumen de enriquecimiento
df_cluster_resumen <- read.xlsx(ruta_clusters, sheet = "Resumen", startRow = 2) # Contiene N_Arboles, N_Con_HGNC, etc.
df_enriquece_resumen <- read.xlsx(ruta_enriquecimiento, sheet = "Resumen", startRow = 2) # Contiene N_Terminos_Sig

# 3. Cruzar la información para evaluar la hipótesis crítica
# Hipótesis: Los clústeres con baja cantidad de especies originales (fuertemente imputados) 
# tienen menor cantidad de términos significativos (FDR < 0.05).

# Si tienes un archivo consolidado o por cluster, unimos el conteo de especies originales por árbol/gen
# Vamos a estructurar una tabla analítica por clúster:
analisis_integrado <- df_cluster_resumen %>%
  left_join(df_enriquece_resumen, by = "Cluster") %>%
  mutate(
    Tiene_Enriquecimiento = ifelse(N_Terminos_Sig > 0, "Con Significado Biológico (Verde)", "Sin Enriquecimiento / Ruido (Rojo)")
  )

cat("Resumen de clústeres procesados:", nrow(analisis_integrado), "\n")

# 4. Generar Estadísticas Descriptivas del Sesgo
cat("\n--- ESTADÍSTICAS DE VALIDACIÓN ---\n")
print(
  analisis_integrado %>%
    group_by(Tiene_Enriquecimiento) %>%
    summarise(
      Conteo_Clusters = n(),
      Media_Arboles = mean(N_Arboles, na.rm = TRUE),
      Mediana_Arboles = median(N_Arboles, na.rm = TRUE)
    )
)

# 5. Generar Gráfico de Validación (Boxplot / Violin Plot) para la Tesis
p_val <- ggplot(analisis_integrado, aes(x = Tiene_Enriquecimiento, y = N_Arboles, fill = Tiene_Enriquecimiento)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black", outlier.shape = 16) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Validación de Resiliencia del Pipeline frente al Sesgo de Injerto",
    subtitle = "Distribución del tamaño del árbol original (cobertura) según significancia funcional",
    x = "Categoría del Clúster",
    y = "Número de Especies / Cobertura Original del Árbol",
    fill = "Estado del Clúster"
  ) +
  scale_fill_manual(values = c("Con Significado Biológico (Verde)" = "#2ecc71", 
                               "Sin Enriquecimiento / Ruido (Rojo)" = "#e74c3c")) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# Guardar el gráfico para el documento de tesis
ggsave(file.path(DIR_RESULTS, "validacion_sesgo_injerto_boxplot.png"), plot = p_val, width = 8, height = 6, dpi = 300)
cat("\nGráfico de validación guardado exitosamente en:", file.path(DIR_RESULTS, "validacion_sesgo_injerto_boxplot.png"), "\n")