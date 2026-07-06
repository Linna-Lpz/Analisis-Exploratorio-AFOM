# =============================================================================
# ANÁLISIS CRUZADO: TÉRMINOS COMPARTIDOS ENTRE CLUSTERS
# =============================================================================
library(openxlsx)
library(ggplot2)
library(dplyr)
library(stringr)

# 1. Cargar los datos
ruta_excel <- file.path(DIR_RESULTS,
                            paste0("enrichment_funcional_", NOMBRE_BDD, ".xlsx"))
df_todos <- read.xlsx(ruta_excel, sheet = "Todos_los_Terminos", startRow = 2)

# 2. Parámetros de filtrado
FUENTE_ANALIZAR <- "GO:BP" # Nos enfocamos en Procesos Biológicos
MIN_CLUSTERS    <- 2       # El término debe aparecer en al menos 2 clusters distintos
TOP_N_TERMINOS  <- 25      # Máximo de términos a graficar (para que no colapse el gráfico)

# 3. Filtrar y procesar datos
df_filtrado <- df_todos %>%
  filter(Fuente == FUENTE_ANALIZAR) %>%
  mutate(log10_FDR = -log10(as.numeric(FDR)))

# Encontrar los términos compartidos (que aparecen en >= MIN_CLUSTERS)
terminos_compartidos <- df_filtrado %>%
  group_by(Termino) %>%
  summarise(
    N_Clusters = n_distinct(Cluster),
    FDR_Medio  = mean(log10_FDR, na.rm = TRUE)
  ) %>%
  filter(N_Clusters >= MIN_CLUSTERS) %>%
  arrange(desc(N_Clusters), desc(FDR_Medio)) %>%
  slice_head(n = TOP_N_TERMINOS) # Nos quedamos con los más compartidos y significativos

# Quedarnos solo con los datos de esos términos top
df_plot <- df_filtrado %>%
  filter(Termino %in% terminos_compartidos$Termino)

# Acortar nombres muy largos para el gráfico
df_plot$Termino_Corto <- str_trunc(df_plot$Termino, width = 60, side = "right")

# Ordenar los términos alfabéticamente
niveles_alfa <- sort(unique(str_trunc(df_plot$Termino, width = 60, side = "right")), 
                     decreasing = TRUE)
df_plot$Termino_Corto <- factor(df_plot$Termino_Corto, levels = niveles_alfa)

# Asegurar que los clusters se grafiquen como factores discretos (1, 2, 3...)
df_plot$Cluster <- factor(df_plot$Cluster, levels = sort(unique(as.numeric(df_plot$Cluster))))

# 4. Generar el Gráfico de Burbujas
p_cruzado <- ggplot(df_plot, aes(x = Cluster, y = Termino_Corto)) +
  geom_point(aes(size = as.numeric(N_Genes_en_Termino), fill = log10_FDR), 
             shape = 21, alpha = 0.8, color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "-log10(FDR)\n(Significancia)") +
  scale_size_continuous(range = c(2, 8), name = "Genes en\nel término") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Convergencia Funcional entre Clusters Evolutivos",
    subtitle = paste("Términos de", FUENTE_ANALIZAR, "presentes en múltiples clusters"),
    x = "ID del Cluster",
    y = "Proceso Biológico (GO:BP)"
  )

# 5. Mostrar y guardar
ruta_plot <- file.path(DIR_RESULTS, paste0("Grafico_Burbujas_Cruzado_",NOMBRE_BDD, ".png"))
ggsave(ruta_plot, plot = p_cruzado, width = 12, height = 8, dpi = 300)