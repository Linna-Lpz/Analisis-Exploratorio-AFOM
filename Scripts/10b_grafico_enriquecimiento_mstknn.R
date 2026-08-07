# =============================================================================
# 10b_grafico_enriquecimiento_mstknn.R
# ANÁLISIS CRUZADO: TÉRMINOS COMPARTIDOS ENTRE CLUSTERS MST-kNN
# =============================================================================
library(openxlsx)
library(ggplot2)
library(dplyr)
library(stringr)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# 1. Cargar datos
ruta_excel <- file.path(DIR_RESULTS, paste0("enrichment_funcional_mstknn_", NOMBRE_BDD, ".xlsx"))
if (!file.exists(ruta_excel)) stop("No se encontro ", ruta_excel)
df_todos <- read.xlsx(ruta_excel, sheet = "Todos_los_Terminos", startRow = 2)

# 2. Parámetros
FUENTE_ANALIZAR <- "GO:BP"
MIN_CLUSTERS    <- 2
TOP_N_TERMINOS  <- 25

# 3. Filtrar
df_filtrado <- df_todos %>%
  filter(Fuente == FUENTE_ANALIZAR) %>%
  mutate(log10_FDR = -log10(as.numeric(FDR)))

terminos_compartidos <- df_filtrado %>%
  group_by(Termino) %>%
  summarise(
    N_Clusters = n_distinct(Cluster),
    FDR_Medio  = mean(log10_FDR, na.rm = TRUE)
  ) %>%
  filter(N_Clusters >= MIN_CLUSTERS) %>%
  arrange(desc(N_Clusters), desc(FDR_Medio)) %>%
  slice_head(n = TOP_N_TERMINOS)

if (nrow(terminos_compartidos) == 0) {
  cat("No hay términos compartidos. Abortando gráfico cruzado.\n")
  quit(save = "no", status = 0)
}

df_plot <- df_filtrado %>%
  filter(Termino %in% terminos_compartidos$Termino) %>%
  mutate(Termino_Corto = str_trunc(Termino, width = 60, side = "right"))

niveles_alfa <- sort(unique(df_plot$Termino_Corto), decreasing = TRUE)
df_plot$Termino_Corto <- factor(df_plot$Termino_Corto, levels = niveles_alfa)
df_plot$Cluster <- factor(df_plot$Cluster, levels = sort(unique(as.numeric(df_plot$Cluster))))

# 4. Gráfico
p_cruzado <- ggplot(df_plot, aes(x = Cluster, y = Termino_Corto)) +
  geom_point(aes(size = as.numeric(N_Genes_en_Termino), fill = log10_FDR), 
             shape = 21, alpha = 0.8, color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "-log10(FDR)\n(Significancia)") +
  scale_size_continuous(range = c(2, 8), name = "Genes en\nel término") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  labs(title = "Convergencia Funcional entre Clusters Evolutivos (MST-kNN)",
       subtitle = paste("Términos de", FUENTE_ANALIZAR, "presentes en múltiples clusters"),
       x = "ID del Cluster", y = "Proceso Biológico (GO:BP)")

ruta_plot <- file.path(DIR_RESULTS, paste0("Grafico_Burbujas_Cruzado_mstknn_", NOMBRE_BDD, ".png"))
ggsave(ruta_plot, plot = p_cruzado, width = 12, height = 8, dpi = 300)
cat("Gráfico guardado en:", ruta_plot, "\n")
