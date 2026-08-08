# =============================================================================
# VALIDACIÓN NULA DE ENRIQUECIMIENTO FUNCIONAL (GO) MEDIANTE PERMUTACIONES
# Genera agrupamientos aleatorios conservando el tamaño de los clústeres
# =============================================================================

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

library(openxlsx)

# Instalar/Cargar paquetes de enriquecimiento si no están
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(doParallel)
  library(foreach)
})

# =============================================================================
# PARÁMETROS
# =============================================================================
N_PERMUTACIONES <- 50
SEMILLA         <- 42
PVAL_CUTOFF     <- 0.05
QVAL_CUTOFF     <- 0.05
UMBRAL_FDR      <- 0.05 # Para comparar significancia ajustada

cat("=== INICIANDO VALIDACIÓN NULA DE ENRIQUECIMIENTO GO ===\n")
cat(sprintf("Permutaciones: %d | Semilla: %d\n", N_PERMUTACIONES, SEMILLA))

# =============================================================================
# 1. CARGAR DATOS REALES
# =============================================================================
ruta_genes <- file.path(DIR_CACHE, "mstknn_iter_genes_hgnc.rds")
if (!file.exists(ruta_genes)) {
  stop("No se encontró el rds de genes: ", ruta_genes)
}

# Leer genes y clusters asignados
df_genes <- readRDS(ruta_genes)
df_validos <- df_genes[!is.na(df_genes$HGNC_Symbol) & !is.na(df_genes$Cluster_Final), ]
genes_universe <- unique(df_validos$HGNC_Symbol)

cat(sprintf("Total de genes con símbolo HGNC para validación: %d\n", length(genes_universe)))

# Guardar tamaños de clústeres reales
tamanos_reales <- table(df_validos$Cluster_Final)

# =============================================================================
# 2. DEFINIR FUNCIÓN DE ENRIQUECIMIENTO
# =============================================================================
calcular_enriquecimiento <- function(df_datos, universo) {
  # Filtrar clusters pequeños para acelerar
  tamanos <- table(df_datos$Cluster_Final)
  clusters_validos <- names(tamanos)[tamanos >= 10]
  df_filtrado <- df_datos[df_datos$Cluster_Final %in% clusters_validos, ]
  
  if (nrow(df_filtrado) == 0) return(0)
  
  tryCatch({
    res <- compareCluster(
      geneCluster   = HGNC_Symbol ~ Cluster_Final,
      data          = df_filtrado,
      fun           = "enrichGO",
      universe      = universo,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = PVAL_CUTOFF,
      qvalueCutoff  = QVAL_CUTOFF,
      readable      = FALSE
    )
    
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      df_res <- as.data.frame(res)
      return(sum(df_res$p.adjust < UMBRAL_FDR))
    }
    return(0)
  }, error = function(e) {
    cat("  [!] Error en compareCluster:", e$message, "\n")
    return(0)
  })
}

# =============================================================================
# 3. ENRIQUECIMIENTO DEL MODELO REAL
# =============================================================================
cat("Calculando enriquecimiento para los clústeres reales...\n")
terminos_reales <- calcular_enriquecimiento(df_validos, genes_universe)
cat(sprintf("Términos GO significativos (FDR < %.2f) en modelo real: %d\n", UMBRAL_FDR, terminos_reales))

# =============================================================================
# 4. ITERACIONES DEL MODELO NULO (EN PARALELO)
# =============================================================================
cat(sprintf("\nCalculando enriquecimiento para %d permutaciones en paralelo...\n", N_PERMUTACIONES))

# Configurar cluster paralelo (usando todos los núcleos menos 1)
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Exportar variables y paquetes necesarios a los workers
terminos_nulos <- foreach(i = 1:N_PERMUTACIONES, .combine = c, .packages = c("clusterProfiler", "org.Hs.eg.db")) %dopar% {
  set.seed(SEMILLA + i)
  # Crear copia y permutar aleatoriamente las etiquetas de los clústeres
  df_permutado <- df_validos
  df_permutado$Cluster_Final <- sample(df_permutado$Cluster_Final)
  
  # Calcular enriquecimiento
  res_iter <- calcular_enriquecimiento(df_permutado, genes_universe)
  res_iter
}

stopCluster(cl)

cat("Permutaciones finalizadas.\n")

# =============================================================================
# 5. RESULTADOS Y EXPORTACIÓN
# =============================================================================
p_valor_empirico <- sum(terminos_nulos >= terminos_reales) / N_PERMUTACIONES

cat("\n=== RESULTADOS FINALES ===\n")
cat(sprintf("Términos reales      : %d\n", terminos_reales))
cat(sprintf("Media términos nulos : %.2f (SD: %.2f)\n", mean(terminos_nulos), sd(terminos_nulos)))
cat(sprintf("Valor p empírico     : %.4f\n", p_valor_empirico))

resultados_df <- data.frame(
  Modelo = c("Real", paste0("Nulo_", 1:N_PERMUTACIONES)),
  Terminos_Significativos = c(terminos_reales, terminos_nulos)
)

ruta_salida <- file.path(DIR_RESULTS, "validacion_nula_GO.csv")
write.csv(resultados_df, ruta_salida, row.names = FALSE)
cat(sprintf("\nResultados guardados en: %s\n", ruta_salida))

# =============================================================================
# 6. GRÁFICO HISTOGRAMA DE DISTRIBUCIÓN NULA
# =============================================================================
library(ggplot2)
p_hist <- ggplot(resultados_df[-1, ], aes(x = Terminos_Significativos)) +
  geom_histogram(binwidth = max(1, round(max(terminos_nulos)/20)), 
                 fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = terminos_reales), color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = terminos_reales, y = Inf, label = paste("Valor Real:", terminos_reales), 
           vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
  labs(title = "Distribución Nula de Términos GO Significativos",
       subtitle = sprintf("N=%d Permutaciones | p-value empírico = %.4f", N_PERMUTACIONES, p_valor_empirico),
       x = "Términos GO Significativos",
       y = "Frecuencia") +
  theme_minimal(base_size = 14)

ruta_plot <- file.path(DIR_RESULTS, "histograma_modelo_nulo_GO.png")
ggsave(ruta_plot, plot = p_hist, width = 8, height = 6, dpi = 300)
cat(sprintf("Histograma guardado en: %s\n", ruta_plot))
