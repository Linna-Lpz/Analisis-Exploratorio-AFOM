# =============================================================================
# VALIDACIÓN NULA DE ENRIQUECIMIENTO FUNCIONAL (GO) MEDIANTE PERMUTACIONES
# Genera agrupamientos aleatorios conservando el tamaño de los clústeres
# =============================================================================

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

library(openxlsx)

# Instalar/Cargar paquetes de enriquecimiento si no están
suppressPackageStartupMessages({
  library(gprofiler2)
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
ruta_pca_excel <- file.path(DIR_RESULTS, paste0("pca_coordenadas_", NOMBRE_BDD, ".xlsx"))

if (!file.exists(ruta_pca_excel)) {
  stop("No se encontró el Excel PCA: ", ruta_pca_excel)
}

hojas <- getSheetNames(ruta_pca_excel)
if (!"Genes_Agrupados" %in% hojas) {
  stop("No se encontró la hoja 'Genes_Agrupados' en el Excel PCA.")
}

df_genes <- read.xlsx(ruta_pca_excel, sheet = "Genes_Agrupados", startRow = 2, colNames = TRUE)

# Determinar columna de cluster
algoritmo <- if(exists("ALGORITMO_DOWNSTREAM")) ALGORITMO_DOWNSTREAM else "AUTO"
COL_CLUSTER <- if (algoritmo == "MST-kNN") {
  "Cluster_MSTKNN_Iter"
} else if (algoritmo == "K-Means") {
  "Cluster_KMeans_Iter"
} else if (algoritmo == "CLARA") {
  "Cluster_CLARA_Iter"
} else if (algoritmo == "PAM") {
  "Cluster_PAM"
} else {
  "Cluster_MSTKNN_Iter"
}

if (!COL_CLUSTER %in% colnames(df_genes)) {
  cols_disponibles <- grep("^Cluster_", colnames(df_genes), value = TRUE)
  if (length(cols_disponibles) > 0) {
    COL_CLUSTER <- cols_disponibles[1]
  } else {
    stop("No hay columnas de cluster disponibles en Genes_Agrupados")
  }
}

df_genes$Cluster_Final <- df_genes[[COL_CLUSTER]]
df_validos <- df_genes[!is.na(df_genes$HGNC_Symbol) & !is.na(df_genes$Cluster_Final), ]
genes_universe <- unique(df_validos$HGNC_Symbol)

cat(sprintf("Total de genes con símbolo HGNC para validación: %d\n", length(genes_universe)))

# Guardar tamaños de clústeres reales
tamanos_reales <- table(df_validos$Cluster_Final)

# =============================================================================
# 2. DEFINIR FUNCIÓN DE ENRIQUECIMIENTO (USANDO gprofiler2)
# =============================================================================
calcular_enriquecimiento <- function(df_datos, universo) {
  # Filtrar clusters pequeños para acelerar (usar n >= 15 igual que el script original de mst-knn si se desea, o 10)
  tamanos <- table(df_datos$Cluster_Final)
  clusters_validos <- names(tamanos)[tamanos >= 10]
  
  df_filtrado <- df_datos[df_datos$Cluster_Final %in% clusters_validos & !is.na(df_datos$HGNC_Symbol), ]
  
  if (nrow(df_filtrado) == 0) return(0)
  
  # Preparar lista de queries (un vector de genes por cada clúster)
  query_list <- split(df_filtrado$HGNC_Symbol, df_filtrado$Cluster_Final)
  
  tryCatch({
    resultado <- gost(
      query              = query_list,
      organism           = "hsapiens",
      sources            = "GO:BP", # Validamos ontología biológica principal
      correction_method  = "fdr",
      user_threshold     = UMBRAL_FDR,
      significant        = TRUE,
      measure_underrepresentation = FALSE,
      evcodes            = FALSE,
      custom_bg          = universo
    )
    
    if (!is.null(resultado) && !is.null(resultado$result) && nrow(resultado$result) > 0) {
      return(nrow(resultado$result))
    }
    return(0)
  }, error = function(e) {
    cat("  [!] Error en gost:", e$message, "\n")
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
terminos_nulos <- foreach(i = 1:N_PERMUTACIONES, .combine = c, .packages = c("gprofiler2")) %dopar% {
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
