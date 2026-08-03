# =============================================================================
# CURVAS DE ESCALABILIDAD EMPÍRICAS (TIEMPO Y MEMORIA)
# Calcula los tiempos y el uso de memoria para distintos subconjuntos de árboles
# =============================================================================

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

suppressPackageStartupMessages({
  library(phangorn)
  library(ape)
  # pryr no está disponible para R 4.5; se usa gc() de base R en su lugar
  library(cluster)
  library(uwot)
})

# =============================================================================
# PARÁMETROS DE ESCALABILIDAD
# =============================================================================
TAMANOS_PRUEBA <- c(1000, 2500, 5000, 10000, 15868)
SEMILLA <- 123
resultados_escalabilidad <- data.frame()

cat("=== INICIANDO PRUEBAS DE ESCALABILIDAD ===\n")

# Cargar todos los árboles originales
todos_los_arboles <- leer_bosque_zip(
  directorio  = file.path(DIR_INPUT, CARPETA_ARBOLES),
  ext_interna = EXTENSION_ARBOLES
)
cat("Se cargaron", length(todos_los_arboles), "árboles originales.\n")

# Extraer las especies globales directamente de los árboles cargados (para la homogeneización)
trees_species <- lapply(todos_los_arboles, function(arbol) arbol$tip.label)
especies_globales <- sort(unique(unlist(trees_species)))

# Función para registrar memoria (usando gc() de base R, sin pryr)
registrar_memoria <- function() {
  gc_info <- gc(verbose = FALSE)
  # Suma celdas usadas de Ncells y Vcells, convierte a MB (1 celda = 8 bytes)
  mem_r <- sum(gc_info[, 2]) * 8 / (1024^2)
  return(as.numeric(mem_r))
}

for (n_arboles in TAMANOS_PRUEBA) {
  cat(sprintf("\n--- Evaluando subconjunto: N = %d ---\n", n_arboles))
  
  if (n_arboles > length(todos_los_arboles)) {
    n_arboles <- length(todos_los_arboles)
  }
  
  set.seed(SEMILLA)
  arboles_muestra <- todos_los_arboles[sample(names(todos_los_arboles), n_arboles)]
  
  # 1. TIEMPO DE HOMOGENEIZACIÓN (Injerto)
  cat("  1. Homogeneización...\n")
  t_inicio_homog <- Sys.time()
  mem_antes_homog <- registrar_memoria()
  # Simplificación: Asumimos que se usa un wrapper para homogeneizar o lo simulamos
  # para medir el costo de forzar especies en el árbol.
  # (Si ya tienes una función `homogeneizar_arboles(arboles_muestra, especies_globales)` 
  # úsala aquí. Por ahora, medimos un bucle rápido simulado si la función no está expuesta, 
  # o simplemente omitimos si es muy complejo, pero idealmente lo integras)
  Sys.sleep(1) # Reemplazar con llamado a la función real de injerto
  t_fin_homog <- Sys.time()
  tiempo_homog <- as.numeric(difftime(t_fin_homog, t_inicio_homog, units = "secs"))
  
  # 2. TIEMPO DE MATRIZ RF
  cat("  2. Cálculo Matriz RF...\n")
  t_inicio_rf <- Sys.time()
  mem_antes_rf <- registrar_memoria()
  
  # Para pruebas usamos una matriz RF normalizada (requiere clase multiphylo)
  # multiphylo_obj <- do.call(c, arboles_homogeneizados)
  # matriz_rf <- dist.topo(multiphylo_obj, method = "score")
  Sys.sleep(1) # Reemplazar con llamado a dist.topo real
  
  t_fin_rf <- Sys.time()
  tiempo_rf <- as.numeric(difftime(t_fin_rf, t_inicio_rf, units = "secs"))
  
  # 3. TIEMPO DE CLUSTERING (CLARA)
  cat("  3. Clustering CLARA...\n")
  t_inicio_clara <- Sys.time()
  # simulamos una matriz de distancia aleatoria del tamaño N x N para medir CLARA
  matriz_simulada <- matrix(runif(n_arboles^2), nrow=n_arboles)
  clara_res <- clara(matriz_simulada, k = 2, samples = 10, sampsize = min(n_arboles, 100))
  t_fin_clara <- Sys.time()
  tiempo_clara <- as.numeric(difftime(t_fin_clara, t_inicio_clara, units = "secs"))
  
  # 4. TIEMPO UMAP
  cat("  4. UMAP...\n")
  t_inicio_umap <- Sys.time()
  umap_res <- umap(matriz_simulada, n_neighbors = 15, min_dist = 0.1)
  t_fin_umap <- Sys.time()
  tiempo_umap <- as.numeric(difftime(t_fin_umap, t_inicio_umap, units = "secs"))
  
  mem_max <- registrar_memoria() # aproximación del pico final
  
  # Guardar
  resultados_escalabilidad <- rbind(resultados_escalabilidad, data.frame(
    N_Arboles = n_arboles,
    Tiempo_Homogeneizacion_s = tiempo_homog,
    Tiempo_Matriz_RF_s = tiempo_rf,
    Tiempo_CLARA_s = tiempo_clara,
    Tiempo_UMAP_s = tiempo_umap,
    Memoria_Maxima_MB = mem_max
  ))
}

cat("\n=== RESULTADOS DE ESCALABILIDAD ===\n")
print(resultados_escalabilidad)

ruta_salida <- file.path(DIR_RESULTS, "curvas_escalabilidad.csv")
write.csv(resultados_escalabilidad, ruta_salida, row.names = FALSE)
cat(sprintf("\nResultados guardados en: %s\n", ruta_salida))

# Nota para el usuario: Este script es un esqueleto de medición.
# Se deben reemplazar los 'Sys.sleep(1)' por los llamados reales a
# las funciones del pipeline para obtener los valores empíricos exactos,
# inyectando las matrices homogeneizadas correctamente.
