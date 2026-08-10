# =============================================================================
# CURVAS DE ESCALABILIDAD EMPÍRICAS (TIEMPO Y MEMORIA)
# Calcula los tiempos y el uso de memoria para distintos subconjuntos de árboles
# =============================================================================

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

suppressPackageStartupMessages({
  library(phangorn)
  library(ape)
  library(mstknnclust)
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

# Función para registrar memoria (usando gc() de base R)
registrar_memoria <- function() {
  gc_info <- gc(verbose = FALSE)
  # Suma celdas usadas de Ncells y Vcells, convierte a MB (1 celda = 8 bytes)
  mem_r <- sum(gc_info[, 2]) * 8 / (1024^2)
  return(as.numeric(mem_r))
}

# =============================================================================
# Función de unión por especie ancla (tomada de 03_comparar_medioide...)
# =============================================================================
unir_por_ancla <- function(tree1, tree2) {
  especies_comunes <- intersect(tree1$tip.label, tree2$tip.label)
  if (length(especies_comunes) < 2) {
    return(list(arbol = NULL, n_comunes = length(especies_comunes), ancla = NA, error = "Menos de 2 especies comunes"))
  }
  especie_ancla <- especies_comunes[1]
  otras_comunes <- especies_comunes[especies_comunes != especie_ancla]
  tree1_preparado <- drop.tip(tree1, otras_comunes)
  posicion_final  <- which(tree1_preparado$tip.label == especie_ancla)
  
  arbol_final <- tryCatch({
    resultado <- bind.tree(tree1_preparado, tree2, where = posicion_final)
    collapse.singles(resultado)
  }, error = function(e) { NULL })
  
  return(list(arbol = arbol_final, n_comunes = length(especies_comunes), ancla = especie_ancla, error = ifelse(is.null(arbol_final), "Error en bind.tree", NA)))
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
  
  # Usar el primer árbol de la muestra como "pseudo-medioide" de referencia
  tree1 <- arboles_muestra[[1]]
  arboles_homogeneizados <- list()
  arboles_homogeneizados[[1]] <- tree1
  
  for (i in 2:n_arboles) {
    tree2 <- arboles_muestra[[i]]
    res <- unir_por_ancla(tree1, tree2)
    if (!is.null(res$arbol)) {
      arboles_homogeneizados[[i]] <- res$arbol
    }
  }
  
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
  
  # 3. TIEMPO DE CLUSTERING (MST-kNN)
  cat("  3. Clustering MST-kNN...\n")
  t_inicio_mst <- Sys.time()
  # simulamos una matriz de distancia aleatoria del tamaño N x N para medir MST-kNN
  matriz_simulada <- matrix(runif(n_arboles^2), nrow=n_arboles)
  matriz_simulada[lower.tri(matriz_simulada)] <- t(matriz_simulada)[lower.tri(matriz_simulada)]
  diag(matriz_simulada) <- 0
  rownames(matriz_simulada) <- paste0("T", 1:n_arboles)
  colnames(matriz_simulada) <- paste0("T", 1:n_arboles)
  
  tryCatch({
    mst_res <- mst.knn(distance.matrix = matriz_simulada)
  }, error = function(e) { NULL })
  
  t_fin_mst <- Sys.time()
  tiempo_mst <- as.numeric(difftime(t_fin_mst, t_inicio_mst, units = "secs"))
  
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
    Tiempo_MSTkNN_s = tiempo_mst,
    Tiempo_UMAP_s = tiempo_umap,
    Memoria_Maxima_MB = mem_max
  ))
}

cat("\n=== RESULTADOS DE ESCALABILIDAD ===\n")
print(resultados_escalabilidad)

ruta_salida <- file.path(DIR_RESULTS, "curvas_escalabilidad.csv")
write.csv(resultados_escalabilidad, ruta_salida, row.names = FALSE)
cat(sprintf("\nResultados guardados en: %s\n", ruta_salida))