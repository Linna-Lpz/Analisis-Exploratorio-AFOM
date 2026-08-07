# =============================================================================
# 16_ANALISIS_SENSIBILIDAD.R
# Análisis de sensibilidad del medioide (Respuesta a Revisor 2.8)
#
# VERSIÓN ACTUALIZADA: Utiliza MST-kNN + subdivisión iterativa (el algoritmo
# principal del estudio) para garantizar que la partición comparada no esté
# penalizada por inestabilidad de algoritmos euclidianos (CLARA/K-Means).
# =============================================================================

library(ape)
library(TreeDist)
library(openxlsx)
library(cluster)
library(mstknnclust)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# FUNCION ARI
# =============================================================================
adjusted_rand_index <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c) / (a + b + c + d)) /
    ((a + b + a + c) / 2 - (a + b) * (a + c) / (a + b + c + d))
  return(ARI)
}

# =============================================================================
# 1. IDENTIFICAR SEGUNDO MEDIOIDE Y PARTICIÓN ORIGINAL (MST-kNN)
# =============================================================================
cat("=== LEYENDO SEGUNDO MEDIOIDE Y PARTICION ORIGINAL ===\n")

ruta_ranking <- file.path(DIR_RESULTS, paste0("ranking_medioide_", NOMBRE_BDD, ".xlsx"))
ranking <- read.xlsx(ruta_ranking, sheet = "Ranking_Medioide", startRow = 2)
segundo_medioide <- ranking$nombre_arbol[ranking$posicion == 2]
cat("Segundo medioide identificado:", segundo_medioide, "\n")

ruta_mstknn_final <- file.path(DIR_RESULTS, "mstknn_subdivision_iterativa.xlsx")
if(!file.exists(ruta_mstknn_final)){
  stop("No se encontró la partición final en ", ruta_mstknn_final)
}

asig_original <- read.xlsx(ruta_mstknn_final, sheet = "Asignaciones_Finales", startRow = 2)
k_original <- length(unique(asig_original$Cluster_Final))
cat("Partición original (MST-kNN) cargada. K final =", k_original, "clústeres.\n")

# =============================================================================
# 2. INJERTO CON EL SEGUNDO MEDIOIDE
# =============================================================================
cat("\n=== INJERTANDO CON SEGUNDO MEDIOIDE ===\n")
arboles <- leer_bosque_zip(directorio = file.path(DIR_INPUT, CARPETA_ARBOLES),
                           ext_interna = EXTENSION_ARBOLES,
                           dir_cache   = DIR_CACHE) 

tree1_sens <- arboles[[segundo_medioide]]
nombres_arboles <- names(arboles)
n_arboles <- length(arboles)
arboles_unidos_sens <- vector("list", n_arboles)
names(arboles_unidos_sens) <- nombres_arboles

unir_por_ancla <- function(tree1, tree2) {
  especies_comunes <- intersect(tree1$tip.label, tree2$tip.label)
  if (length(especies_comunes) < 2) return(NULL)
  especie_ancla <- especies_comunes[1]
  otras_comunes <- especies_comunes[especies_comunes != especie_ancla]
  tree1_preparado <- drop.tip(tree1, otras_comunes)
  posicion_final  <- which(tree1_preparado$tip.label == especie_ancla)
  arbol_final <- tryCatch({
    resultado <- bind.tree(tree1_preparado, tree2, where = posicion_final)
    collapse.singles(resultado)
  }, error = function(e) NULL)
  return(arbol_final)
}

for (i in seq_len(n_arboles)) {
  tree2 <- arboles[[i]]
  nombre_tree2 <- nombres_arboles[i]
  if (nombre_tree2 == segundo_medioide) {
    arboles_unidos_sens[[i]] <- tree1_sens
    next
  }
  res <- unir_por_ancla(tree1_sens, tree2)
  if (!is.null(res)) {
    arboles_unidos_sens[[i]] <- res
  }
}

arboles_unidos_sens <- arboles_unidos_sens[!sapply(arboles_unidos_sens, is.null)]
class(arboles_unidos_sens) <- "multiPhylo"
cat("Árboles injertados exitosamente:", length(arboles_unidos_sens), "\n")

# =============================================================================
# 3. MATRIZ RF DEL NUEVO BOSQUE
# =============================================================================
cat("\n=== CALCULANDO MATRIZ RF (SENSIBILIDAD) ===\n")
matriz_distancias_sens <- RobinsonFoulds(arboles_unidos_sens, normalize = TRUE)
matriz_cuadrada_sens   <- as.matrix(matriz_distancias_sens)
rownames(matriz_cuadrada_sens) <- names(arboles_unidos_sens)
colnames(matriz_cuadrada_sens) <- names(arboles_unidos_sens)
cat("Matriz calculada:", nrow(matriz_cuadrada_sens), "x", ncol(matriz_cuadrada_sens), "\n")

# =============================================================================
# 4. CLUSTERING CON MST-kNN + Iterativo
# =============================================================================
cat("\n=== CLUSTERING MST-kNN INICIAL ===\n")
set.seed(2)
mst_sens_res <- mst.knn(distance.matrix = matriz_cuadrada_sens, suggested.k = 10)

asig_inicial_sens <- as.integer(mst_sens_res$cluster)
names(asig_inicial_sens) <- names(mst_sens_res$cluster)

cat("\n=== SUBDIVISION ITERATIVA (PAM) ===\n")
TAMANO_MAX      <- 200
TAMANO_MIN_DURO <- 15
DELTA_SIL       <- -999
K_SUBDIVISION   <- 2
MAX_ITERACIONES <- 30
SEED            <- 2

subdividir_con_pam <- function(ids_arboles, matriz_completa, k = 2, seed = 2, min_n = 3) {
  n <- length(ids_arboles)
  if (n < max(k + 1, min_n)) {
    return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA, razon = "n_insuficiente"))
  }
  ids_validos <- ids_arboles[ids_arboles %in% rownames(matriz_completa)]
  if (length(ids_validos) < k + 1) {
    return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA, razon = "ids_no_en_matriz"))
  }
  submatriz <- matriz_completa[ids_validos, ids_validos]
  dist_sub  <- as.dist(submatriz)
  set.seed(seed)
  pam_res <- tryCatch(pam(dist_sub, k = k, diss = TRUE), error = function(e) NULL)
  if (is.null(pam_res)) return(list(asig = setNames(rep(1L, n), ids_arboles), sil_avg = NA, razon = "pam_error"))
  sil_obj <- silhouette(pam_res$clustering, dist_sub)
  sil_avg <- mean(sil_obj[, 3])
  list(asig = setNames(as.integer(pam_res$clustering), ids_validos), sil_avg = round(sil_avg, 4), razon = "ok")
}

estado_actual <- asig_inicial_sens
iteracion <- 0L

repeat {
  iteracion <- iteracion + 1L
  if (iteracion > MAX_ITERACIONES) break
  
  tamanos_actuales <- table(estado_actual)
  clusters_grandes <- as.integer(names(tamanos_actuales[tamanos_actuales > TAMANO_MAX & tamanos_actuales >= TAMANO_MIN_DURO]))
  
  if (length(clusters_grandes) == 0) break
  
  proximo_id <- max(estado_actual) + 1L
  nuevo_estado <- estado_actual
  hubo_cambio <- FALSE
  
  for (cid in clusters_grandes) {
    arboles_cid <- names(estado_actual[estado_actual == cid])
    resultado_sub <- subdividir_con_pam(arboles_cid, matriz_cuadrada_sens, k = K_SUBDIVISION, seed = SEED, min_n = TAMANO_MIN_DURO)
    
    if (is.na(resultado_sub$sil_avg) || resultado_sub$razon != "ok") next
    if (resultado_sub$sil_avg < DELTA_SIL) next
    
    asig_nueva <- resultado_sub$asig
    sub_ids <- sort(unique(asig_nueva))
    mapa_ids <- setNames(c(cid, seq(proximo_id, proximo_id + length(sub_ids) - 2L)), sub_ids)
    asig_reetiquetada <- mapa_ids[as.character(asig_nueva)]
    
    nuevo_estado[names(asig_nueva)] <- asig_reetiquetada
    proximo_id <- proximo_id + length(sub_ids) - 1L
    hubo_cambio <- TRUE
  }
  estado_actual <- nuevo_estado
  if (!hubo_cambio) break
}

ids_originales <- sort(unique(estado_actual))
mapa_final <- setNames(seq_along(ids_originales), ids_originales)
cluster_final <- mapa_final[as.character(estado_actual)]

asig_sens_df <- data.frame(
  Arbol = names(estado_actual),
  Cluster_Sensibilidad = as.integer(cluster_final),
  stringsAsFactors = FALSE
)
cat("Nuevos clústeres generados:", length(unique(asig_sens_df$Cluster_Sensibilidad)), "\n")

# =============================================================================
# 5. COMPARACION ARI
# =============================================================================
cat("\n=== CALCULANDO INDICE DE RAND AJUSTADO (ARI) ===\n")
comparacion <- merge(asig_original, asig_sens_df, by="Arbol", all=FALSE)
ari_score <- adjusted_rand_index(comparacion$Cluster_Final, comparacion$Cluster_Sensibilidad)

cat(sprintf("-> ARI Obtenido (MST-kNN): %.4f\n", ari_score))

if (ari_score > 0.7) {
  conclusion <- "Alto (La elección del andamiaje no altera significativamente la topología global. Objeción cerrada)."
} else if (ari_score > 0.4) {
  conclusion <- "Moderado (Existen diferencias, pero se conserva cierta estructura global)."
} else {
  conclusion <- "Bajo (El resultado sigue siendo dependiente del andamiaje, el injerto afecta la variabilidad estructural)."
}
cat("Conclusión ARI:", conclusion, "\n")

resumen_sensibilidad <- data.frame(
  Metrica = c("Segundo Medioide Usado", "Arboles Comparados", "K Original MST-kNN", "ARI Obtenido", "Conclusion"),
  Valor = c(segundo_medioide, nrow(comparacion), k_original, round(ari_score, 4), conclusion)
)

wb <- createWorkbook()
wb <- agregar_hoja_formateada(wb, "Resumen_Sensibilidad", "Análisis de Sensibilidad de Medioide", resumen_sensibilidad)
wb <- agregar_hoja_formateada(wb, "Comparacion_Clusters", "Clusters Originales vs Sensibilidad", comparacion)

ruta_salida_sens <- file.path(DIR_RESULTS, paste0("analisis_sensibilidad_medioide_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_salida_sens, overwrite = TRUE)

cat("\nReporte guardado en:", ruta_salida_sens, "\n")
cat("=== FIN DEL SCRIPT ===\n")
