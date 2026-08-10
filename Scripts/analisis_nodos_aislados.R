# =============================================================================
# analisis_nodos_aislados.R
# PROPÓSITO: Extraer y caracterizar la distribución de cobertura original
#            de los 67 nodos aislados en el grafo MST-kNN
# =============================================================================

library(openxlsx)
library(here)

source(here::here("Scripts", "config.R"))

cat("=== ANÁLISIS DE NODOS AISLADOS MST-kNN ===\n\n")

# 1. Leer la lista de TODOS los árboles desde la matriz RF (15.868)
cat("--- [1] Cargando nombres de todos los árboles (15.868) desde matriz RF ---\n")
ruta_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")
if (!file.exists(ruta_matriz)) stop("No se encontró matriz_rf.rds")

# Solo necesitamos los rownames, no la matriz completa en memoria
# Cargar en un entorno temporal
e <- new.env()
e$mat <- readRDS(ruta_matriz)
todos_arboles <- rownames(e$mat)
n_total <- length(todos_arboles)
rm(e); gc(verbose = FALSE)
cat(sprintf("  Total de árboles en corpus (matriz RF): %d\n", n_total))

# 2. Leer los árboles CON asignación MST-kNN (15.801)
cat("\n--- [2] Cargando asignaciones MST-kNN iterativas ---\n")
ruta_asig <- file.path(DIR_CACHE, "mstknn_iter_asignaciones.rds")
if (!file.exists(ruta_asig)) stop("No se encontró mstknn_iter_asignaciones.rds")

asig_df <- readRDS(ruta_asig)
arboles_asignados <- asig_df$Arbol
n_asignados <- length(arboles_asignados)
cat(sprintf("  Árboles con asignación (conectados): %d\n", n_asignados))

# 3. Identificar los nodos aislados (setdiff)
nodos_aislados <- setdiff(todos_arboles, arboles_asignados)
n_aislados <- length(nodos_aislados)
cat(sprintf("  Nodos aislados (sin asignación): %d\n", n_aislados))
cat(sprintf("  Verificación: %d + %d = %d (esperado: %d)\n",
            n_asignados, n_aislados, n_asignados + n_aislados, n_total))

# 4. Cargar cobertura original (bosque sin injertar = archivo más pequeño)
cat("\n--- [3] Cargando cobertura original de los árboles ---\n")
archivos_bosque <- list.files(DIR_CACHE, pattern = "^bosque_.*\\.rds$", full.names = TRUE)
if (length(archivos_bosque) == 0) stop("No hay bosque en cache.")

bosque_rds <- archivos_bosque[which.min(file.size(archivos_bosque))]
cat(sprintf("  Usando: %s (%.1f MB)\n",
            basename(bosque_rds), file.size(bosque_rds) / 1e6))

bosque <- readRDS(bosque_rds)
n_hojas_vec <- sapply(bosque, function(t) length(t$tip.label))
rm(bosque); gc(verbose = FALSE)

cat(sprintf("  Árboles en bosque: %d\n", length(n_hojas_vec)))

# Crear data.frame de cobertura
cobertura_df <- data.frame(
  Arbol     = names(n_hojas_vec),
  Cobertura = as.integer(n_hojas_vec),
  stringsAsFactors = FALSE
)

# 5. Distribución de cobertura GLOBAL (todos los 15.868)
cat("\n--- [4] Distribución de cobertura GLOBAL (15.868 árboles) ---\n")
cob_global <- cobertura_df$Cobertura
cat(sprintf("  Media   : %.2f\n", mean(cob_global)))
cat(sprintf("  Mediana : %.0f\n", median(cob_global)))
cat(sprintf("  Mínimo  : %d\n", min(cob_global)))
cat(sprintf("  Máximo  : %d\n", max(cob_global)))
cat(sprintf("  SD      : %.2f\n", sd(cob_global)))

# 6. Distribución de cobertura de los NODOS AISLADOS (67)
cat("\n--- [5] Distribución de cobertura NODOS AISLADOS (67 árboles) ---\n")
cob_aislados <- cobertura_df$Cobertura[cobertura_df$Arbol %in% nodos_aislados]
n_encontrados_aislados <- length(cob_aislados)

if (n_encontrados_aislados == 0) {
  cat("  ERROR: No se encontraron coberturas para los nodos aislados.\n")
  cat("  IDs de nodos aislados (primeros 5):", paste(head(nodos_aislados, 5), collapse=", "), "\n")
  cat("  IDs en bosque (primeros 5):", paste(head(cobertura_df$Arbol, 5), collapse=", "), "\n")
} else {
  cat(sprintf("  Nodos aislados encontrados: %d de %d\n", n_encontrados_aislados, n_aislados))
  cat(sprintf("  Media   : %.2f\n", mean(cob_aislados)))
  cat(sprintf("  Mediana : %.0f\n", median(cob_aislados)))
  cat(sprintf("  Mínimo  : %d\n", min(cob_aislados)))
  cat(sprintf("  Máximo  : %d\n", max(cob_aislados)))
  cat(sprintf("  SD      : %.2f\n", sd(cob_aislados)))
  
  # Distribución por rangos de cobertura (cuántos tienen < 50, 50-100, 100-150, 150-190)
  cat("\n  Distribución por rango de cobertura original:\n")
  rangos <- cut(cob_aislados,
                breaks = c(0, 50, 100, 150, 180, 190),
                labels = c("0-50", "51-100", "101-150", "151-180", "181-190"),
                include.lowest = TRUE)
  print(table(rangos))
  
  # Porcentaje de injerto promedio (hojas injertadas / 190)
  pct_injerto_aislados <- (190 - cob_aislados) / 190 * 100
  cat(sprintf("\n  Porcentaje de injerto promedio (nodos aislados): %.2f%%\n", mean(pct_injerto_aislados)))
  cat(sprintf("  vs. porcentaje de injerto promedio GLOBAL        : %.2f%%\n", (190 - mean(cob_global)) / 190 * 100))
  
  # Cuántos tienen cobertura baja (< 100 especies, > 47% injertado)
  n_baja_cobertura <- sum(cob_aislados < 100)
  n_alta_cobertura <- sum(cob_aislados >= 150)
  cat(sprintf("\n  Nodos aislados con cobertura < 100 esp. (>47%% injertado): %d (%.1f%%)\n",
              n_baja_cobertura, 100 * n_baja_cobertura / n_encontrados_aislados))
  cat(sprintf("  Nodos aislados con cobertura >= 150 esp. (<21%% injertado): %d (%.1f%%)\n",
              n_alta_cobertura, 100 * n_alta_cobertura / n_encontrados_aislados))
  cat(sprintf("  Nodos aislados con cobertura = 190 esp. (sin injerto)     : %d (%.1f%%)\n",
              sum(cob_aislados == 190), 100 * sum(cob_aislados == 190) / n_encontrados_aislados))
}

# 7. Prueba de Wilcoxon: cobertura aislados vs. conectados
cat("\n--- [6] Prueba de Wilcoxon: aislados vs. conectados ---\n")
cob_conectados <- cobertura_df$Cobertura[cobertura_df$Arbol %in% arboles_asignados]
cat(sprintf("  Media cobertura conectados: %.2f\n", mean(cob_conectados, na.rm=TRUE)))
cat(sprintf("  Media cobertura aislados  : %.2f\n", mean(cob_aislados,   na.rm=TRUE)))

if (length(cob_aislados) >= 2 && length(cob_conectados) >= 2) {
  wt <- wilcox.test(cob_aislados, cob_conectados)
  cat(sprintf("  Prueba Wilcoxon p-value: %.4f\n", wt$p.value))
  if (wt$p.value < 0.05) {
    cat("  → Diferencia SIGNIFICATIVA en cobertura entre aislados y conectados.\n")
    if (mean(cob_aislados) < mean(cob_conectados)) {
      cat("  → Los nodos aislados tienen MENOR cobertura original (más injertados).\n")
      cat("  → IMPLICACIÓN: el aislamiento puede ser parcialmente un artefacto del injerto.\n")
    } else {
      cat("  → Los nodos aislados tienen MAYOR o similar cobertura original.\n")
      cat("  → IMPLICACIÓN: son genuinamente divergentes topológicamente.\n")
    }
  } else {
    cat("  → No hay diferencia significativa. Ambos grupos tienen cobertura similar.\n")
    cat("  → IMPLICACIÓN: la exclusión de los 67 no introduce sesgo de cobertura.\n")
  }
}

# 8. Guardar tabla de nodos aislados con su cobertura
cat("\n--- [7] Guardando tabla de nodos aislados ---\n")
aislados_con_cobertura <- cobertura_df[cobertura_df$Arbol %in% nodos_aislados, ]
aislados_con_cobertura$Pct_Injerto <- round((190 - aislados_con_cobertura$Cobertura) / 190 * 100, 2)
aislados_con_cobertura <- aislados_con_cobertura[order(aislados_con_cobertura$Cobertura), ]

ruta_csv_out <- file.path(DIR_RESULTS, "nodos_aislados_cobertura.csv")
write.csv(aislados_con_cobertura, ruta_csv_out, row.names = FALSE)
cat(sprintf("  CSV guardado: %s\n", ruta_csv_out))

# 9. Resumen final
cat("\n=== RESUMEN PARA EL INFORME ===\n")
cat(sprintf("Corpus total OrthoMaM:          %d árboles\n", n_total))
cat(sprintf("Nodos en grafo MST-kNN (k=5):   %d árboles\n", n_asignados))
cat(sprintf("Nodos aislados (excluidos):      %d árboles (%.2f%% del corpus)\n",
            n_aislados, 100 * n_aislados / n_total))
if (n_encontrados_aislados > 0) {
  cat(sprintf("\nCobertura original nodos aislados:\n"))
  cat(sprintf("  Media = %.1f esp., Mediana = %.0f esp., Mín = %d, Máx = %d\n",
              mean(cob_aislados), median(cob_aislados), min(cob_aislados), max(cob_aislados)))
  cat(sprintf("  Porcentaje promedio de hojas injertadas: %.1f%%\n", mean(pct_injerto_aislados)))
  cat(sprintf("\nCobertura original corpus general:\n"))
  cat(sprintf("  Media = %.1f esp., Mediana = %.0f esp.\n", mean(cob_global), median(cob_global)))
}
cat("\n=== FIN DEL ANÁLISIS ===\n")
