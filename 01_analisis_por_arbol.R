# =============================================================================
# 1-ANÁLISIS DE ESPECIES POR ÁRBOL
# Genera: species_analysis_X.csv con ranking, resumen y duplicados
# =============================================================================

library(ape)

# --- Cargar configuración y funciones globales ---
source("config.R")
source(here::here("scripts", "00_funciones_globales.R"))

# --- Lectura ---
bosque_global <- leer_bosque_zip(directorio = DIR_INPUT, ext_interna = EXTENSION_ARBOLES)

# --- Extraer metadata ---
cat("Extrayendo etiquetas de especies...\n")
trees_species <- lapply(bosque_global, function(arbol) arbol$tip.label)

# =============================================================================
# Ranking de especies por frecuencia (orden alfabético)
# =============================================================================
cat("\n=== PASO 2: RANKING DE ESPECIES ===\n")

todas_especies <- sort(unique(unlist(trees_species)))

frecuencias <- sapply(todas_especies, function(sp) {
  sum(sapply(trees_species, function(tips) sp %in% tips))
})

species_counts_df <- data.frame(
  Especie    = todas_especies,
  Presencias = as.integer(frecuencias),
  stringsAsFactors = FALSE
)
species_counts_df <- species_counts_df[order(species_counts_df$Especie), ]

cat("Total especies únicas:", nrow(species_counts_df), "\n")

# =============================================================================
# Número de especies por árbol
# =============================================================================
n_especies <- sapply(trees_species, length)

tree_counts_df <- data.frame(
  Arbol         = names(n_especies),
  N_de_especies = as.integer(n_especies),
  stringsAsFactors = FALSE
)
tree_counts_df <- tree_counts_df[order(tree_counts_df$N_de_especies, decreasing = TRUE), ]
rownames(tree_counts_df) <- NULL

# =============================================================================
# Resumen de frecuencias de conteos
# =============================================================================
resumen_tabla <- as.data.frame(table(tree_counts_df$N_de_especies))
colnames(resumen_tabla) <- c("N_de_especies", "N_de_arboles")
resumen_tabla$N_de_especies <- as.integer(as.character(resumen_tabla$N_de_especies))
resumen_tabla <- resumen_tabla[order(resumen_tabla$N_de_especies, decreasing = TRUE), ]
rownames(resumen_tabla) <- NULL

# =============================================================================
# Árboles duplicados (mismo conjunto de especies)
# =============================================================================
claves <- sapply(trees_species, function(tips) paste(sort(tips), collapse = "|"))

grupos <- split(names(claves), claves)

duplicados_list <- list()
for (clave in names(grupos)) {
  arboles_grupo <- grupos[[clave]]
  if (length(arboles_grupo) > 1) {
    n_sp <- length(strsplit(clave, "\\|")[[1]])
    duplicados_list[[length(duplicados_list) + 1]] <- data.frame(
      Arboles_duplicados = paste(arboles_grupo, collapse = ", "),
      N_de_especies      = as.integer(n_sp),
      N_de_arboles       = as.integer(length(arboles_grupo)),
      stringsAsFactors   = FALSE
    )
  }
}

if (length(duplicados_list) > 0) {
  duplicados_df <- do.call(rbind, duplicados_list)
  duplicados_df <- duplicados_df[order(duplicados_df$N_de_especies, decreasing = TRUE), ]
  rownames(duplicados_df) <- NULL
} else {
  duplicados_df <- data.frame(
    Arboles_duplicados = "NINGUNO",
    N_de_especies      = NA,
    N_de_arboles       = NA,
    stringsAsFactors   = FALSE
  )
}

# =============================================================================
# Exportar a CSV
# =============================================================================
cat("\n=== PASO 6: EXPORTANDO CSV ===\n")

write.csv(species_counts_df,
          file = file.path(DIR_RESULTS, paste0("species_ranking_",    ORIGEN_DATOS, ".csv")),
          row.names = FALSE)

write.csv(resumen_tabla,
          file = file.path(DIR_RESULTS, paste0("species_resumen_",    ORIGEN_DATOS, ".csv")),
          row.names = FALSE)

write.csv(duplicados_df,
          file = file.path(DIR_RESULTS, paste0("species_duplicados_", ORIGEN_DATOS, ".csv")),
          row.names = FALSE)

cat("Archivos guardados:\n")
cat(" -", paste0("species_ranking_",    ORIGEN_DATOS, ".csv"), "\n")
cat(" -", paste0("species_resumen_",    ORIGEN_DATOS, ".csv"), "\n")
cat(" -", paste0("species_duplicados_", ORIGEN_DATOS, ".csv"), "\n")
cat("\n=== COMPLETADO ===\n")