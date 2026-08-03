# =============================================================================
# 05b_reduccion_dimensional_nMDS.R
# nMDS (no métrico) sobre la matriz RF normalizada — con caché .rds
# Salidas: coordenadas en Excel, gráfico ggplot2, Shepard plot, scree plot, caché .rds
# =============================================================================

library(vegan)
library(ggplot2)
library(patchwork)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# CARGAR MATRIZ RF DESDE CACHÉ (generada en script 04)
# =============================================================================
cat("=== CARGANDO MATRIZ RF ===\n")

ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

if (!file.exists(ruta_cache_matriz)) {
  stop("No se encontró la matriz RF en caché: ", ruta_cache_matriz,
       "\nEjecuta primero el script 04_calcular_matriz.R")
}

matriz_cuadrada <- readRDS(ruta_cache_matriz)
dist_rf         <- as.dist(matriz_cuadrada)

cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# CALCULAR nMDS — con caché condicional
# =============================================================================
cat("\n=== CALCULANDO nMDS ===\n")

ruta_cache_nmds <- file.path(DIR_CACHE, "nmds_coords.rds")

# --- Hiperparámetros ---
N_COMPONENTS <- 3      # dimensiones de salida
TRYMAX       <- 10     # máximo de arranques aleatorios
MAXIT        <- 300    # máximo de iteraciones por arranque
SEED         <- 42

# Referencia de stress (Kruskal 1964)
interpretar_stress <- function(s) {
  dplyr::case_when(
    s < 0.05 ~ "Excelente",
    s < 0.10 ~ "Bueno",
    s < 0.20 ~ "Aceptable",
    TRUE     ~ "Deficiente — considerar más dimensiones"
  )
}

if (file.exists(ruta_cache_nmds)) {
  cat("nMDS encontrado en caché. Cargando...\n")
  tiempo_nmds <- system.time({
    cache_data <- readRDS(ruta_cache_nmds)
    # Compatibilidad con versión anterior que no guardaba salida_consola
    if (is.list(cache_data) && "nmds" %in% names(cache_data)) {
      nmds_resultado <- cache_data$nmds
      salida_consola <- cache_data$log
    } else {
      nmds_resultado <- cache_data
      salida_consola <- character(0)
    }
  })
  cat("Cargado desde caché.\n")
  
} else {
  cat("Calculando nMDS (puede tardar varios minutos con", nrow(matriz_cuadrada), "árboles)...\n")
  cat("trymax =", TRYMAX, "| maxit =", MAXIT, "| k =", N_COMPONENTS, "\n\n")
  
  tiempo_nmds <- system.time({
    set.seed(SEED)
    salida_consola <- capture.output({
      nmds_resultado <- metaMDS(
        comm          = dist_rf,   # matriz de distancias precomputada
        k             = N_COMPONENTS,
        trymax        = TRYMAX,
        maxit         = MAXIT,
        autotransform = FALSE,     # desactivar: no son datos de comunidad
        wascores      = FALSE,     # desactivar: no hay scores de especies
        noshare       = FALSE,     # desactivar: stepacross no aplica con dist precomputada
        trace         = 1          # mostrar progreso en consola (0 = silencioso)
      )
    }, type = "output")
    # Guardar objeto y log en caché
    saveRDS(list(nmds = nmds_resultado, log = salida_consola), file = ruta_cache_nmds)
  })
  cat(paste(salida_consola, collapse = "\n"))
  cat("\nnMDS calculado y guardado en caché.\n")
}

# =============================================================================
# GRÁFICO DE CONVERGENCIA nMDS (Stress por Run)
# =============================================================================
if (length(salida_consola) > 0) {
  cat("\n=== GENERANDO GRÁFICO DE CONVERGENCIA ===\n")
  lineas_run <- grep("^Run [0-9]+ stress ", salida_consola, value = TRUE)
  if (length(lineas_run) > 0) {
    runs <- as.numeric(sub("^Run ([0-9]+) stress .*", "\\1", lineas_run))
    stresses <- as.numeric(sub("^Run [0-9]+ stress ([0-9.]+)", "\\1", lineas_run))
    conv_df <- data.frame(Run = runs, Stress = stresses)
    
    p_conv <- ggplot(conv_df, aes(x = Run, y = Stress)) +
      geom_line(color = "forestgreen", linewidth = 0.8) +
      geom_point(size = 2, color = "forestgreen") +
      labs(
        title = paste0("Convergencia nMDS — Stress por Run (", NOMBRE_BDD, ")"),
        subtitle = paste0("Stress final del mejor modelo: ", round(min(stresses), 4)),
        x = "Run",
        y = "Stress",
        caption = "Un stress más bajo indica una mejor preservación de las distancias originales."
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"),
        plot.caption  = element_text(color = "gray60", size = 8)
      )
    
    ruta_conv <- file.path(DIR_RESULTS, paste0("nmds_convergencia_", NOMBRE_BDD, ".png"))
    ggsave(ruta_conv, plot = p_conv, width = 8, height = 5, dpi = 300)
    cat("Gráfico de convergencia guardado:", ruta_conv, "\n")
  }
}

# =============================================================================
# EVALUAR CALIDAD — STRESS
# =============================================================================
cat("\n=== EVALUACIÓN DE CALIDAD ===\n")

stress_val    <- nmds_resultado$stress
interpretacion <- interpretar_stress(stress_val)

cat(sprintf("Stress final  : %.4f\n", stress_val))
cat(sprintf("Interpretación: %s\n",   interpretacion))
cat(sprintf("Convergencia  : %s\n",
            ifelse(nmds_resultado$converged, "SÍ", "NO — aumentar trymax o maxit")))
cat(sprintf("Tiempo        : %.1f segundos\n", tiempo_nmds["elapsed"]))

# =============================================================================
# CONSTRUIR DATAFRAME DE COORDENADAS
# =============================================================================
coords_nmds <- scores(nmds_resultado, display = "sites")

coords_df <- data.frame(
  Arbol  = rownames(matriz_cuadrada),
  NMDS_1 = coords_nmds[, 1],
  NMDS_2 = coords_nmds[, 2],
  stringsAsFactors = FALSE
)

# --- Liberar matriz de memoria (ya tenemos las coordenadas) ---
rm(matriz_cuadrada, coords_nmds)
gc(verbose = FALSE)

cat("\nPrimeras filas de coordenadas:\n")
print(head(coords_df, 5))

# =============================================================================
# LEER Y UNIR ETIQUETAS DE CLUSTERING (k óptimo + k=10, k=15)
# =============================================================================
resultado    <- unir_etiquetas_clustering(coords_df, DIR_RESULTS)
coords_df    <- resultado$coords_df
k_optimos    <- resultado$k_optimos
k_extra      <- resultado$k_extra

# =============================================================================
# VERIFICAR QUE HAY ETIQUETAS DE CLUSTERING
# =============================================================================
cols_cluster <- grep("^Cluster_", colnames(coords_df), value = TRUE)
if (length(cols_cluster) == 0) {
  stop(
    "[nMDS] No se encontraron asignaciones de clustering en Resultados/.\n",
    "Ejecuta primero los scripts de clustering antes de la visualización:\n",
    "  - clustering_kmeans.R\n",
    "  - clustering_pam.R\n",
    "  - clustering_clara.R\n",
    "Los archivos esperados en Resultados/ son:\n",
    "  kmeans_subdivision_iterativa.xlsx, pam_resultados.xlsx, clara_subdivision_iterativa.xlsx"
  )
}
cat(sprintf("  [OK] Columnas de clustering cargadas (%d): %s\n",
            length(cols_cluster), paste(cols_cluster, collapse = ", ")))

k_km  <- ifelse(is.na(k_optimos["KMeans"]), "?", k_optimos["KMeans"])
k_pam <- ifelse(is.na(k_optimos["PAM"]),    "?", k_optimos["PAM"])
k_cl  <- ifelse(is.na(k_optimos["CLARA"]),  "?", k_optimos["CLARA"])

# =============================================================================
# FUNCIÓN PARA GRAFICAR CLUSTERING EN nMDS
# =============================================================================
graficar_clustering_nmds <- function(df, col_cluster, titulo, subtitulo = "") {
  if (!col_cluster %in% colnames(df)) {
    warning(sprintf("Columna '%s' no encontrada. Saltando gráfico.", col_cluster))
    return(NULL)
  }
  n_clusters <- length(unique(na.omit(df[[col_cluster]])))
  
  # Paleta dinámica: fija hasta 20 clusters, generada dinámicamente para más
  paleta <- if (n_clusters <= 20) {
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
      "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
      "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A")[seq_len(n_clusters)]
  } else {
    hcl.colors(n_clusters, palette = "Dynamic")
  }
  
  ggplot(df, aes(x = NMDS_1, y = NMDS_2, color = .data[[col_cluster]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = paleta,
                       name = "Cluster", na.value = "grey80") +
    labs(title = titulo, subtitle = subtitulo,
         x = "nMDS 1", y = "nMDS 2") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(color = "gray40", size = 9),
      # oculta la leyenda si hay demasiados clusters — evita que domine el panel
      legend.position = if (n_clusters > 15) "none" else "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 5))
}

# =============================================================================
# GRÁFICOS CLUSTERING nMDS — K ÓPTIMO
# =============================================================================
cat("\n=== GENERANDO GRÁFICOS CLUSTERING nMDS (K ÓPTIMO) ===\n")

p_km <- graficar_clustering_nmds(coords_df, "Cluster_KMeans", "K-Means",
                                 sprintf("k = %s  |  Stress = %.4f  |  n = %d", k_km, stress_val, nrow(coords_df)))
p_pa <- graficar_clustering_nmds(coords_df, "Cluster_PAM", "PAM (K-Medoids)",
                                 sprintf("k = %s  |  Stress = %.4f  |  n = %d", k_pam, stress_val, nrow(coords_df)))
p_cl <- graficar_clustering_nmds(coords_df, "Cluster_CLARA", "CLARA",
                                 sprintf("k = %s  |  Stress = %.4f  |  n = %d", k_cl, stress_val, nrow(coords_df)))

# =============================================================================
# GUARDAR GRÁFICOS INDIVIDUALES POR ALGORITMO (K ÓPTIMO)
# =============================================================================
cat("\n=== GUARDANDO GRÁFICOS INDIVIDUALES (K ÓPTIMO) ===\n")

for (g in list(
  list(plot = p_km, nombre = "nmds_kmeans"),
  list(plot = p_pa, nombre = "nmds_pam"),
  list(plot = p_cl, nombre = "nmds_clara")
)) {
  if (!is.null(g$plot)) {
    ruta_g <- file.path(DIR_RESULTS, paste0(g$nombre, "_", NOMBRE_BDD, ".png"))
    ggsave(ruta_g, plot = g$plot, width = 8, height = 6, dpi = 300)
    cat("Gráfico guardado:", ruta_g, "\n")
  }
}
rm(p_km, p_pa, p_cl)
gc(verbose = FALSE)

# =============================================================================
# GRÁFICOS CLUSTERING — K EXTRA (k=10, k=15, etc.)
# =============================================================================
for (ke in k_extra) {
  cat(sprintf("\n=== GRÁFICOS nMDS (K = %d) ===\n", ke))
  
  col_km  <- paste0("Cluster_KMeans_K", ke)
  col_pam <- paste0("Cluster_PAM_K",    ke)
  col_cl  <- paste0("Cluster_CLARA_K",  ke)
  
  pe_km <- graficar_clustering_nmds(coords_df, col_km, "K-Means",
                                    sprintf("k = %d  |  Stress = %.4f  |  n = %d", ke, stress_val, nrow(coords_df)))
  pe_pa <- graficar_clustering_nmds(coords_df, col_pam, "PAM (K-Medoids)",
                                    sprintf("k = %d  |  Stress = %.4f  |  n = %d", ke, stress_val, nrow(coords_df)))
  pe_cl <- graficar_clustering_nmds(coords_df, col_cl, "CLARA",
                                    sprintf("k = %d  |  Stress = %.4f  |  n = %d", ke, stress_val, nrow(coords_df)))
  
  for (g in list(
    list(plot = pe_km, nombre = paste0("nmds_kmeans_K", ke)),
    list(plot = pe_pa, nombre = paste0("nmds_pam_K",    ke)),
    list(plot = pe_cl, nombre = paste0("nmds_clara_K",  ke))
  )) {
    if (!is.null(g$plot)) {
      ruta_g <- file.path(DIR_RESULTS, paste0(g$nombre, "_", NOMBRE_BDD, ".png"))
      ggsave(ruta_g, plot = g$plot, width = 8, height = 6, dpi = 300)
      cat("Gráfico guardado:", ruta_g, "\n")
    }
  }
  rm(pe_km, pe_pa, pe_cl)
}

# =============================================================================
# SHEPARD PLOT — diagnóstico de ajuste
# Compara distancias originales RF vs distancias en espacio nMDS
# Un buen ajuste = puntos cercanos a la línea escalonada
# =============================================================================
shepard_data <- data.frame(
  dist_original = as.vector(dist_rf),
  dist_nmds     = as.vector(dist(scores(nmds_resultado)))
)

nmds_converged <- nmds_resultado$converged
# --- Liberar dist_rf y resultado nMDS (ya extrajimos lo necesario) ---
rm(dist_rf, nmds_resultado)
gc(verbose = FALSE)

p_shepard <- ggplot(shepard_data, aes(x = dist_original, y = dist_nmds)) +
  geom_point(alpha = 0.1, size = 0.8, color = "gray40") +
  geom_smooth(method = "loess", se = FALSE,
              color = "firebrick", linewidth = 0.8) +
  labs(
    title    = paste0("Shepard Plot — nMDS (", NOMBRE_BDD, ")"),
    subtitle = paste0("Stress = ", round(stress_val, 4),
                      "  |  Ajuste ideal = relación monótona"),
    x       = "Distancia RF original",
    y       = "Distancia en espacio nMDS",
    caption = "Línea roja: ajuste LOESS — dispersión baja indica buen ajuste"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    plot.caption  = element_text(color = "gray60", size = 8)
  )

ruta_shepard <- file.path(DIR_RESULTS, paste0("nmds_shepard_", NOMBRE_BDD, ".png"))
ggsave(ruta_shepard, plot = p_shepard, width = 8, height = 6, dpi = 300)
cat("Shepard plot guardado:", ruta_shepard, "\n")

# =============================================================================
# SCREE PLOT nMDS — stress vs. número de dimensiones
# Análogo al scree del PCA: muestra cómo cae el stress al aumentar k.
# Permite elegir el número mínimo de dimensiones con stress aceptable.
# Referencia Kruskal (1964): <0.05 excelente, <0.10 bueno, <0.20 aceptable.
# =============================================================================
cat("\n=== GENERANDO SCREE PLOT nMDS (dimensiones vs. stress) ===\n")

# Recargar dist_rf (fue liberada antes)
scree_dist <- as.dist(readRDS(file.path(DIR_CACHE, "matriz_rf.rds")))

# Dimensiones a evaluar (k=1..6 es suficiente para este análisis)
k_seq <- 1:6

cat(sprintf("Evaluando %d valores de k (dimensiones): %s\n",
            length(k_seq), paste(k_seq, collapse = ", ")))

stress_seq <- numeric(length(k_seq))

for (i in seq_along(k_seq)) {
  kd <- k_seq[i]
  cat(sprintf("  k = %d ...", kd))
  res_tmp <- tryCatch({
    set.seed(SEED)
    metaMDS(
      comm          = scree_dist,
      k             = kd,
      trymax        = 5,      # mínimo para exploración rápida del scree
      maxit         = 100,
      autotransform = FALSE,
      wascores      = FALSE,
      noshare       = FALSE,
      trace         = 0       # silencioso
    )
  }, error = function(e) NULL)
  if (!is.null(res_tmp)) {
    stress_seq[i] <- res_tmp$stress
    cat(sprintf(" stress = %.4f  [%s]\n",
                stress_seq[i], interpretar_stress(stress_seq[i])))
  } else {
    stress_seq[i] <- NA
    cat(" [error — omitido]\n")
  }
}

rm(scree_dist)
gc(verbose = FALSE)

scree_nmds_df <- data.frame(
  k      = k_seq,
  stress = stress_seq,
  stringsAsFactors = FALSE
)
scree_nmds_df$usado <- scree_nmds_df$k == N_COMPONENTS

p_scree_nmds <- ggplot(scree_nmds_df, aes(x = k, y = stress)) +
  # Líneas de referencia Kruskal
  geom_hline(yintercept = c(0.05, 0.10, 0.20),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  annotate("text", x = max(k_seq), y = c(0.052, 0.102, 0.202),
           label = c("Excelente (<0.05)", "Bueno (<0.10)", "Aceptable (<0.20)"),
           hjust = 1, size = 3, color = "gray40") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(aes(shape = usado, size = usado,
                 color = ifelse(usado, "usado", "otro"))) +
  scale_color_manual(values = c("usado" = "firebrick", "otro" = "steelblue"),
                     guide  = "none") +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16),
                     labels = c(`TRUE` = paste0("k = ", N_COMPONENTS,
                                                " (valor usado)"),
                                `FALSE` = "Otros valores"),
                     name   = "") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
  labs(
    title    = paste0("Scree Plot nMDS — Dimensiones vs. Stress (", NOMBRE_BDD, ")"),
    subtitle = paste0(
      "Stress de Kruskal según el número de dimensiones del nMDS.\n",
      "Elegir el k mínimo que alcanza stress aceptable (<0.20)."
    ),
    x       = "Número de dimensiones (k)",
    y       = "Stress de Kruskal",
    caption = paste0("Rombo rojo: k = ", N_COMPONENTS,
                     " (valor usado en el análisis principal)")
  ) +
  scale_x_continuous(breaks = k_seq) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40", size = 9),
    plot.caption  = element_text(color = "gray60", size = 8),
    legend.position = "bottom"
  )

ruta_scree_nmds <- file.path(DIR_RESULTS,
                             paste0("nmds_scree_dimensiones_", NOMBRE_BDD, ".png"))
ggsave(ruta_scree_nmds, plot = p_scree_nmds, width = 9, height = 6, dpi = 300)
cat("Scree plot nMDS guardado:", ruta_scree_nmds, "\n")

# =============================================================================
# EXPORTAR COORDENADAS Y DIAGNÓSTICO A EXCEL
# =============================================================================
cat("\n=== EXPORTANDO RESULTADOS ===\n")

parametros_df <- data.frame(
  Parametro = c("metodo", "engine", "k", "trymax", "maxit", "seed",
                "stress", "interpretacion_stress", "convergencia",
                "n_arboles", "tiempo_calculo_s"),
  Valor     = c("nMDS no métrico", "monoMDS (vegan)", N_COMPONENTS,
                TRYMAX, MAXIT, SEED,
                round(stress_val, 4), interpretacion,
                ifelse(nmds_converged, "Sí", "No"),
                nrow(coords_df), round(tiempo_nmds["elapsed"], 1)),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()

agregar_hoja_formateada(wb, "Coordenadas_nMDS",
                        paste0("Coordenadas nMDS + Etiquetas de Clustering — ", NOMBRE_BDD),
                        coords_df,
                        anchos_col = "auto")

agregar_hoja_formateada(wb, "Diagnostico",
                        "Diagnóstico y Parámetros nMDS",
                        parametros_df,
                        anchos_col = c(30, 35))

ruta_excel <- file.path(DIR_RESULTS, paste0("nmds_coordenadas_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
rm(wb)
gc(verbose = FALSE)
cat("Coordenadas guardadas:", ruta_excel, "\n")

cat("\n=== COMPLETADO ===\n")