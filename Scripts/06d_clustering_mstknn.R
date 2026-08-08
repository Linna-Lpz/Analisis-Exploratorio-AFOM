# =============================================================================
# CLUSTERING MST-kNN SOBRE MATRIZ RF — versión final
# =============================================================================
library(mstknnclust)
library(igraph)
library(cluster)
library(clValid)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# Operador %||% — debe definirse ANTES del loop que lo usa
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# 1. LEER LA MATRIZ RF
# =============================================================================
cat("Leyendo matriz RF...\n")
ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

if (file.exists(ruta_cache_matriz)) {
  matriz_cuadrada <- readRDS(ruta_cache_matriz)
  cat("Matriz cargada desde caché RDS.\n")
} else {
  cat("Caché RDS no encontrado. Leyendo desde CSV...\n")
  matriz_cuadrada <- as.matrix(
    read.table(file.path(DIR_RESULTS, "matriz_rf_conjunto.csv"),
               sep = ";", header = TRUE, row.names = 1, check.names = FALSE)
  )
  saveRDS(matriz_cuadrada, ruta_cache_matriz)
  cat("Matriz guardada en caché RDS para uso futuro.\n")
}

cat("Matriz cargada:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")

# =============================================================================
# 2. FUNCIÓN AUXILIAR: calcular índices solo con nodos presentes en solución
# =============================================================================
calcular_indices <- function(mst_resultado, matriz_completa) {
  
  k <- mst_resultado$cnumber
  
  if (k <= 1) {
    cat("  AVISO: k=1, índices de calidad no aplicables.\n")
    return(list(Dunn = NA, Connectivity = NA, Silhouette = NA,
                n_nodos = NA, k_efectivo = NA))
  }
  
  nodos_solucion <- names(mst_resultado$cluster)
  n_total        <- nrow(matriz_completa)
  n_solucion     <- length(nodos_solucion)
  
  if (n_solucion < n_total) {
    cat(sprintf("  AVISO: MST-kNN resolvió %d de %d nodos. Subconjuntando matriz.\n",
                n_solucion, n_total))
  }
  
  nodos_validos <- nodos_solucion[nodos_solucion %in% rownames(matriz_completa)]
  
  if (length(nodos_validos) < n_solucion) {
    cat(sprintf("  AVISO: %d nodos de la solución no están en la matriz. Se omiten.\n",
                n_solucion - length(nodos_validos)))
  }
  
  submatriz  <- matriz_completa[nodos_validos, nodos_validos]
  asignacion <- as.integer(mst_resultado$cluster[nodos_validos])
  k_efectivo <- length(unique(asignacion))
  
  if (k_efectivo < 2) {
    cat("  AVISO: k efectivo < 2 tras filtrar nodos. Índices no calculables.\n")
    return(list(Dunn = NA, Connectivity = NA, Silhouette = NA,
                n_nodos = length(nodos_validos), k_efectivo = k_efectivo))
  }
  
  puntaje_Dunn <- tryCatch(
    dunn(distance = submatriz, asignacion),
    error = function(e) { cat("  ERROR Dunn:", e$message, "\n"); NA }
  )
  
  puntaje_Connectivity <- tryCatch(
    connectivity(distance = submatriz, asignacion),
    error = function(e) { cat("  ERROR Connectivity:", e$message, "\n"); NA }
  )
  
  puntaje_Silhouette <- tryCatch({
    sil <- silhouette(asignacion, as.dist(submatriz))
    mean(sil[, 3])
  }, error = function(e) { cat("  ERROR Silhouette:", e$message, "\n"); NA })
  
  return(list(
    Dunn         = round(puntaje_Dunn, 3),
    Connectivity = round(puntaje_Connectivity, 3),
    Silhouette   = round(puntaje_Silhouette, 3),
    n_nodos      = length(nodos_validos),
    k_efectivo   = k_efectivo
  ))
}

# =============================================================================
# 3. FUNCIÓN AUXILIAR: guardar gráfico de red igraph como PNG
# =============================================================================
guardar_red_png <- function(mst_resultado, ruta_png, titulo) {
  
  tryCatch({
    
    png(filename = ruta_png, width = 1800, height = 1400, res = 150)
    
    igraph::V(mst_resultado$network)$label.cex <- 0.5
    membresia <- igraph::clusters(mst_resultado$network)$membership
    
    # Paleta de colores distinguible por cluster
    paleta <- rainbow(max(membresia))
    
    plot(
      mst_resultado$network,
      vertex.size   = 5,
      vertex.color  = paleta[membresia],
      vertex.label  = NA,
      vertex.frame.color = NA,
      edge.width    = 0.8,
      edge.color    = "gray60",
      layout        = igraph::layout.fruchterman.reingold(
        mst_resultado$network, niter = 10000),
      main          = titulo
    )
    
    # Leyenda de clusters
    legend("bottomleft",
           legend = paste("Cluster", seq_len(max(membresia)),
                          paste0("(n=", mst_resultado$csize, ")")),
           fill   = paleta,
           border = NA,
           cex    = 0.7,
           bty    = "n")
    
    dev.off()
    cat(sprintf("  Red guardada: %s\n", basename(ruta_png)))
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("  ERROR al guardar PNG '%s': %s\n", basename(ruta_png), e$message))
  })
}

# =============================================================================
# 4. EJECUTAR ESCENARIOS MST-kNN
# =============================================================================
set.seed(2)

escenarios <- list(
  list(nombre      = "Auto (sin suggested.k)",
       nombre_corto = "auto",
       suggested_k = NULL),
  list(nombre      = "suggested.k = 3",
       nombre_corto = "k3",
       suggested_k = 3),
  list(nombre      = "suggested.k = 5",
       nombre_corto = "k5",
       suggested_k = 5),
  list(nombre      = "suggested.k = 7",
       nombre_corto = "k7",
       suggested_k = 7),
  list(nombre      = "suggested.k = 10",
       nombre_corto = "k10",
       suggested_k = 10)
)

`%||%` <- function(a, b) if (!is.null(a)) a else b
resultados_escenarios <- list()
res_calidad_todos     <- data.frame()

for (esc in escenarios) {
  
  cat(sprintf("\n--- Ejecutando escenario: %s ---\n", esc$nombre))
  
  tiempo_inicio <- proc.time()
  
  mst_resultado <- tryCatch({
    if (is.null(esc$suggested_k)) {
      mst.knn(distance.matrix = matriz_cuadrada)
    } else {
      mst.knn(distance.matrix = matriz_cuadrada,
              suggested.k     = esc$suggested_k)
    }
  }, error = function(e) {
    cat(sprintf("  ERROR en mst.knn: %s\n", e$message))
    NULL
  })
  
  tiempo_esc <- proc.time() - tiempo_inicio
  
  if (is.null(mst_resultado)) {
    res_calidad_todos <- rbind(res_calidad_todos, data.frame(
      Escenario    = esc$nombre,
      Method       = "MST-kNN",
      k_encontrado = NA,
      n_nodos      = NA,
      k_efectivo   = NA,
      Dunn         = NA,
      Connectivity = NA,
      Silhouette   = NA,
      Tiempo_s     = round(as.numeric(tiempo_esc)[3], 2)
    ))
    next
  }
  
  k_encontrado <- mst_resultado$cnumber
  cat(sprintf("  Clusters encontrados : %d\n", k_encontrado))
  cat(sprintf("  Nodos en solución    : %d de %d\n",
              length(names(mst_resultado$cluster)), nrow(matriz_cuadrada)))
  cat(sprintf("  Tamaño de clusters   : %s\n",
              paste(mst_resultado$csize, collapse = " | ")))
  
  # Calcular índices
  indices <- calcular_indices(mst_resultado, matriz_cuadrada)
  
  fila_calidad <- data.frame(
    Escenario    = esc$nombre,
    Method       = "MST-kNN",
    k_encontrado = k_encontrado,
    n_nodos      = indices$n_nodos %||% length(names(mst_resultado$cluster)),
    k_efectivo   = indices$k_efectivo %||% k_encontrado,
    Dunn         = indices$Dunn,
    Connectivity = indices$Connectivity,
    Silhouette   = indices$Silhouette,
    Tiempo_s     = round(as.numeric(tiempo_esc)[3], 2)
  )
  
  res_calidad_todos <- rbind(res_calidad_todos, fila_calidad)
  
  # Guardar red PNG de ESTE escenario
  ruta_png_esc <- file.path(
    DIR_RESULTS,
    paste0("mstknn_red_", esc$nombre_corto, "_k", k_encontrado, ".png")
  )
  
  titulo_png <- paste0(
    "MST-kNN | Escenario: ", esc$nombre, "\n",
    "Clusters: ", k_encontrado,
    " | Nodos en solución: ", length(names(mst_resultado$cluster)),
    " de ", nrow(matriz_cuadrada)
  )
  
  cat("  Guardando red PNG...\n")
  guardar_red_png(mst_resultado, ruta_png_esc, titulo_png)
  
  resultados_escenarios[[esc$nombre]] <- list(
    resultado    = mst_resultado,
    calidad      = fila_calidad,
    tiempo       = tiempo_esc,
    ruta_png     = ruta_png_esc
  )
}

cat("\n\nResumen de escenarios MST-kNN:\n")
print(res_calidad_todos)

# =============================================================================
# 5. SELECCIONAR MEJOR ESCENARIO
# =============================================================================
calidad_valida <- res_calidad_todos[
  !is.na(res_calidad_todos$Silhouette) & res_calidad_todos$k_encontrado > 1, ]

if (nrow(calidad_valida) > 0) {
  idx_mejor    <- which.max(calidad_valida$Silhouette)
  mejor_nombre <- calidad_valida$Escenario[idx_mejor]
} else {
  mejor_nombre <- escenarios[[1]]$nombre
  cat("AVISO: Ningún escenario produjo índices válidos. Usando escenario automático.\n")
}

cat(sprintf("\nMejor escenario seleccionado: %s\n", mejor_nombre))

mst_optimo <- resultados_escenarios[[mejor_nombre]]$resultado
k_optimo   <- mst_optimo$cnumber

# =============================================================================
# 6. TABLAS DE RESULTADOS DEL MEJOR ESCENARIO
# =============================================================================
nodos_validos <- names(mst_optimo$cluster)
nodos_validos <- nodos_validos[nodos_validos %in% rownames(matriz_cuadrada)]

asignaciones_df <- data.frame(
  Arbol   = nodos_validos,
  Cluster = as.integer(mst_optimo$cluster[nodos_validos])
)
asignaciones_df <- asignaciones_df[order(asignaciones_df$Cluster), ]
rownames(asignaciones_df) <- NULL

tamanos_df <- data.frame(
  Cluster = seq_len(k_optimo),
  Tamano  = as.integer(mst_optimo$csize)
)

nodos_excluidos <- setdiff(rownames(matriz_cuadrada), names(mst_optimo$cluster))
excluidos_df    <- if (length(nodos_excluidos) > 0) {
  cat(sprintf("\nNodos excluidos de la solución MST-kNN: %d\n", length(nodos_excluidos)))
  data.frame(Arbol = nodos_excluidos, Cluster = NA)
} else {
  NULL
}

cat("\nDistribución de árboles por cluster (mejor escenario):\n")
print(tamanos_df)

# =============================================================================
# 9. REPORTE EXCEL
# =============================================================================
cat("Generando reporte Excel...\n")

wb <- createWorkbook()

# Hoja 1: Tiempos
tiempos_df <- data.frame(
  Escenario      = res_calidad_todos$Escenario,
  Tiempo_Total_s = res_calidad_todos$Tiempo_s
)

wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Tiempos",
                              titulo_tabla = "Tiempos de Ejecución por Escenario (Segundos)",
                              datos        = tiempos_df,
                              anchos_col   = "auto")

# Hoja 2: Calidad de los 3 escenarios
wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Calidad_Clusters",
                              titulo_tabla = paste0("Índices de Calidad MST-kNN (3 escenarios) — Mejor: ",
                                                    mejor_nombre),
                              datos        = res_calidad_todos,
                              anchos_col   = "auto")

# Hoja 3: Asignaciones del mejor escenario
wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Asignaciones",
                              titulo_tabla = paste0("Asignaciones — k=", k_optimo,
                                                    " | ", mejor_nombre,
                                                    " | Nodos: ", nrow(asignaciones_df)),
                              datos        = asignaciones_df,
                              anchos_col   = "auto")

# Hoja 4: Tamaños de clusters
wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Tamanos_Clusters",
                              titulo_tabla = paste0("Tamaño de Clusters — ", mejor_nombre),
                              datos        = tamanos_df,
                              anchos_col   = "auto")

# Hoja 5: Nodos excluidos (si los hay)
if (!is.null(excluidos_df)) {
  wb <- agregar_hoja_formateada(wb           = wb,
                                nombre_hoja  = "Nodos_Excluidos",
                                titulo_tabla = paste0("Nodos excluidos de la solución MST-kNN (",
                                                      nrow(excluidos_df), " árboles)"),
                                datos        = excluidos_df,
                                anchos_col   = "auto")
}

# Hoja 7: Índice de redes PNG generadas
if (length(resultados_escenarios) > 0) {
  redes_generadas_df <- do.call(rbind, lapply(names(resultados_escenarios), function(nm) {
    data.frame(
      Escenario    = nm,
      k_encontrado = resultados_escenarios[[nm]]$calidad$k_encontrado,
      Archivo_PNG  = basename(resultados_escenarios[[nm]]$ruta_png)
    )
  }))
  
  wb <- agregar_hoja_formateada(wb           = wb,
                                nombre_hoja  = "Redes_PNG",
                                titulo_tabla = "Archivos PNG de Redes Generados por Escenario",
                                datos        = redes_generadas_df,
                                anchos_col   = "auto")
} else {
  cat("AVISO: No se generaron redes PNG (todos los escenarios fallaron).\n")
}

ruta_excel <- file.path(DIR_RESULTS, "mstknn_resultados.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)

cat("Resultados completos guardados en:", ruta_excel, "\n")
cat(sprintf("\n=== Resumen final MST-kNN ===\n"))
cat(sprintf("  Mejor escenario  : %s\n", mejor_nombre))
cat(sprintf("  k encontrado     : %d\n", k_optimo))
cat(sprintf("  Nodos en solución: %d de %d\n",
            length(nodos_validos), nrow(matriz_cuadrada)))
if (nrow(calidad_valida) > 0) {
  cat(sprintf("  Silhouette       : %.3f\n",
              max(calidad_valida$Silhouette, na.rm = TRUE)))
}
cat(sprintf("\nPNGs generados:\n"))
for (nm in names(resultados_escenarios)) {
  cat(sprintf("  - %s\n", basename(resultados_escenarios[[nm]]$ruta_png)))
}
