# =============================================================================
# COMPARACIÓN MEDIOIDE VS TODOS LOS ÁRBOLES
# =============================================================================
# Recibe como insumo:
#   - rf_matrix_normalizada.rds  (matriz de distancias precomputada)
#   - ranking_medioide.xlsx      (para identificar el nombre del medioide)
#   - los archivos ZIP originales (para releer los árboles)
# =============================================================================

library(ape)
library(TreeDist)
library(openxlsx)

# --- Cargar configuración y funciones globales ---
source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# =============================================================================
# Recuperar el árbol medioide desde los resultados anteriores
# =============================================================================
cat("=== CARGANDO MEDIOIDE ===\n")

# -- Leer el ranking para obtener el nombre del medioide
ranking <- read.xlsx(file.path(DIR_RESULTS, paste0("ranking_medioide_", NOMBRE_BDD, ".xlsx")), sheet = "Ranking_Medioide", startRow = 2)
nombre_medioide <- ranking$nombre_arbol[ranking$posicion == 1]
cat("Árbol medioide identificado:", nombre_medioide, "\n")

# -- Releer todos los árboles desde los ZIP
arboles <- leer_bosque_zip(directorio = file.path(DIR_INPUT, CARPETA_ARBOLES),
                           ext_interna = EXTENSION_ARBOLES,
                           dir_cache   = DIR_CACHE
                           ) 

# -- Extraer el medioide como tree1 (fijo para todas las comparaciones)
tree1 <- arboles[[nombre_medioide]]
cat("Tips en medioide:", length(tree1$tip.label), "\n")

# =============================================================================
# Función de unión por especie ancla
# =============================================================================
tiempo_funcion_ancla <- system.time({
  unir_por_ancla <- function(tree1, tree2) {
    
    especies_comunes <- intersect(tree1$tip.label, tree2$tip.label)
    
    if (length(especies_comunes) < 2) {
      return(list(arbol = NULL,
                  n_comunes = length(especies_comunes),
                  ancla = NA,
                  error = "Menos de 2 especies comunes"))
    }
    
    especie_ancla <- especies_comunes[1]
    otras_comunes <- especies_comunes[especies_comunes != especie_ancla]
    
    tree1_preparado <- drop.tip(tree1, otras_comunes)
    posicion_final  <- which(tree1_preparado$tip.label == especie_ancla)
    
    arbol_final <- tryCatch({
      resultado <- bind.tree(tree1_preparado, tree2, where = posicion_final)
        collapse.singles(resultado)
    }, error = function(e) {
      return(NULL)
    })
    
    return(list(
      arbol     = arbol_final,
      n_comunes = length(especies_comunes),
      ancla     = especie_ancla,
      error     = ifelse(is.null(arbol_final), "Error en bind.tree", NA)
    ))
  }
})


# =============================================================================
# Comparación medioide vs todos los árboles
# =============================================================================
cat("\n=== COMPARACIÓN MEDIOIDE VS TODOS LOS ÁRBOLES ===\n")

nombres_arboles <- names(arboles)
n_arboles       <- length(arboles)

# Estructuras para guardar resultados
arboles_unidos  <- vector("list", n_arboles)
names(arboles_unidos) <- nombres_arboles

log_resultados <- data.frame(
  posicion       = seq_len(n_arboles),
  nombre         = nombres_arboles,
  es_medioide    = nombres_arboles == nombre_medioide,
  n_tips_tree2   = NA_integer_,
  n_comunes      = NA_integer_,
  ancla          = NA_character_,
  tips_union     = NA_integer_,
  estado         = NA_character_,
  stringsAsFactors = FALSE
)

tiempo_comparacion <- system.time({
  for (i in seq_len(n_arboles)) {
    
    tree2       <- arboles[[i]]
    nombre_tree2 <- nombres_arboles[i]
    
    # El medioide comparado consigo mismo se registra pero no se une
    if (nombre_tree2 == nombre_medioide) {
      log_resultados$n_tips_tree2[i] <- length(tree2$tip.label)
      log_resultados$n_comunes[i]    <- length(tree2$tip.label)
      log_resultados$estado[i]       <- "MEDIOIDE (referencia)"
      arboles_unidos[[i]]            <- tree1
      cat(" [", i, "/", n_arboles, "]", nombre_tree2, "→ es el medioide (referencia)\n")
      next
    }
    
    resultado <- unir_por_ancla(tree1, tree2)
    
    log_resultados$n_tips_tree2[i] <- length(tree2$tip.label)
    log_resultados$n_comunes[i]    <- resultado$n_comunes
    log_resultados$ancla[i]        <- resultado$ancla
    log_resultados$estado[i]       <- ifelse(is.na(resultado$error), "OK", resultado$error)
    
    if (!is.null(resultado$arbol)) {
      arboles_unidos[[i]]          <- resultado$arbol
      log_resultados$tips_union[i] <- length(resultado$arbol$tip.label)
    } else {
      cat(" [", i, "/", n_arboles, "]", nombre_tree2,
          "→ ERROR:", resultado$error, "\n")
    }
  }
})

# =============================================================================
# Resumen y validación
# =============================================================================
cat("\n=== RESUMEN ===\n")

n_ok     <- sum(log_resultados$estado == "OK", na.rm = TRUE)
n_error  <- sum(!is.na(log_resultados$error) & log_resultados$estado != "OK" &
                  log_resultados$estado != "MEDIOIDE (referencia)")

cat("Comparaciones exitosas :", n_ok, "\n")
cat("Comparaciones con error:", n_error, "\n")
cat("Tips promedio en unión  :",
    round(mean(log_resultados$tips_union, na.rm = TRUE), 1), "\n")
cat("Tips mínimo en unión    :", min(log_resultados$tips_union, na.rm = TRUE), "\n")
cat("Tips máximo en unión    :", max(log_resultados$tips_union, na.rm = TRUE), "\n")


# =============================================================================
# Guardar cada árbol unido como .rootree.zip individual
# =============================================================================

cat("\n=== GUARDANDO ÁRBOLES PODADOS ===\n")

dir_arboles <- file.path(DIR_PROCESSED, "arboles_podados")
if (!dir.exists(dir_arboles)) {
  dir.create(dir_arboles, recursive = TRUE)
  cat("Carpeta creada:", dir_arboles, "\n")
}

guardados <- 0
fallidos_guardado <- c()

tiempo_guardado <- system.time({
  for (i in seq_len(n_arboles)) {
    
    if (is.null(arboles_unidos[[i]])) next
    
    nombre_base <- nombres_arboles[i]
    arbol_actual <- arboles_unidos[[i]]
    
    # Ruta del .rootree temporal y del .zip final
    ruta_rootree <- file.path(dir_arboles, paste0(nombre_base))
    ruta_zip     <- file.path(dir_arboles, paste0(nombre_base, ".zip"))
    
    exito <- tryCatch({
      # Escribir árbol en formato Newick como .rootree
      write.tree(arbol_actual, file = ruta_rootree)
      
      # Comprimir en ZIP manteniendo el nombre interno como .rootree
      zip(zipfile = ruta_zip,
          files   = ruta_rootree,
          flags   = "-j")  # -j: no incluir rutas, solo el archivo
      
      # Eliminar el .rootree temporal (ya está dentro del ZIP)
      file.remove(ruta_rootree)
      
      TRUE
    }, error = function(e) {
      message("  Error guardando ", nombre_base, ": ", e$message)
      FALSE
    })
    
    if (exito) {
      guardados <- guardados + 1
    } else {
      fallidos_guardado <- c(fallidos_guardado, nombre_base)
    }
  }
})


cat("\n=== RESUMEN GUARDADO ===\n")
cat("Árboles guardados :", guardados, "\n")
cat("Fallos            :", length(fallidos_guardado), "\n")
cat("Carpeta           :", dir_arboles, "\n")

if (length(fallidos_guardado) > 0) {
  cat("No guardados:\n")
  cat(paste(" ", fallidos_guardado, collapse = "\n"), "\n")
}

# =============================================================================
# Guardar resultados en Excel
# =============================================================================
cat("\n=== GUARDANDO RESULTADOS ===\n")

# ── Hoja para Tiempos ───────────────────────

# Preparar el dataframe base
t_ancla <- as.numeric(tiempo_funcion_ancla)
t_comp  <- as.numeric(tiempo_comparacion)
t_guardado <- as.numeric(tiempo_guardado)

tiempos_df <- data.frame(
  Proceso          = c("2. Definición función ancla",
                       "3. Comparación medioide vs todos",
                       "4. Guardado de árboles podados"),
  Tiempo_Usuario_s = c(t_ancla[1], t_comp[1], t_guardado[1]),
  Tiempo_Sistema_s = c(t_ancla[2], t_comp[2], t_guardado[2]),
  Tiempo_Total_s   = c(t_ancla[3], t_comp[3], t_guardado[3])
)

# Añadir la fila de totales directamente en R para que la función global la formatee
fila_total <- data.frame(
  Proceso          = "TOTAL",
  Tiempo_Usuario_s = sum(tiempos_df$Tiempo_Usuario_s),
  Tiempo_Sistema_s = sum(tiempos_df$Tiempo_Sistema_s),
  Tiempo_Total_s   = sum(tiempos_df$Tiempo_Total_s)
)
tiempos_df <- rbind(tiempos_df, fila_total)

# ── Hoja para Log de Comparaciones ──────────────────────────
# (Asumiendo que el libro se crea aquí en el entorno, si es un archivo
# físico existente que se quiere leer, cambiar createWorkbook() por 
# loadWorkbook("ruta/al/archivo.xlsx"))

wb <- createWorkbook()

# Aplicar la función global
agregar_hoja_formateada(wb, "Log_Comparaciones",
                        paste("Log de Comparaciones Medioide -", NOMBRE_BDD),
                        log_resultados
                        )

agregar_hoja_formateada(wb, "Tiempos_Ejecucion",
                        paste("Tiempos de Ejecución Script 03 -", NOMBRE_BDD),
                        tiempos_df
                        )

ruta_comparaciones_output <- file.path(DIR_RESULTS, paste0("comparaciones_medioide_", NOMBRE_BDD, ".xlsx"))
saveWorkbook(wb, ruta_comparaciones_output, overwrite = TRUE)

cat("Log de comparaciones guardado en:", ruta_comparaciones_output, "\n")

cat("\n=== COMPLETADA ===\n")
#_____________________________________ TEST ___________________________________________
# ── Cambiar este valor para explorar distintos árboles ────────────
ARBOL_A_VER <- "56139_NT_AL.rootree"   # índice numérico, o por nombre: "nombre_del_arbol"

# ── Obtener los árboles ───────────────────────────────────────────
tree2_viz <- arboles[[ARBOL_A_VER]]
nombre_viz <- nombres_arboles[ARBOL_A_VER]

resultado_viz <- unir_por_ancla(tree1, tree2_viz)

# ── Visualización lado a lado ─────────────────────────────────────
par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))

# Panel 1: medioide (tree1)
plot(tree1,
     type           = "fan",
     show.tip.label = FALSE,
     edge.color     = "steelblue",
     edge.width     = 0.8,
     main           = paste0("MEDIOIDE\n(", length(tree1$tip.label), " tips)"),
     cex.main       = 0.9)

# Panel 2: árbol a comparar (tree2)
plot(tree2_viz,
     type           = "fan",
     show.tip.label = FALSE,
     edge.color     = "darkorange",
     edge.width     = 0.8,
     main           = paste0(nombre_viz, "\n(", length(tree2_viz$tip.label), " tips)"),
     cex.main       = 0.9)

# Panel 3: árbol unido
if (!is.null(resultado_viz$arbol)) {
  plot(resultado_viz$arbol,
       type           = "fan",
       show.tip.label = FALSE,
       edge.color     = "forestgreen",
       edge.width     = 0.8,
       main           = paste0("UNIÓN\n(", length(resultado_viz$arbol$tip.label), " tips)",
                               "\nAncla: ", resultado_viz$ancla),
       cex.main       = 0.9)
} else {
  plot.new()
  text(0.5, 0.5, paste("ERROR:\n", resultado_viz$error), cex = 1.2, col = "red")
}

par(mfrow = c(1, 1))

# ── Resumen en consola ────────────────────────────────────────────
cat("=== ÁRBOL:", nombre_viz, "===\n")
cat("Tips tree2  :", length(tree2_viz$tip.label), "\n")
cat("Tips unión  :", ifelse(!is.null(resultado_viz$arbol),
                            length(resultado_viz$arbol$tip.label), "ERROR"), "\n")
cat("N° comunes  :", resultado_viz$n_comunes, "\n")
cat("Especie ancla:", resultado_viz$ancla, "\n")

cat("\n=== COMPLETADA 2 ===\n")