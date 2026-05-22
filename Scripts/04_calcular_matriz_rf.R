# =============================================================================
# CÁLCULO DE MATRIZ RF NORMALIZADO — CON CACHÉ .rds
# =============================================================================
library(ape)
library(TreeDist)
library(openxlsx)

source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# --- Directorio de caché ---
DIR_CACHE <- file.path(DIR_PROCESSED, "cache")
if (!dir.exists(DIR_CACHE)) dir.create(DIR_CACHE, recursive = TRUE)

ruta_cache_matriz <- file.path(DIR_CACHE, "matriz_rf.rds")

# =============================================================================
# LEER ÁRBOLES
# =============================================================================
cat("Iniciando lectura de árboles...\n")
tiempo_lectura <- system.time({
  arboles_conjunto <- leer_bosque_zip(
    directorio   = file.path(DIR_PROCESSED, "arboles_podados"),
    ext_interna  = EXTENSION_ARBOLES,
    dir_cache   = DIR_CACHE
  )
})

# =============================================================================
# CALCULAR MATRIZ RF — con caché condicional
# =============================================================================
if (file.exists(ruta_cache_matriz)) {
  cat("Matriz RF encontrada en caché. Cargando...\n")
  tiempo_matriz <- system.time({
    matriz_cuadrada <- readRDS(ruta_cache_matriz)
  })
  cat("Cargada desde caché:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")
  
} else {
  cat("Calculando matriz Robinson-Foulds...\n")
  tiempo_matriz <- system.time({
    matriz_distancias <- RobinsonFoulds(arboles_conjunto, normalize = TRUE)
    matriz_cuadrada   <- as.matrix(matriz_distancias)
    # Asignar nombres de arboles para la matriz
    nombres <- names(arboles_conjunto)
    rownames(matriz_cuadrada) <- nombres
    colnames(matriz_cuadrada) <- nombres
    saveRDS(matriz_cuadrada, file = ruta_cache_matriz)            # <- Guardar caché
  })
  cat("Matriz calculada y guardada en caché:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")
}

# =============================================================================
# EXPORTAR TAMBIÉN COMO CSV (opcional, para inspección externa)
# =============================================================================
#ruta_csv <- file.path(DIR_RESULTS, "matriz_rf_conjunto.csv")
#write.table(matriz_cuadrada,
#            file      = ruta_csv,
#            sep       = ";",
#            row.names = TRUE,
#            col.names = NA,
#            quote     = FALSE)
#cat("Copia CSV guardada en:", ruta_csv, "\n")

# =============================================================================
# REPORTE DE TIEMPOS EN EXCEL
# =============================================================================
t_lectura <- as.numeric(tiempo_lectura)
t_matriz  <- as.numeric(tiempo_matriz)

tiempos_df <- data.frame(
  Proceso          = c("1. Lectura de árboles", "2. Cálculo/carga matriz RF", "TOTAL"),
  Tiempo_Usuario_s = c(t_lectura[1], t_matriz[1], t_lectura[1] + t_matriz[1]),
  Tiempo_Sistema_s = c(t_lectura[2], t_matriz[2], t_lectura[2] + t_matriz[2]),
  Tiempo_Total_s   = c(t_lectura[3], t_matriz[3], t_lectura[3] + t_matriz[3])
)

wb <- createWorkbook()
wb <- agregar_hoja_formateada(wb           = wb,
                              nombre_hoja  = "Tiempos",
                              titulo_tabla = "Reporte de Tiempos de Ejecución (Segundos)",
                              datos        = tiempos_df,
                              anchos_col   = "auto")

ruta_excel <- file.path(DIR_RESULTS, "tiempos_matriz_rf.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Tiempos guardados en:", ruta_excel, "\n")