# =============================================================================
# CÁLCULO DE MATRIZ RF NORMALIZADO A TODO EL CONJUNTO DE ÁRBOLES
# =============================================================================

library(ape)
library(TreeDist)
library(openxlsx)

# --- Cargar configuración y funciones globales ---
# Se asume que config.R define DIR_RAW y DIR_RESULTS
source(here::here("Scripts", "config.R"))
source(here::here("Scripts", "00_funciones_globales.R"))

# 1. Leer árboles usando la función global
cat("Iniciando lectura de árboles...\n")
tiempo_lectura <- system.time({
  arboles_conjunto <- leer_bosque_zip(directorio = file.path(DIR_PROCESSED, "arboles_podados"),
                                      ext_interna = EXTENSION_ARBOLES
                                      )
})

# 2. Calcular la matriz de distancias Robinson-Foulds normalizada
cat("Calculando matriz Robinson-Foulds...\n")
tiempo_matriz <- system.time({
  matriz_distancias <- RobinsonFoulds(arboles_conjunto, normalize = TRUE)
  matriz_cuadrada   <- as.matrix(matriz_distancias) # Convertir para exportar a CSV
})

# 3. Exportar matriz como CSV (más liviano para matrices grandes)
ruta_csv <- file.path(DIR_RESULTS, "matriz_rf_conjunto.csv")
write.table(matriz_cuadrada,
            file      = ruta_csv,
            sep       = ";",
            row.names = TRUE,
            col.names = NA,
            quote     = FALSE)
cat("Matriz guardada como CSV:", nrow(matriz_cuadrada), "x", ncol(matriz_cuadrada), "\n")


# 4. Reporte de Tiempos en Excel usando la función global
cat("Generando reporte de tiempos...\n")
t_lectura <- as.numeric(tiempo_lectura)
t_matriz  <- as.numeric(tiempo_matriz)

# Armamos el dataframe e incluimos la fila TOTAL calculada nativamente en R
tiempos_df <- data.frame(
  Proceso          = c("1. Lectura de árboles", "2. Cálculo matriz RF", "TOTAL"),
  Tiempo_Usuario_s = c(t_lectura[1], t_matriz[1], t_lectura[1] + t_matriz[1]),
  Tiempo_Sistema_s = c(t_lectura[2], t_matriz[2], t_lectura[2] + t_matriz[2]),
  Tiempo_Total_s   = c(t_lectura[3], t_matriz[3], t_lectura[3] + t_matriz[3])
)

# Inicializar workbook y aplicar la función de formato global
wb <- createWorkbook()
wb <- agregar_hoja_formateada(wb           = wb, 
                              nombre_hoja  = "Tiempos", 
                              titulo_tabla = "Reporte de Tiempos de Ejecución (Segundos)", 
                              datos        = tiempos_df, 
                              anchos_col   = "auto")

# Guardar el Excel
ruta_excel <- file.path(DIR_RESULTS, "tiempos_matriz_rf.xlsx")
saveWorkbook(wb, ruta_excel, overwrite = TRUE)
cat("Tiempos guardados exitosamente en:", ruta_excel, "\n")