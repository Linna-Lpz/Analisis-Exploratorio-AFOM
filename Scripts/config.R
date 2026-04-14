# ==============================================================================
# ARCHIVO: config.R
# PROPÓSITO: Variables globales, rutas y setup de infraestructura
# ==============================================================================
library(here) 

# 1. Variable de Control Principal
ORIGEN_DATOS <- "OrthoMaM"
EXTENSION_ARBOLES <- ".rootree"

# 2. Definición de Rutas Base 
DIR_INPUT     <- here("Datos", "Original")
DIR_PROCESSED <- here("Datos", "Procesados")     # Matrices y coordenadas
DIR_RESULTS   <- here("Resultados", "Figuras")    # Gráficos
DIR_SCRIPTS   <- here("Scripts")

# 3. Generación de Rutas Dinámicas 
ruta_arboles_input <- file.path(DIR_INPUT, paste0(ORIGEN_DATOS, "_arboles", EXTENSION_ARBOLES))
ruta_matriz_output <- file.path(DIR_PROCESSED, paste0("matriz_RF_", ORIGEN_DATOS, ".rds"))
ruta_coord_output  <- file.path(DIR_PROCESSED, paste0("coordenadas_MDS_", ORIGEN_DATOS, ".rds"))

# ------------------------------------------------------------------------------
# 4. Validación e Infraestructura (Self-healing)
# ------------------------------------------------------------------------------
carpetas_necesarias <- c(DIR_INPUT, DIR_PROCESSED, DIR_RESULTS, DIR_SCRIPTS)

for (carpeta in carpetas_necesarias) {
  if (!dir.exists(carpeta)) {
    dir.create(carpeta, recursive = TRUE)
    message("Infraestructura: Carpeta creada -> ", carpeta)
  }
}

# Mensaje de confirmación final
message("Configuración cargada y entorno preparado para: ", ORIGEN_DATOS)