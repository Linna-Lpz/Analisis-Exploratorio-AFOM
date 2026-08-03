# ==============================================================================
# ARCHIVO: config.R
# PROPÓSITO: Variables globales, rutas y setup de infraestructura
# ==============================================================================
library(here) 

# 1. Variable de Control Principal
NOMBRE_BDD <- "OrthoMaM" # Nombre de la base de datos
EXTENSION_ARBOLES <- ".rootree"
CARPETA_ARBOLES <- "archiveTreesV12" # Nombre original de la carpeta con árboles

# 2. Definición de Rutas Base 
DIR_INPUT     <- here("Datos", "Original")
DIR_PROCESSED <- here("Datos", "Procesados")     
DIR_RESULTS   <- here("Resultados")
DIR_SCRIPTS   <- here("Scripts")

DIR_CACHE <- file.path(DIR_RESULTS, "cache")

# ------------------------------------------------------------------------------
# 3. Validación e Infraestructura
# ------------------------------------------------------------------------------
carpetas_necesarias <- c(DIR_PROCESSED, DIR_RESULTS, DIR_CACHE)

for (carpeta in carpetas_necesarias) {
  if (!dir.exists(carpeta)) {
    dir.create(carpeta, recursive = TRUE)
    message("Infraestructura: Carpeta creada -> ", carpeta)
  }
}

# Mensaje de confirmación final
message("Configuración cargada y entorno preparado para: ", NOMBRE_BDD)