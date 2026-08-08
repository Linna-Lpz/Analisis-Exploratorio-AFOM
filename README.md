# Ejecución del Pipeline Analítico

Este documento detalla las instrucciones y condiciones necesarias para ejecutar el pipeline completo de análisis filogenómico a través del script maestro `run_pipeline.R`.

## PIPELINE MAESTRO — Ejecuta todos los scripts en orden secuencial

### USO
Puedes ejecutar el pipeline de dos maneras:
- **Opción A (RStudio Server/Desktop):** Abre el archivo `run_pipeline.R` en el IDE y pulsa el botón **"Source"**.
- **Opción B (Terminal/Línea de comandos):** Ejecuta el siguiente comando:
  ```bash
  Rscript Scripts/run_pipeline.R
  ```

### REQUISITOS PREVIOS
- El proyecto debe abrirse **estrictamente desde el archivo `.Rproj`** (para garantizar que el paquete `here()` resuelva las rutas relativas correctamente).
- Los datos originales (el corpus de árboles filogenéticos) deben estar descomprimidos y ubicados en:
  `Datos/Original/archiveTreesV12/`
- Todos los paquetes de R requeridos deben estar instalados. Puedes verificar o instalar dependencias corriendo el script auxiliar:
  `Scripts/instalar_paquetes.R`

### CONTROL DE EJECUCIÓN
El script está diseñado de manera modular y se divide en dos fases: los **Pasos Base** (reducción dimensional y clustering general) y los **Pasos Downstream** (enriquecimiento funcional y validación).
En la sección de configuración al inicio de `run_pipeline.R` puedes ajustar:

1. **Pasos a ejecutar (`PASOS_EJECUTAR`)**: Una lista booleana (`TRUE`/`FALSE`). Cambia a `FALSE` los pasos que **NO** quieras re-ejecutar (útil si ya procesaste esa etapa y cuentas con los archivos en caché).
2. **Selección Dinámica (`ALGORITMO_DOWNSTREAM`)**: 
   - `"AUTO"`: El pipeline evalúa empíricamente todos los agrupamientos base y escoge automáticamente el mejor (el de mayor índice Silhouette) para continuar con el enriquecimiento.
   - `"MST-kNN"`, `"CLARA"`, etc.: Fuerza la ejecución del downstream para un algoritmo específico.
3. **Manejo de Errores (`DETENER_EN_ERROR`)**: Si está en `TRUE` y un script falla por cualquier motivo, la orquestación se detiene de inmediato y se imprime un informe del error para facilitar la depuración, protegiendo así los resultados previos.

### GESTIÓN DE MEMORIA
El pipeline está optimizado para el manejo de memoria:
- Cada script individual se ejecuta en un **entorno local aislado** (`eval()` en un nuevo `env`). Esto garantiza que los objetos pesados se destruyan y el recolector de basura (`gc()`) libere la RAM al finalizar cada etapa, previniendo crashes por memoria agotada (OOM) en el servidor.
- Se registra de manera automatizada el uso de RAM antes y después de cada paso.
- El log completo de consumo de memoria y tiempos de ejecución se exporta automáticamente en formatos Excel y CSV dentro de la carpeta `Resultados/`.
