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
- Se registra de manera automatizada el uso de RAM antes y después de cada paso.
- El log completo de consumo de memoria y tiempos de ejecución se exporta automáticamente en formatos Excel y CSV dentro de la carpeta `Resultados/`.

## 4. Descripción de los Scripts (Pipeline)

A continuación se detalla la función de cada script del directorio `Scripts/`, junto con sus insumos clave (datos de entrada) y artefactos generados (datos de salida).

| Script | Descripción | Datos de Entrada | Datos de Salida |
|--------|-------------|------------------|-----------------|
| `00_funciones_globales.R` | Funciones utilitarias (exportar a excel, RAM, logs). | N/A | Funciones globales en entorno |
| `01_analisis_por_arbol.R` | Lee, preprocesa y filtra los árboles filogenéticos. | Archivos Newick originales | Metadatos y variables de árboles |
| `02_generar_arbol_medioide.R` | Calcula el árbol medioide como centroide global. | Árboles preprocesados | Árbol Medioide |
| `03_comparar_medioide_vs_arboles.R`| Cuantifica distancias de cada árbol al medioide. | Árboles y Medioide | Ranking y distancias al medioide |
| `04_calcular_matriz_rf.R` | Computa la matriz densa de distancias Robinson-Foulds. | Árboles injertados | Matriz RF, caché local (`.rds`) |
| `05a_reduccion_dimensional_umap.R` | Proyección no lineal UMAP del espacio topológico. | Matriz RF o Coordenadas PCA | Coordenadas UMAP, Visualización |
| `05c_reduccion_dimensional_tSNE.R` | Proyección no lineal t-SNE del espacio. | Matriz RF o Coordenadas PCA | Coordenadas t-SNE, Visualización |
| `05d_reduccion_dimensional_PCA.R` | Análisis de Componentes Principales clásico (PCA). | Matriz RF | Coordenadas PCA, Scree Plot |
| `06a_clustering_kmeans.R` | Agrupamiento particional K-Means sobre PCA/MDS. | Coordenadas reducidas | Resultados de partición (K-Means) |
| `06b_clustering_pam.R` | Agrupamiento PAM directo sobre distancias métricas. | Matriz RF | Resultados de partición (PAM) |
| `06c_clustering_clara.R` | Agrupamiento escalable CLARA sobre MDS/PCA. | Coordenadas reducidas | Resultados de partición (CLARA) |
| `06d_clustering_mstknn.R` | Agrupamiento topológico en grafos (MST-kNN). | Matriz RF | Red y particiones MST-kNN |
| `06e_comparativa_silhouette_rf.R` | Re-evalúa Silueta y Dunn de algoritmos previos. | Clústeres y Matriz RF | Tabla comparativa de métricas |
| `07_iterar_[algoritmo].R` | Ejecuta subdivisión recursiva de clústeres. | Clústeres base | Clústeres finales e iterados |
| `08_etiquetas.R` | Mapeo de identificadores NCBI a símbolos HGNC. | Clústeres finales | Archivos de genes por clúster |
| `09_enriquecimiento.R` | Análisis funcional biológico frente a Gene Ontology. | Listas HGNC | Tabla de términos enriquecidos |
| `10_grafico_enriquecimiento.R` | Visualización de términos GO más significativos. | Tabla de enriquecimiento | Gráfico de burbujas GO |
| `11_mapa_de_calor.R` | Visualización integradora funcional (Heatmap). | Tabla de enriquecimiento | Heatmap cruzado de clústeres |
| `12_validacion_nula_GO.R` | Modelo nulo para validación estocástica. | Listas HGNC de genes | Tabla y p-valor empírico, Histograma |
| `13_escalabilidad_subconjuntos.R` | Mide tiempo y uso de memoria del pipeline sobre subconjuntos crecientes de árboles, para evaluar escalabilidad. | Árboles originales (bosque completo) | `curvas_escalabilidad.csv` |
| `14_verificacion_PC1_vs_hojas.R` | Verifica si PC1 refleja variabilidad topológica real o solo diferencias en cobertura taxonómica (número de hojas) por árbol. | Resultado PCA (`pca_resultado.rds`), Matriz RF (`matriz_rf.rds`), bosque en caché | Gráfico de correlación PC1 vs. n° hojas, CSV con scores y análisis |
| `15_control_324_arboles_completos.R` | Repite el pipeline de clustering y enriquecimiento restringido al subconjunto "core set" de los árboles con cobertura taxonómica completa, como control frente al sesgo del injerto. | Conjunto core set (`conjunto_core_set.rds`) | Matriz RF subconjunto × subconjunto×324, clustering, gráfico de burbujas GO, tablas de enriquecimiento y clústeres (`.xlsx`, `.png`) |
| `16_analisis_sensibilidad.R` | Analiza la sensibilidad del medioide y de la partición, injertando a partir del segundo medioide y comparando (vía ARI) contra la partición original. | Ranking del medioide (`ranking_medioide_*.xlsx`), partición final (`metodoClustering_subdivision_iterativa.xlsx`), árboles originales | `analisis_sensibilidad_medioide_*.xlsx` |
| `17_comparacion_ari_smith.R` | Compara la distancia Robinson-Foulds con la métrica de información `ClusteringInfoDistance`, reagrupando el core set con MST-kNN y calculando el ARI entre ambas particiones. | Conjunto core set (`conjunto_core_set.rds`), clustering de control (`control_324_clustering_*.xlsx`) | `control_324_smith_ari_*.xlsx` |
| `18_calcular_BHI.R` | Índice de Homogeneidad Biológica empírico (BHI). | Clústeres e identificadores | Valor cuantitativo de BHI |
| `config.R` | Define variables globales del pipeline (nombre de la base de datos, extensión y carpeta de árboles), rutas base (`Datos/Original`, `Datos/Procesados`, `Resultados`, `Resultados/cache`, `Scripts`) y crea automáticamente los directorios necesarios si no existen. | N/A | Variables y rutas globales en el entorno; estructura de carpetas creada |
| `evaluacion_sesgo.R` | Analiza el impacto analítico del injerto en árboles. | Topologías pre y post injerto | Gráfico de caja (sesgo evaluado) |
| `instalar_paquetes.R` | Verifica e instala todas las dependencias del pipeline (paquetes CRAN y Bioconductor: árboles filogenéticos, clustering, reducción dimensional, visualización, Excel, enriquecimiento funcional, paralelización). Se ejecuta una sola vez al configurar el entorno. | N/A | Paquetes R instalados en el entorno |
| `restaurar_cache_kmeans` | Reconstruye manualmente el archivo de caché de calidad de K-Means (`kmeans_calidad_res.rds`) a partir de valores registrados en consola, evitando repetir hrs de cómputo en el servidor. | Valores de calidad (Dunn, Connectivity, Silhouette) y tiempos de ejecución transcritos del log de consola | `kmeans_calidad_res.rds` (caché de resultados y tiempos de K-Means) |
| `run_pipeline.R` | Orquestador principal del KDD y dependencias. | `config.R`, Scripts | Pipeline completo procesado |
