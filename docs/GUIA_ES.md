# Guía de uso en español

## 1. Instalación y pruebas con datos sintéticos

1. Descomprime el archivo y abre en RStudio la carpeta que contiene `run_analysis.R` y `README.md`.
2. Ejecuta los comandos siguientes desde la raíz del proyecto:

```r
source("scripts/install_dependencies.R")
source("tests/run_tests.R")
source("tests/run_annotation_tests.R")
source("tests/run_visualization_tests.R")
source("run_analysis.R")
source("scripts/run_functional_analysis.R")
```

3. Comprueba que las pruebas concluyen sin errores y revisa los CSV, Excel y gráficos generados en `results/`. Los datos del ejemplo son sintéticos.

## 2. Análisis de datos propios

- Consulta los CSV de `data/example/` y prepara los datos con las mismas columnas.
- Modifica las rutas y parámetros en `config/settings.R` y, si procede, la configuración del preprocesamiento PEAKS.
- Mantén datos originales y resultados de investigaciones no publicadas fuera del repositorio compartido.
- Comprueba que los identificadores de muestras de `sample_metadata.csv` coinciden con las columnas de abundancia.
- Consulta [la metodología](SCIENTIFIC_METHODS.md), [los ejemplos](EXECUTION_EXAMPLES.md) y [la validación](VALIDATION.md) para conocer el alcance estadístico y sus limitaciones.

## 3. Interpretación y citación

- Inspecciona los resultados, diagnósticos y figuras antes de formular conclusiones biológicas.
- Los resúmenes de categorías funcionales son descriptivos; no equivalen a análisis de enriquecimiento.
- Los metadatos de citación se encuentran en `CITATION.cff` y la licencia en `LICENSE`.
- Para obtener la información de la versión archivada, consulta las Releases del repositorio y su DOI correspondiente cuando esté disponible.
