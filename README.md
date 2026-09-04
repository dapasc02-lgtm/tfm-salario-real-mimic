# Pipeline de estimación del salario real mediante modelo MIMIC

Código del Trabajo de Fin de Máster *"Uso inteligente de datos administrativos: 
un enfoque eficiente de estimación del salario real con modelos de ecuaciones 
estructurales"* — Máster en Ingeniería Matemática, Universidad Complutense de Madrid.

Autor: Daniel Pascual Martínez
Tutor académico: Daniel Vélez Serrano
Tutor de empresa: Guillermo Gorgas Ruiz (INE)

## Descripción

Pipeline reproducible para la estimación del salario real de los trabajadores 
asalariados de la EPA, integrando información de la Agencia Tributaria y la 
Seguridad Social mediante un modelo de variable latente MIMIC (Multiple 
Indicators, Multiple Causes). El pipeline se organiza en tres etapas:

1. **Estimación** del modelo MIMIC y cálculo de las puntuaciones factoriales 
   de Bartlett.
2. **Depuración** de observaciones atípicas mediante un residuo estudentizado 
   heterocedástico con control de la tasa de falsos descubrimientos (FDR).
3. **Imputación** del salario para observaciones sin información administrativa 
   suficiente, mediante regresión robusta de Huber y remuestreo no paramétrico 
   de residuos por estratos.

## Estructura del repositorio

| Archivo | Contenido |
|---|---|
| `config.R` | Configuración centralizada: semilla, parámetros de Huber, nivel FDR, etc. |
| `funciones_modelo.R` | Funciones auxiliares para la estimación del modelo MIMIC. |
| `funciones_depuracion.R` | Funciones auxiliares para la detección de atípicos. |
| `estimacion.R` | Script de la etapa 1: ajuste del modelo MIMIC y puntuaciones de Bartlett. |
| `depuracion.R` | Script de la etapa 2: detección de atípicos. |
| `imputacion.R` | Script de la etapa 3: imputación del salario. |

## Ejecución

Los scripts se ejecutan en orden: `estimacion.R` → `depuracion.R` → `imputacion.R`, 
cargando previamente `config.R`.

## Nota sobre los datos

Este repositorio contiene únicamente el código del pipeline. Los datos 
administrativos empleados (EPA, Agencia Tributaria, Seguridad Social) son 
confidenciales y no se incluyen, al proceder de registros del Instituto 
Nacional de Estadística.
