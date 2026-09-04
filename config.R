# =============================================================================
# CONFIG.R — CONFIGURACIÓN CENTRALIZADA DEL PIPELINE MIMIC
# =============================================================================
#
# Punto único de configuración: rutas, nombres de fichero y parámetros del
# modelo. Los scripts de producción y evaluación cargan este fichero al
# inicio (source("config.R")) en lugar de definir estos valores por su
# cuenta, de modo que un cambio aquí se propaga a todo el pipeline.
#
# Si cambia el nombre del fichero de entrada cada mes, o se mueve la
# ubicación de alguna carpeta, solo hay que tocar este fichero.
# =============================================================================


# ── Versión del pipeline ──────────────────────────────────────────────────
# Incrementar manualmente cuando se modifique la especificación del modelo
# (variables incluidas, estimador, criterios de filtrado). Queda registrada
# en los metadatos de cada ejecución para trazabilidad.
version_pipeline <- "1.0.0"

# ── Fichero de entrada ────────────────────────────────────────────────────
# En la raíz del proyecto (here()).
ruta_entrada <- "~/SGEML_PRATICAS/EPASAL/ASALARIADOS_ANONIMIZADOS.xlsx"

# ── Rutas de salida ───────────────────────────────────────────────────────
# outputs/ y logs/ viven directamente en la raíz del proyecto (EPASAL),
# como carpetas hermanas de pipeline/ y atipicos/.
ruta_output    <- here("outputs")
ruta_log       <- here("logs")
ruta_grafico   <- file.path(ruta_output, "graficos")
ruta_metadata  <- file.path(ruta_output, "metadata")

# Fichero donde se guardan los modelos MIMIC ya ajustados (ajuste_mim,
# ajuste_ind), para que script_evaluacion_stp3.R pueda cargarlos sin
# reajustar, y para que script_produccion.R no reajuste si ya existen
# y no se solicita explícitamente lo contrario (reajustar_modelo).
ruta_modelos_ajust <- file.path(ruta_output, "modelos_ajustados.rds")

crear_carpetas_salida <- function() {
  if(!(dir.exists(ruta_output) & dir.exists(ruta_log) 
       & dir.exists(ruta_grafico) & dir.exists(ruta_metadata))){
    dir.create(ruta_output,    recursive = TRUE)
    dir.create(ruta_log,       recursive = TRUE)
    dir.create(ruta_grafico,   recursive = TRUE)
    dir.create(ruta_metadata,  recursive = TRUE)
  }
}

# ── Columnas mínimas exigidas en el fichero de entrada ──────────────────────
# Si el fichero de origen cambia de formato (columna renombrada o eliminada),
# el pipeline debe fallar de forma explícita en lugar de producir resultados
# silenciosamente incorrectos.
columnas_esperadas <- c(
  "AT_SALARIO", "SALARIO_SS_ANO", "SALARIO_SS_MDR",
  "EDAD1", "SEXO1", "DCOM", "PARCO1",
  "TOPE_SS_ANO", "TOPE_SS_MDR",
  "FACTOR_SS", "FACTOR_EPA", "STP3",
  "HORASH", "NAC1", "NFORMA", "OCUP", "ACT",
  "NGRUP_MDR", "cntrato"
)

# ── Códigos de contrato a tiempo parcial ─────────────────────────────────────
contratos_parciales <- c("200", "209", "239", "250", "289",
                         "500", "501", "502", "506", "508",
                         "510", "518", "520")

# ── Especificación del modelo (sintaxis lavaan) ──────────────────────────────
medida_M1 <- '
  Salario_Real =~ 1*ln_AT_mens + ln_SS_mens + ln_MDR
'

causas_estructurales <- '
  Salario_Real ~ Jornada_Completa + Mujer_Dummy + Edad_Z + Edad_Z2 +
                 Exp_Z + Exp_Z2 + Extranjero + Horas_Z +
                 Educ_Secundaria + Educ_FP_Bach + Educ_Superior +
                 Ocup_Directivos + Ocup_TecSup + Ocup_TecAp +
                 Ocup_Admin + Ocup_Servicios + Ocup_Agricultura +
                 Ocup_Artesanos + Ocup_Operadores +
                 Sector_Agri + Sector_Ind + Sector_Const
'

causas_ind_ss <- '
  ln_SS_mens ~ Contrato_Parc + GC_alta + GC_baja
  ln_MDR     ~ Contrato_Parc + GC_alta + GC_baja
'

# ── Variables de control de ejecución ────────────────────────────────────
# Centralizadas aquí para que estimacion.R pueda ejecutarse de forma
# completamente independiente, sin pasar por ningún script orquestador
# externo. Si se quiere forzar un comportamiento distinto en una sesión
# concreta, basta con definir la variable correspondiente ANTES de hacer
# source(config.R) o de ejecutar estimacion.R — esa definición previa
# tiene prioridad porque el chequeo de abajo es if(!exists(...)).
#
#   reajustar_modelo <- TRUE   ajusta de nuevo el MIMIC y sobrescribe
#                               outputs/modelos_ajustados.rds
#   reajustar_modelo <- FALSE  (defecto) carga el ajuste guardado si
#                               existe; si no existe, ajusta igualmente
if (!exists("reajustar_modelo")) reajustar_modelo <- FALSE


# =============================================================================
# CONFIGURACIÓN DE LA FASE DE DEPURACIÓN (depuracion.R)
# =============================================================================
# Esta sección añade lo necesario para la detección de atípicos sobre la
# predicción ya generada por estimacion.R. No afecta a la fase de
# estimación en sí.

# ── Causas estructurales como VECTOR de nombres ──────────────────────────
# A diferencia de causas_estructurales (arriba), que es la sintaxis lavaan
# completa usada por el MIMIC, aquí se necesita un vector simple de
# nombres de columnas para ajustar la regresión robusta y construir la
# fórmula del IQR/Eurostat. Debe mantenerse coherente con las variables
# de causas_estructurales y causas_ind_ss.
causas_base <- c(
  "Edad_Z", "Edad_Z2", "Exp_Z", "Exp_Z2", "Horas_Z",
  "Jornada_Completa", "Mujer_Dummy", "Extranjero",
  "Educ_Superior", "Educ_FP_Bach", "Educ_Secundaria",
  "Ocup_Directivos", "Ocup_TecSup", "Ocup_TecAp",
  "Ocup_Admin", "Ocup_Servicios", "Ocup_Artesanos", "Ocup_Operadores",
  "Sector_Agri", "Sector_Ind", "Sector_Const",
  "GC_alta", "GC_baja", "Contrato_Parc"
)

# ── Parámetros de la regresión robusta (Huber) ────────────────────────────
# k=1.345 da 95% de eficiencia bajo normalidad (Huber, 1981). Aquí se usa
# para AJUSTAR el modelo robusto (rlm), no como criterio de detección —
# la detección final usa FDR sobre el p-valor real del residuo
# studentizado z (ver calcular_residuo_heterocedastico en
# funciones_depuracion.R).
huber_k                 <- 1.345
huber_k_extremo         <- 2.5
huber_alpha_bonferroni  <- 0.05   # alpha de partida para Bonferroni y FDR

# ── Factores del criterio IQR (Tukey / Eurostat) ──────────────────────────
iqr_factor_normal  <- 1.5
iqr_factor_extremo <- 3.0

# ── Rutas de salida de la fase de depuración ──────────────────────────────
ruta_out_depura  <- file.path(ruta_output, "depuracion")

crear_carpetas_depuracion <- function() {
  dir.create(ruta_out_depura,  showWarnings = FALSE, recursive = TRUE)
}

# ── CRITERIO CONFIGURABLE DE ATÍPICO FINAL ────────────────────────────────
# Define en un único punto qué combinación de métodos determina
# Flag_atipico_final, para poder cambiarlo sin tocar la lógica de
# depuracion.R. Opciones disponibles:
#
#   "solo_robusta"        — Flag_atipico_final = Flag_Robusta_bin
#                            (criterio actual: la robusta es necesaria
#                            y suficiente; el IQR solo aporta el nivel
#                            de severidad/refuerzo, no decide por sí solo)
#   "robusta_o_iqr"        — Flag_atipico_final = 1 si la robusta O el
#                            IQR detectan (unión de ambos conjuntos)
#   "robusta_y_iqr"        — Flag_atipico_final = 1 solo si AMBOS
#                            coinciden (intersección, más conservador)
#
# Para cambiar el criterio de detección de atípicos en el futuro, basta
# con cambiar el valor de esta variable — no es necesario modificar
# depuracion.R.
if (!exists("criterio_atipico_final")) criterio_atipico_final <- "solo_robusta"


# =============================================================================
# CONFIGURACIÓN DE LA FASE DE IMPUTACIÓN (imputacion.R)
# =============================================================================
# Esta sección añade lo necesario para imputar el salario de las
# observaciones atípicas y de las que no obtienen predicción MIMIC, de modo
# que el dataset final tenga un valor de salario para las 27.363
# observaciones del dataset original. No afecta a estimación ni depuración.

# ── Semilla para reproducibilidad del ruido no paramétrico ─────────────────
# El muestreo de residuos (donantes) es aleatorio; fijar la semilla
# garantiza que la imputación sea reproducible entre ejecuciones.
semilla_imputacion <- 20260629

# ── Nº de tramos salariales del estrato de varianza ───────────────────────
# El estrato del ruido es Ocupación × Jornada_Completa × tramo de predicción.
# 4 = cuartiles de la predicción robusta.
n_tramos_estrato <- 4

# ── Nº máximo de intentos del muestreo por rechazo ────────────────────────
# En cada intento se prueba pred + residuo aleatorio; si cae fuera del
# rango admisible se reintenta. Tras agotar los intentos, se usa el
# fallback determinista (predicción acotada al rango).
max_intentos_rechazo <- 20

# ── Rutas de salida de la fase de imputación ──────────────────────────────
ruta_out_imputa <- file.path(ruta_output, "imputacion")

crear_carpetas_imputacion <- function() {
  dir.create(ruta_out_imputa, showWarnings = FALSE, recursive = TRUE)
}

# ── Parámetros para el suelo del salario imputado ─────────────────────────
# El SMI mensual se usa como suelo del salario imputado: un asalariado no
# debería aparecer con un salario imputado por debajo del mínimo legal que
# le corresponde según su jornada. NO se usan las bases de cotización
# (mínima por grupo ni máxima BMAX) como límites del salario, porque son
# topes de COTIZACIÓN, no de salario real: el salario real puede estar
# legítimamente por encima de la base máxima (p.ej. un directivo cotiza por
# el tope pero cobra más) o, en la parte de cotización mínima, no es un
# suelo del salario en sí. Por eso solo se aplica un suelo (SMI) y NO se
# aplica ningún techo.
smi_mensual_2024 <- 1134.00

# ── Límite inferior admisible del salario imputado (euros/mes) ─────────────
# Suelo SOLO para jornada completa: un asalariado a jornada completa no
# debería aparecer con un salario imputado por debajo del SMI. Para jornada
# parcial NO se aplica suelo: el SMI parcial no existe como tal (dependería
# de las horas exactas trabajadas, y no se dispone de una variable de horas
# 100% fiable), por lo que prorratear el SMI sería una aproximación gruesa;
# es preferible no imponer suelo en ese caso. `fila` es una fila del
# data.frame con la columna Jornada_Completa.
smi_prorrateado <- function(fila) {
  jc <- if (!is.null(fila$Jornada_Completa)) fila$Jornada_Completa else NA
  if (!is.na(jc) && jc == 1) {
    smi_mensual_2024        # jornada completa → suelo = SMI
  } else {
    0                       # jornada parcial (o desconocida) → sin suelo
  }
}

# ── Sin techo ──────────────────────────────────────────────────────────────
# No se impone límite superior al salario imputado: no existe un tope real
# de salario (BMAX es tope de cotización, no de salario). El ruido se toma
# de residuos reales de gente del mismo estrato, por lo que los valores
# altos imputados son plausibles por construcción.
max_admisible <- function(fila) {
  Inf
}