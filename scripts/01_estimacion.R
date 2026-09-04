# =============================================================================
# SCRIPT_PRODUCCION.R — GENERACIÓN DE PREDICCIONES + DIAGNÓSTICO DEL MODELO
# =============================================================================
#
# Ejecuta el pipeline completo de producción: carga, preprocesamiento,
# ajuste de los 2 modelos MIMIC, diagnóstico completo del modelo (cargas,
# parámetros estructurales, auditoría de residuos) y generación de
# predicciones. Todo sin usar STP3 en ningún momento.
#
# Modelos ajustados:
#   ajuste_mim : M1 MIMIC      (medición + causas estructurales)
#   ajuste_ind : M1_ind MIMIC  (medición + causas + causas en indicadores SS)
#
# Control de reajuste (útil si los datos no han cambiado):
#   reajustar_modelo <- TRUE   ajusta de nuevo y sobrescribe el guardado
#   reajustar_modelo <- FALSE  carga el ajuste previo si existe (por defecto)
#
# Estructura:
#   1. Carga y validación
#   2. Diagnóstico de pérdida de observaciones (sobre el dataset original)
#   3. Preprocesamiento
#   4. Ajuste de los 2 modelos (o carga del ajuste guardado)
#   5. Parámetros de medición (cargas, varianzas, Heywood, índices de ajuste)
#   6. Parámetros estructurales (efecto de cada causa)
#   7. Auditoría: residuos e índices de modificación (con corrección por
#      comparaciones múltiples, dado el elevado número de contrastes
#      simultáneos que implica revisar todos los índices de modificación)
#   8. Predicción (estimador de Bartlett, con verificación de missings y
#      corrección de sigma en la retransformación logarítmica)
#   9. Gráficos de diagnóstico sin STP3 (G3-G9)
#  10. Checks manuales y guardado del dataset final
#
# Salidas:
#   outputs/dataset_predicciones_<fecha>.rds / .csv
#   outputs/modelos_ajustados.rds
#   outputs/metadata/metadata_<fecha>.rds
#   outputs/graficos/...
#   logs/log_produccion_<fecha>.txt
# =============================================================================
library(tidyverse)
library(readxl)
library(lavaan)
library(here)


source(here("config.R"))
source(here("funciones_modelo.R"))

crear_carpetas_salida()

fecha_ejecucion    <- Sys.time()
tag_fecha          <- format(fecha_ejecucion, "%Y%m%d_%H%M")
sink(file.path(ruta_log, paste0("log_produccion_", tag_fecha, ".txt")),
     split = TRUE, append = FALSE)

cat("=============================================================\n")
cat(" PIPELINE DE PRODUCCION — SALARIO REAL (MIMIC)\n")
cat(" Version:", version_pipeline, "\n")
cat(" Fecha:  ", format(fecha_ejecucion), "\n")
cat("=============================================================\n\n")

cat("INDICE\n")
cat("  1. Carga y validación   — lee el Excel y comprueba columnas esperadas\n")
cat("  2. Diagnóstico pérdida  — por qué se descarta cada obs. del original\n")
cat("  3. Preprocesamiento     — indicadores log, mensualización, causas\n")
cat("  4. Ajuste de modelos    — M1 MIMIC y M1_ind (o carga del ajuste guardado)\n")
cat("  5. Parámetros medición  — cargas λ, varianzas θ, Heywood, RMSEA/CFI/SRMR\n")
cat("  6. Params. estructurales— efecto de cada causa sobre el salario real\n")
cat("  7. Auditoría            — residuos e índices de modificación\n")
cat("                            (corrección por comparaciones múltiples)\n")
cat("  8. Predicción           — Bartlett, verificación missings, corrección sigma\n")
cat("  9. Gráficos (sin STP3) — G3-G9\n")
cat(" 10. Checks manuales y dataset final\n")
cat("\n")


# * 1. CARGA Y VALIDACIÓN ----

cat("── 1. CARGA Y VALIDACION ─────────────────────────────────\n")

if (!file.exists(ruta_entrada)) {
  stop("No se encuentra el fichero de entrada en: ", ruta_entrada)
}

data_salarios <- read_excel(ruta_entrada)
cat("Fichero leido:", ruta_entrada, "\n")
cat("N filas:   ", nrow(data_salarios), "\n")
cat("N columnas:", ncol(data_salarios), "\n")

validar_columnas(data_salarios, columnas_esperadas)
cat("Validacion de columnas: OK\n\n")


# * 2. DIAGNÓSTICO DE PÉRDIDA DE OBSERVACIONES ----
# Calcula, sobre el dataset original completo, cuántas observaciones se
# pierden en cada etapa del preprocesamiento y por qué. Esto permite
# conocer el alcance real del pipeline (25.257 sobre 27.363) y es
# indispensable para diseñar la estrategia de imputación del salario para
# las observaciones que no entran en el MIMIC.
#
# Grupos de pérdida (mutuamente excluyentes, por orden de prioridad):
#   A: Sin ninguna fuente salarial bruta (AT=0, SS=0, MDR=0)
#      → No hay información administrativa con la que estimar. Solo
#        imputación por perfil demográfico/laboral.
#   B: Tienen fuente bruta pero pierden todos los ln_* tras aplicar
#      topes y mensualización:
#        B1: Sin FACTOR_final (ni FACTOR_SS ni FACTOR_EPA válido)
#            → No se puede mensualizar AT ni SS.
#        B2: Solo SS anual y está topada
#            → Única fuente eliminada por el criterio de topes.
#        B3: SS + MDR ambas topadas simultáneamente (sin AT disponible)
#            → Las dos fuentes de SS caen; AT no está presente.
#        B4: Otros (AT + otras fuentes todas topadas o con ln < 0)
#   C: Tienen algún ln_* válido pero HORASH = 9900 (código "no procede"
#      en la EPA) → HORASH_real queda NA y Horas_Z no se puede calcular.
#      Estas observaciones SÍ tienen fuente salarial válida; pendiente de
#      decisión con el tutor sobre cómo recuperarlas (imputar HORASH por
#      mediana de grupo Ocupación×Jornada, o MIMIC alternativo sin Horas_Z).

cat("── 2. DIAGNOSTICO DE PERDIDA DE OBSERVACIONES ────────────\n")

diagnostico_perdida(data_salarios)

cat("\n")


# * 3. PREPROCESAMIENTO ----

cat("── 3. PREPROCESAMIENTO ────────────────────────────────────\n")

datos          <- preprocesar_datos(data_salarios, contratos_parciales)
data_modelo    <- datos$base
data_modelo_ss <- datos$ss

cat("N modelo base:  ", nrow(data_modelo),    "\n")
cat("N modelo ind.SS:", nrow(data_modelo_ss), "\n\n")

if (nrow(data_modelo_ss) < 100) {
  stop("N tras preprocesamiento demasiado bajo (<100). ",
       "Revisar el fichero de entrada.")
}


# * 3B. IMPUTACIÓN DE HORASH Y RECUPERACIÓN DE OBSERVACIONES ----
# Tras el preprocesamiento, una parte de las observaciones se pierde por
# tener HORASH = 9900 ("no procede" en la EPA) pese a disponer de fuente
# salarial válida: como Horas_Z es una CAUSA del modelo y lavaan
# (fixed.x = TRUE) elimina por listwise cualquier caso con causa faltante,
# esas observaciones no entran en el MIMIC.
#
# El orden es deliberado: PRIMERO se ejecuta el preprocesamiento estándar
# (que deja constancia de cuántas observaciones se pierden por HORASH),
# y DESPUÉS se imputa HORASH por la mediana de su grupo Ocupación ×
# Jornada y se reincorporan esas observaciones. De este modo:
#   (a) queda documentada la pérdida real por HORASH en el flujo normal;
#   (b) las observaciones recuperadas entran en el MIMIC desde el inicio,
#       y por tanto pasan también por la posterior detección de atípicos
#       (depuracion.R) — evitando que un posible outlier entre esas
#       observaciones se cuele sin filtrar.
#
# La imputación de HORASH es determinista (mediana de grupo): más
# transparente y auditable que un método estocástico para una variable
# auxiliar de baja dispersión intragrupo.

cat("── 3B. IMPUTACIÓN DE HORASH Y RECUPERACIÓN ────────────────\n")

n_antes <- nrow(data_modelo_ss)

# Imputar HORASH por mediana de grupo (función compartida) y repreprocesar,
# de modo que las observaciones antes descartadas por HORASH = 9900 se
# reincorporen al MIMIC.
data_salarios_imp <- imputar_horash(data_salarios)

datos          <- preprocesar_datos(data_salarios_imp, contratos_parciales)
data_modelo    <- datos$base
data_modelo_ss <- datos$ss

cat("N modelo ind.SS antes de imputar HORASH:", n_antes, "\n")
cat("N modelo ind.SS tras imputar HORASH:    ", nrow(data_modelo_ss),
    sprintf("(+%d recuperadas)\n\n", nrow(data_modelo_ss) - n_antes))


# * 4. AJUSTE DE LOS 2 MODELOS (o carga del ajuste guardado) ----

cat("── 4. AJUSTE DE LOS MODELOS ───────────────────────────────\n")

if (!reajustar_modelo && file.exists(ruta_modelos_ajust)) {
  
  cat("Cargando modelos ya ajustados desde:\n  ", ruta_modelos_ajust, "\n",
      "(usar reajustar_modelo <- TRUE para forzar un nuevo ajuste)\n\n",
      sep = "")
  
  modelos_guardados <- readRDS(ruta_modelos_ajust)
  ajuste_mim        <- modelos_guardados$ajuste_mim
  ajuste_ind        <- modelos_guardados$ajuste_ind
  
} else {
  
  cat("Ajustando los 2 modelos (puede tardar varios minutos)...\n\n")
  
  ajuste_mim <- ajustar_mimic(data_modelo,    medida_M1, causas_estructurales)
  ajuste_ind <- ajustar_mimic(data_modelo_ss, medida_M1, causas_estructurales,
                              causas_ind_ss)
  
  saveRDS(list(ajuste_mim   = ajuste_mim,
               ajuste_ind   = ajuste_ind,
               fecha_ajuste = Sys.time()),
          ruta_modelos_ajust)
  cat("Modelos ajustados y guardados en:\n  ", ruta_modelos_ajust, "\n\n",
      sep = "")
}

# Se calculan una sola vez y se reutilizan en el resto del script
# (cat() de esta sección, sección 10 y metadata) sin volver a llamar a
# fitMeasures(), que con índices robustos es una operación costosa sobre
# un ajuste FIML con ~25.000 observaciones.

rmsea_mim <- fitMeasures(ajuste_mim, "rmsea.robust")
cfi_mim   <- fitMeasures(ajuste_mim, "cfi.robust")
rmsea_ind <- fitMeasures(ajuste_ind, "rmsea.robust")
cfi_ind   <- fitMeasures(ajuste_ind, "cfi.robust")

cat("M1 MIMIC   — RMSEA rob:", round(rmsea_mim, 4),
    "| CFI rob:", round(cfi_mim, 4), "\n")
cat("M1_ind     — RMSEA rob:", round(rmsea_ind, 4),
    "| CFI rob:", round(cfi_ind, 4), "\n\n")


# * 5. PARÁMETROS DE MEDICIÓN ----

cat("── 5. PARAMETROS DE MEDICION ──────────────────────────────\n\n")

imprimir_medicion(ajuste_mim, "M1 MIMIC (con causas)")
imprimir_medicion(ajuste_ind, "M1_ind MIMIC (con causas + ind.SS)")


# * 6. PARÁMETROS ESTRUCTURALES ----

cat("── 6. PARAMETROS ESTRUCTURALES ────────────────────────────\n\n")

imprimir_estructural(ajuste_mim, "M1 MIMIC — causas del factor latente")
imprimir_estructural(ajuste_ind, "M1_ind MIMIC — causas del factor latente")
imprimir_causas_ind_ss(ajuste_ind)


# * 7. AUDITORÍA: RESIDUOS E ÍNDICES DE MODIFICACIÓN ----
# imprimir_auditoria() ya reporta los índices de modificación (mi) junto
# con su epc. Aquí se añade la corrección por comparaciones múltiples:
# al revisar simultáneamente decenas de índices de modificación (cada uno
# es, en esencia, un contraste de razón de verosimilitudes con 1 grado de
# libertad), evaluarlos uno a uno con el umbral habitual de mi>3.84
# (equivalente a p<0.05) infla la tasa de falsos positivos. Se reporta
# también el umbral corregido por Bonferroni para que la interpretación
# de "qué índice es realmente relevante" sea más conservadora.

cat("── 7. RESIDUOS E INDICES DE MODIFICACION ──────────────────\n\n")

imprimir_auditoria(ajuste_mim, "M1 MIMIC")
imprimir_auditoria(ajuste_ind, "M1_ind MIMIC")

reportar_correccion_multiple <- function(ajuste, nombre) {
  mi_todos <- modindices(ajuste, sort. = TRUE, minimum.value = 0)
  n_tests  <- nrow(mi_todos)
  alpha_bonferroni <- 0.05 / n_tests
  mi_critico_bonferroni <- qchisq(1 - alpha_bonferroni, df = 1)
  
  n_sig_individual  <- sum(mi_todos$mi > qchisq(0.95, df = 1))
  n_sig_bonferroni   <- sum(mi_todos$mi > mi_critico_bonferroni)
  
  cat(">>> ", nombre, " — correccion por comparaciones multiples <<<\n", sep = "")
  cat("  N indices de modificacion evaluados:", n_tests, "\n")
  cat("  mi critico SIN corregir (p<0.05):     ", round(qchisq(0.95, df = 1), 3),
      " -> ", n_sig_individual, " indices superan el umbral\n", sep = "")
  cat("  mi critico CON Bonferroni (alpha/N):  ", round(mi_critico_bonferroni, 3),
      " -> ", n_sig_bonferroni, " indices superan el umbral\n\n", sep = "")
}

reportar_correccion_multiple(ajuste_mim, "M1 MIMIC")
reportar_correccion_multiple(ajuste_ind, "M1_ind MIMIC")


# * 8. PREDICCIÓN (ESTIMADOR DE BARTLETT) ----
# Dos verificaciones añadidas respecto a una predicción Bartlett simple:
#
# (a) Comprobación de missings: lavPredict() con FIML no debería generar
#     NA salvo para observaciones que no tengan NINGÚN indicador
#     disponible. Se comprueba explícitamente que no se estén imputando
#     o perdiendo más observaciones de las esperadas.
#
# (b) Corrección de sigma en la retransformación logarítmica: exp(mu) es
#     la MEDIANA de la distribución lognormal, no la media. Si se quiere
#     reportar la media esperada del salario, hay que aplicar
#     E[Y] = exp(mu + sigma^2/2), donde sigma^2 es la varianza del error
#     de la puntuación factorial de Bartlett. Se calculan y reportan
#     ambas versiones (mediana y media corregida) para que se pueda
#     decidir cuál usar según el propósito del informe.

cat("── 8. PREDICCION (BARTLETT) ──────────────────────────────\n")

# (a) Verificación de missings antes de predecir
verificar_missings_bartlett <- function(data, ajuste, nombre) {
  n_total <- nrow(data)
  n_sin_ningun_indicador <- data %>%
    filter(is.na(ln_AT_mens), is.na(ln_SS_mens), is.na(ln_MDR)) %>%
    nrow()
  
  pred_temp <- predecir_bartlett(ajuste, data)
  n_pred_na <- sum(is.na(pred_temp))
  
  cat(">>> ", nombre, " — verificacion de missings en Bartlett <<<\n", sep = "")
  cat("  N total:                                  ", n_total, "\n")
  cat("  N sin NINGUN indicador disponible:         ", n_sin_ningun_indicador, "\n")
  cat("  N con prediccion NA tras lavPredict:       ", n_pred_na, "\n")
  
  if (n_pred_na != n_sin_ningun_indicador) {
    cat("  *** ATENCION: el numero de NA en la prediccion no coincide con\n")
    cat("      el numero de observaciones sin ningun indicador. Revisar si\n")
    cat("      lavPredict esta imputando de mas o perdiendo observaciones\n")
    cat("      que si tenian al menos un indicador disponible. ***\n")
  } else {
    cat("  OK: la prediccion solo es NA cuando no hay ningun indicador.\n")
    cat("      No se estan introduciendo missings nuevos ni imputaciones\n")
    cat("      de mas alla de lo que permite FIML.\n")
  }
  cat("\n")
  invisible(pred_temp)
}

verificar_missings_bartlett(data_modelo,    ajuste_mim, "M1 MIMIC")
verificar_missings_bartlett(data_modelo_ss, ajuste_ind, "M1_ind MIMIC")

# (b) Varianza del error de Bartlett, para la corrección de sigma
calcular_var_error_bartlett <- function(ajuste) {
  params <- parameterEstimates(ajuste)
  
  cargas <- params %>%
    filter(op == "=~", lhs == "Salario_Real") %>%
    select(rhs, est) %>%
    deframe()
  
  theta <- params %>%
    filter(op == "~~", lhs == rhs,
           lhs %in% names(cargas)) %>%
    select(lhs, est) %>%
    deframe()
  
  # Var(error Bartlett) = 1 / sum(lambda_i^2 / theta_i)
  1 / sum(cargas^2 / theta[names(cargas)])
}

var_error_mim <- calcular_var_error_bartlett(ajuste_mim)
var_error_ind <- calcular_var_error_bartlett(ajuste_ind)

cat("Varianza del error de Bartlett — M1 MIMIC:", round(var_error_mim, 4), "\n")
cat("Varianza del error de Bartlett — M1_ind:  ", round(var_error_ind, 4), "\n\n")

# Predicción: se calculan ambas versiones (mediana sin corregir, y media
# corregida con sigma^2/2) para dejar constancia de la diferencia
data_modelo <- data_modelo %>%
  mutate(
    Score_mim           = predecir_bartlett(ajuste_mim, data_modelo, devolver_score = TRUE),
    Pred_M1_mim_mediana = exp(Score_mim),
    Pred_M1_mim_B       = exp(Score_mim + var_error_mim / 2)
  )

data_modelo_ss <- data_modelo_ss %>%
  mutate(
    Score_ind              = predecir_bartlett(ajuste_ind, data_modelo_ss, devolver_score = TRUE),
    Score_mim_ss           = predecir_bartlett(ajuste_mim, data_modelo_ss, devolver_score = TRUE),
    Pred_M1_ind_mediana    = exp(Score_ind),
    Pred_M1_ind_B          = exp(Score_ind    + var_error_ind / 2),
    Pred_M1_mim_mediana_ss = exp(Score_mim_ss),
    Pred_M1_mim_B          = exp(Score_mim_ss + var_error_mim / 2)
  )

cat("M1 MIMIC   — Media (mediana lognormal):", round(mean(data_modelo$Pred_M1_mim_mediana, na.rm=TRUE), 0), "\n")
cat("M1 MIMIC   — Media (corregida sigma):  ", round(mean(data_modelo$Pred_M1_mim_B,       na.rm=TRUE), 0), "\n")
cat("M1_ind     — Media (mediana lognormal):", round(mean(data_modelo_ss$Pred_M1_ind_mediana, na.rm=TRUE), 0), "\n")
cat("M1_ind     — Media (corregida sigma):  ", round(mean(data_modelo_ss$Pred_M1_ind_B,       na.rm=TRUE), 0),
    "| SD:", round(sd(data_modelo_ss$Pred_M1_ind_B, na.rm=TRUE), 0), "\n\n")


# * 9. GRÁFICOS SIN STP3 ----
# G3: distribución de predicciones
# G4: predicción por ocupación (boxplot)
# G5: predicción por educación y jornada
# G6: M1_ind vs AT mensual observado
# G7: M1_ind vs SS mensual observado
# G8: M1_ind vs MDR mensual observado
# G9: diferencia entre M1_ind y M1 MIMIC

cat("── 9. GRAFICOS (sin STP3) ───────────────────────────────\n")
t_inicio_graficos <- Sys.time()

p99_ind <- quantile(data_modelo_ss$Pred_M1_ind_B, 0.99, na.rm = TRUE)

# G3: distribución de predicciones con líneas de referencia
g3 <- data_modelo_ss %>%
  filter(!is.na(Pred_M1_ind_B), Pred_M1_ind_B < p99_ind) %>%
  ggplot(aes(x = Pred_M1_ind_B)) +
  geom_histogram(bins = 80, fill = "#534AB7", color = "white", alpha = 0.8) +
  geom_vline(xintercept = c(1323, 4720), color = "#D85A30",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 1323, y = Inf, label = "SMI (1.323 euros)",
           vjust = 2, hjust = -0.1, color = "#D85A30", size = 3.5) +
  annotate("text", x = 4720, y = Inf, label = "Tope MDR (4.720 euros)",
           vjust = 2, hjust = -0.1, color = "#D85A30", size = 3.5) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Prediccion Bartlett (euros/mes)", y = "Frecuencia",
       title = "Distribucion de las predicciones del salario real") +
  theme_minimal()
guardar_grafico(g3, ruta_grafico, "G3_distribucion_predicciones")

# G4: boxplot de predicciones por ocupación, ordenado por mediana
g4 <- data_modelo_ss %>%
  filter(!is.na(Pred_M1_ind_B), !is.na(Ocupacion), Pred_M1_ind_B < p99_ind) %>%
  mutate(Ocupacion = reorder(Ocupacion, Pred_M1_ind_B, median)) %>%
  ggplot(aes(x = Ocupacion, y = Pred_M1_ind_B, fill = Ocupacion)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_viridis_d(guide = "none") +
  labs(x = NULL, y = "Prediccion Bartlett (euros/mes)",
       title = "Distribucion de predicciones por ocupacion") +
  coord_flip() +
  theme_minimal()
guardar_grafico(g4, ruta_grafico, "G4_prediccion_por_ocupacion")

# G5: predicciones por nivel educativo y jornada
g5 <- data_modelo_ss %>%
  filter(!is.na(Pred_M1_ind_B), Pred_M1_ind_B < p99_ind) %>%
  mutate(
    Jornada_lab = ifelse(Jornada_Completa == 1, "Completa", "Parcial"),
    Educ_lab    = case_when(
      Educ_Superior   == 1 ~ "Superior",
      Educ_FP_Bach    == 1 ~ "FP/Bach",
      Educ_Secundaria == 1 ~ "Secundaria",
      TRUE                 ~ "Basica"
    ),
    Educ_lab = factor(Educ_lab,
                      levels = c("Basica", "Secundaria", "FP/Bach", "Superior"))
  ) %>%
  ggplot(aes(x = Educ_lab, y = Pred_M1_ind_B, fill = Jornada_lab)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.4, outlier.alpha = 0.3) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("Completa" = "#1D9E75", "Parcial" = "#D85A30")) +
  labs(x = "Nivel educativo", y = "Prediccion Bartlett (euros/mes)",
       fill = "Jornada",
       title = "Predicciones por nivel educativo y tipo de jornada") +
  theme_minimal() +
  theme(legend.position = "bottom")
guardar_grafico(g5, ruta_grafico, "G5_prediccion_educacion_jornada")

g6 <- data_modelo_ss %>%
  filter(!is.na(ln_AT_mens), !is.na(Pred_M1_ind_B)) %>%
  mutate(AT_obs = exp(ln_AT_mens)) %>%
  filter(AT_obs < quantile(AT_obs, 0.99, na.rm = TRUE),
         Pred_M1_ind_B < p99_ind) %>%
  ggplot(aes(x = AT_obs, y = Pred_M1_ind_B)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#1D9E75") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_smooth(method = "lm", color = "#D85A30", se = FALSE, linewidth = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "AT mensual observado (euros/mes)", y = "M1_ind Bartlett (euros/mes)",
       title = "M1_ind vs AT mensual observado") +
  theme_minimal()
guardar_grafico(g6, ruta_grafico, "G6_scatter_vs_AT")

g7 <- data_modelo_ss %>%
  filter(!is.na(ln_SS_mens), !is.na(Pred_M1_ind_B)) %>%
  mutate(SS_obs = exp(ln_SS_mens)) %>%
  filter(SS_obs < quantile(SS_obs, 0.99, na.rm = TRUE),
         Pred_M1_ind_B < p99_ind) %>%
  ggplot(aes(x = SS_obs, y = Pred_M1_ind_B)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#534AB7") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_smooth(method = "lm", color = "#D85A30", se = FALSE, linewidth = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "SS mensual observado (euros/mes)", y = "M1_ind Bartlett (euros/mes)",
       title = "M1_ind vs SS mensual observado") +
  theme_minimal()
guardar_grafico(g7, ruta_grafico, "G7_scatter_vs_SS")

g8 <- data_modelo_ss %>%
  filter(!is.na(ln_MDR), !is.na(Pred_M1_ind_B)) %>%
  mutate(MDR_obs = exp(ln_MDR)) %>%
  filter(MDR_obs < quantile(MDR_obs, 0.99, na.rm = TRUE),
         Pred_M1_ind_B < p99_ind) %>%
  ggplot(aes(x = MDR_obs, y = Pred_M1_ind_B)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#BA7517") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_smooth(method = "lm", color = "#D85A30", se = FALSE, linewidth = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "MDR mensual observado (euros/mes)", y = "M1_ind Bartlett (euros/mes)",
       title = "M1_ind vs MDR mensual observado") +
  theme_minimal()
guardar_grafico(g8, ruta_grafico, "G8_scatter_vs_MDR")

# G9: diferencia entre M1_ind y M1 MIMIC — efecto de las causas en ind. SS
g9 <- data_modelo_ss %>%
  filter(!is.na(Pred_M1_mim_B), !is.na(Pred_M1_ind_B)) %>%
  mutate(
    Dif_modelos = Pred_M1_ind_B - Pred_M1_mim_B,
    Jornada_lab = ifelse(Jornada_Completa == 1, "Completa", "Parcial")
  ) %>%
  filter(abs(Dif_modelos) < quantile(abs(Dif_modelos), 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Dif_modelos, fill = Jornada_lab)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Completa" = "#1D9E75", "Parcial" = "#D85A30")) +
  labs(x = "M1_ind vs M1 MIMIC (euros/mes)", y = "Frecuencia",
       fill = "Jornada",
       title = "Diferencia entre M1_ind y M1 MIMIC por jornada",
       subtitle = "Efecto neto de anadir las causas en indicadores SS") +
  theme_minimal() +
  theme(legend.position = "bottom")
guardar_grafico(g9, ruta_grafico, "G9_diferencia_modelos")

t_fin_graficos <- Sys.time()
cat("\nTiempo en generar los 7 gráficos (G3-G9):",
    round(difftime(t_fin_graficos, t_inicio_graficos, units = "secs"), 1), "segundos\n")

cat("\n")


# * 10. CHECKS MANUALES Y DATASET FINAL ----
# Sustituye la comparación automática contra la ejecución anterior por un
# resumen de valores clave para que el usuario los revise y decida si la
# ejecución es válida antes de darla por buena. Pensado para un script
# que se ejecuta con poca frecuencia, donde una comparación automática
# contra "la ejecución anterior" (que puede ser de hace meses) aporta
# poco.

cat("── 10. CHECKS MANUALES Y DATASET FINAL ─────────────────────\n")
t_inicio_seccion10 <- Sys.time()

media_pred <- mean(data_modelo_ss$Pred_M1_ind_B, na.rm = TRUE)
sd_pred    <- sd(data_modelo_ss$Pred_M1_ind_B,   na.rm = TRUE)

cat("Revisar los siguientes valores antes de aceptar esta ejecucion:\n\n")

cat("  1. Tamano de muestra:\n")
cat("     N filas en el Excel de entrada:  ", nrow(data_salarios), "\n")
cat("     N tras preprocesamiento (ind.SS):", nrow(data_modelo_ss), "\n")
cat("     N perdidos (ver seccion 2 para detalle por causa):",
    nrow(data_salarios) - nrow(data_modelo_ss), "\n\n")

cat("  2. Ajuste del modelo:\n")
cat("     RMSEA robusto (M1_ind): ", round(rmsea_ind, 4),
    " (deseable < 0.06)\n")
cat("     CFI robusto   (M1_ind): ", round(cfi_ind, 4),
    " (deseable > 0.95)\n\n")

cat("  3. Magnitud de la prediccion:\n")
cat("     Media (corregida sigma):", round(media_pred, 0), "euros/mes\n")
cat("     SD:                     ", round(sd_pred, 0), "euros/mes\n\n")

cat("  4. Revision visual recomendada:\n")
cat("     - G3 (distribucion): ¿forma razonable, sin picos extraños?\n")
cat("     - G9 (diferencia entre modelos): ¿el efecto de las causas SS\n")
cat("       tiene el signo y magnitud esperados?\n\n")

dataset_final <- data_modelo_ss %>%
  select(
    idx,
    any_of(c("IDENCPERHOG", "NVIVI", "NPERS")),
    Pred_M1_ind_B
  )

saveRDS(dataset_final,
        file.path(ruta_output, paste0("dataset_predicciones_", tag_fecha, ".rds")))
write_csv(dataset_final,
          file.path(ruta_output, paste0("dataset_predicciones_", tag_fecha, ".csv")))

# La metadata se conserva como registro histórico de cada ejecución (útil
# para trazabilidad y para comparar manualmente si se desea en el futuro),
# pero ya no se compara automáticamente ni genera alertas.
metadata_actual <- list(
  version          = version_pipeline,
  fecha_ejecucion  = fecha_ejecucion,
  fichero_entrada  = ruta_entrada,
  n_filas_entrada  = nrow(data_salarios),
  n_modelo         = nrow(data_modelo_ss),
  rmsea_mim        = rmsea_mim,
  cfi_mim          = cfi_mim,
  rmsea_ind        = rmsea_ind,
  cfi_ind          = cfi_ind,
  media_prediccion = media_pred,
  sd_prediccion    = sd_pred
)
saveRDS(metadata_actual,
        file.path(ruta_metadata, paste0("metadata_", tag_fecha, ".rds")))

cat("Dataset final guardado: ", nrow(dataset_final), "filas\n")
cat("Metadata guardada (registro historico, sin comparacion automatica).\n\n")

t_fin_seccion10 <- Sys.time()
cat("Tiempo en la sección 10 (checks + guardado):",
    round(difftime(t_fin_seccion10, t_inicio_seccion10, units = "secs"), 1), "segundos\n\n")

cat("=============================================================\n")
cat(" PRODUCCION COMPLETADA\n")
cat("=============================================================\n")

sink()
