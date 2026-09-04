# =============================================================================
# IMPUTACION.R — IMPUTACIÓN DEL SALARIO REAL
# =============================================================================
#
# Tercer y último script principal del pipeline (tras estimacion.R y
# depuracion.R). Genera un valor de salario para TODAS las observaciones
# del dataset original (27.363), imputando:
#
#   (A) Los atípicos detectados por depuracion.R (Flag_atipico_final = 1),
#       cuyo valor observado es poco fiable.
#   (B) Las observaciones sin predicción MIMIC:
#         - Las que tienen HORASH = 9900: se imputa HORASH por mediana de
#           grupo y se recuperan para el MIMIC (dejan de necesitar
#           imputación de salario al obtener predicción).
#         - Las que no tienen ningún indicador administrativo: no tienen
#           predicción MIMIC posible; se imputan con la regresión robusta
#           (Mincer sobre perfil) + ruido no paramétrico.
#
# MÉTODO DE IMPUTACIÓN DE SALARIO (acordado con el tutor):
#   Regresión robusta única (Huber) de la predicción sobre el perfil
#   sociodemográfico, + ruido no paramétrico por remuestreo de residuos
#   reales de casos limpios (ni atípicos ni topados), estratificados por
#   Ocupación × Jornada_Completa × cuartil de predicción. Muestreo por
#   rechazo para garantizar que el valor imputado sea admisible (≥ suelo
#   de cotización/SMI prorrateado y ≤ base máxima de cotización).
#
# PRINCIPIOS (indicados por el tutor):
#   - Se parte SIEMPRE del dataset original (27.363). Los filtros se
#     renombran, nunca se pisa el original. `idx` es el identificador
#     estable en todo momento.
#   - Toda observación imputada queda marcada con una flag (Imputado) y
#     el método usado (Metodo_imputacion).
#
# Estructura:
#   1. Carga (dataset original + modelo + atípicos)
#   2. Preprocesamiento con HORASH imputada (predice con el modelo ya
#      ajustado en estimacion.R; NO reajusta ni imputa horas de nuevo)
#   3. Regresión robusta (Mincer) sobre perfil + estratos de varianza
#   4. Imputación por ruido no paramétrico (atípicos + sin indicador)
#   5. Ensamblado del dataset final (27.363, con flags)
#   6. Checks de validez
#   7. Summary y gráficos
#   8. Guardado
# =============================================================================

library(tidyverse)
library(readxl)
library(lavaan)
library(MASS)
library(here)

select <- dplyr::select

source(here("config.R"))
source(here("funciones_modelo.R"))
source(here("funciones_depuracion.R"))

crear_carpetas_salida()
crear_carpetas_depuracion()
crear_carpetas_imputacion()

fecha_ejecucion <- Sys.time()
tag_fecha       <- format(fecha_ejecucion, "%Y%m%d_%H%M")
sink(file.path(ruta_log, paste0("log_imputacion_", tag_fecha, ".txt")),
     split = TRUE, append = FALSE)

set.seed(semilla_imputacion)

cat("=============================================================\n")
cat(" IMPUTACIÓN — SALARIO REAL PARA EL DATASET COMPLETO\n")
cat(" Versión:", version_pipeline, "| Fecha:", format(fecha_ejecucion), "\n")
cat(" Semilla:", semilla_imputacion, "\n")
cat("=============================================================\n\n")


# * 1. CARGA (DATASET ORIGINAL + MODELO + ATÍPICOS) ----
# Se carga el dataset original completo (27.363), que se conserva intacto
# bajo el nombre `data_original` durante todo el script. Cualquier filtrado
# posterior se hace sobre copias renombradas, nunca sobre este objeto.
# `idx` (row_number sobre el original) es el identificador que permite
# reensamblar todo al final. Se cargan también el modelo MIMIC ajustado
# por estimacion.R y el dataset de atípicos generado por depuracion.R.

cat("── 1. CARGA ──────────────────────────────────────────────\n")

if (!file.exists(ruta_entrada))       stop("No se encuentra el fichero de entrada: ", ruta_entrada)
if (!file.exists(ruta_modelos_ajust)) stop("No existe modelos_ajustados.rds. Ejecutar estimacion.R.")

data_original <- read_excel(ruta_entrada) %>% mutate(idx = row_number())
validar_columnas(data_original, columnas_esperadas)
cat("Dataset original cargado:", nrow(data_original), "observaciones (intacto).\n")

modelos_guardados <- readRDS(ruta_modelos_ajust)
ajuste_ind        <- modelos_guardados$ajuste_ind

# Dataset de atípicos más reciente en la carpeta de depuración
ficheros_atip <- list.files(ruta_out_depura, pattern = "^dataset_atipicos_.*\\.rds$",
                            full.names = TRUE)
if (length(ficheros_atip) == 0) stop("No hay dataset de atípicos. Ejecutar depuracion.R.")
ruta_atip_reciente <- ficheros_atip[which.max(file.mtime(ficheros_atip))]
dataset_atipicos   <- readRDS(ruta_atip_reciente)
cat("Atípicos cargados desde:", basename(ruta_atip_reciente),
    "|", sum(dataset_atipicos$Flag_atipico_final, na.rm = TRUE), "atípicos\n\n")


# * 2. PREPROCESAMIENTO CON HORASH IMPUTADA ----
# La imputación de HORASH y el ajuste del MIMIC con las observaciones
# recuperadas ya se realizan en estimacion.R (sección 3B). Aquí solo se
# reconstruye el MISMO dataset preprocesado —imputando HORASH con idéntico
# criterio (mediana de grupo Ocupación × Jornada)— para poder generar la
# predicción con el modelo ya ajustado y guardado. NO se reajusta el
# modelo: se usa el ajuste_ind cargado, que estimacion.R ya ajustó sobre
# el conjunto con las 477 observaciones recuperadas.

cat("── 2. PREPROCESAMIENTO CON HORASH IMPUTADA ───────────────\n")

# Imputar HORASH con idéntico criterio que estimacion.R (función compartida)
# y predecir con el modelo YA AJUSTADO (sin reajustar).
data_para_mimic <- imputar_horash(data_original)

datos_mimic    <- preprocesar_datos(data_para_mimic, contratos_parciales)
data_modelo_ss <- datos_mimic$ss %>%
  mutate(
    Pred_Bartlett = predecir_bartlett(ajuste_ind, .),
    ln_Pred       = log(Pred_Bartlett)
  )

cat("N con predicción MIMIC:", nrow(data_modelo_ss),
    "(incluye las observaciones recuperadas por HORASH)\n\n")


# * 3. REGRESIÓN ROBUSTA (MINCER) + ESTRATOS DE VARIANZA ----
# Se ajusta UNA sola regresión robusta de Huber de la predicción (ln_Pred)
# sobre el perfil sociodemográfico —SIN indicadores ni variables de
# cotización, que las observaciones sin indicador apenas tienen (GC: 5.3%,
# Contrato: 40.7%; perfil demográfico: 100%)—. Al depender solo del perfil,
# este modelo puede predecir para CUALQUIER observación, incluidas las que
# no tienen indicador.
#
# El pool de residuos para el ruido se toma SOLO de casos limpios (ni
# atípicos ni topados), estratificados por Ocupación × Jornada_Completa ×
# cuartil de predicción. Así el ruido añadido es coherente con el perfil
# y el nivel salarial del caso a imputar, y no está contaminado por
# anomalías ni por el efecto de los topes de cotización.

cat("── 3. REGRESIÓN ROBUSTA (MINCER) + ESTRATOS ──────────────\n")

# Marcar topados y atípicos sobre data_modelo_ss (vía idx). Se traen las
# tres banderas de depuracion.R: la final, la de la robusta y la del IQR,
# para poder incluirlas en el dataset de salida.
data_modelo_ss <- data_modelo_ss %>%
  left_join(dataset_atipicos %>%
              select(idx, Flag_atipico_final, Flag_Robusta_bin, Flag_IQR_bin),
            by = "idx") %>%
  mutate(
    Flag_atipico_final = coalesce(Flag_atipico_final, 0L),
    Flag_Robusta_bin   = coalesce(Flag_Robusta_bin, 0L),
    Flag_IQR_bin       = coalesce(Flag_IQR_bin, 0L),
    Topado = Tope_ANO %in% c("1. BMAX SALARIO_SS_ANO", "2. BMIN SALARIO_SS_ANO") |
             Tope_MDR %in% c("1. BMAX SALARIO_SS_MDR", "2. BMIN SALARIO_SS_MDR"),
    Limpio = (Flag_atipico_final == 0) & !Topado
  )

# Variables de perfil para la Mincer (solo las que las observaciones sin
# indicador tienen al 100%: perfil demográfico puro, sin GC ni contrato)
causas_mincer <- c("Edad_Z", "Edad_Z2", "Exp_Z", "Exp_Z2",
                   "Jornada_Completa", "Mujer_Dummy", "Extranjero",
                   "Educ_Secundaria", "Educ_FP_Bach", "Educ_Superior",
                   "Ocup_Directivos", "Ocup_TecSup", "Ocup_TecAp",
                   "Ocup_Admin", "Ocup_Servicios", "Ocup_Artesanos",
                   "Ocup_Operadores", "Sector_Agri", "Sector_Ind", "Sector_Const")

formula_mincer <- as.formula(paste("ln_Pred ~", paste(causas_mincer, collapse = " + ")))

datos_limpios <- data_modelo_ss %>%
  filter(Limpio, complete.cases(select(., ln_Pred, all_of(causas_mincer))))

cat("N casos limpios para ajustar la Mincer:", nrow(datos_limpios), "\n")

rr <- rlm(formula_mincer, data = datos_limpios, psi = psi.huber,
          k = huber_k, maxit = 100)

# Predicción robusta para TODA la muestra con predicción MIMIC disponible
data_modelo_ss <- data_modelo_ss %>%
  mutate(pred_robusta = as.numeric(predict(rr, newdata = .)))

# Estratos de varianza: Ocupación × Jornada_Completa × tramo de predicción
# (n_tramos_estrato cuantiles, definido en config.R; 4 = cuartiles)
cortes_pred <- quantile(data_modelo_ss$pred_robusta,
                        probs = 0:n_tramos_estrato / n_tramos_estrato,
                        na.rm = TRUE)
data_modelo_ss <- data_modelo_ss %>%
  mutate(
    tramo_pred = cut(pred_robusta, breaks = cortes_pred,
                     include.lowest = TRUE, labels = FALSE),
    estrato    = interaction(Ocupacion, Jornada_Completa, tramo_pred, drop = TRUE)
  )

# Pool de residuos por estrato, SOLO de casos limpios
residuos_limpios <- data_modelo_ss %>%
  filter(Limpio) %>%
  mutate(residuo = ln_Pred - pred_robusta)

res_pool <- split(residuos_limpios$residuo, residuos_limpios$estrato)

cat("N estratos con pool de residuos:", length(res_pool), "\n")
cat("Tamaño mediano de pool por estrato:",
    round(median(lengths(res_pool)), 0), "residuos\n\n")


# * 4. IMPUTACIÓN POR RUIDO NO PARAMÉTRICO ----
# Para cada observación a imputar se suma a su predicción robusta un
# residuo tomado al azar del pool de su estrato (donante real). Se aplica
# muestreo por rechazo (hasta max_intentos_rechazo intentos): el valor
# imputado debe caer entre el suelo admisible (SMI/base mínima de
# cotización prorrateados por jornada) y la base máxima de cotización. Si
# tras agotar los intentos ninguno cae en rango, se usa como respaldo la
# predicción robusta acotada al rango admisible (fallback determinista).
#
# Casos a imputar:
#   (A) Atípicos con predicción MIMIC (Flag_atipico_final = 1)
#   (B) Observaciones sin indicador (nunca entraron al MIMIC): misma Mincer
#       + ruido usando solo su perfil.

cat("── 4. IMPUTACIÓN POR RUIDO NO PARAMÉTRICO ────────────────\n")

# Suelo admisible en escala log (SMI prorrateado por jornada). No se aplica
# techo: max_admisible devuelve Inf (ver config.R), por lo que la condición
# superior siempre se cumple y solo se rechaza por debajo del suelo.
imputar_una <- function(fila) {
  pool <- res_pool[[ as.character(fila$estrato) ]]
  li   <- log(smi_prorrateado(fila))
  ls   <- log(max_admisible(fila))   # Inf → sin techo

  if (!is.null(pool) && length(pool) > 0) {
    for (intento in seq_len(max_intentos_rechazo)) {
      val <- fila$pred_robusta + sample(pool, 1)
      if (val >= li && val <= ls) return(val)
    }
  }
  # Fallback determinista: predicción robusta con el suelo como mínimo
  max(fila$pred_robusta, li)
}

# --- (A) Atípicos con predicción MIMIC ---
idx_atipicos <- which(data_modelo_ss$Flag_atipico_final == 1)
sal_imp_atip <- vapply(idx_atipicos, function(i) {
  imputar_una(data_modelo_ss[i, ])
}, numeric(1))

data_modelo_ss$ln_sal_imputado <- NA_real_
data_modelo_ss$ln_sal_imputado[idx_atipicos] <- sal_imp_atip

cat("Atípicos imputados (con predicción MIMIC):", length(idx_atipicos), "\n")

# --- (B) Observaciones sin indicador ---
# Se reconstruye su perfil con preprocesar_datos_causas() (que NO filtra por
# indicador) y se quedan las que NO están en data_modelo_ss (no obtuvieron
# predicción MIMIC) y tienen el perfil de la Mincer completo.
datos_todos   <- preprocesar_datos_causas(data_para_mimic, contratos_parciales)
idx_con_mimic <- data_modelo_ss$idx
data_sin_ind  <- datos_todos %>%
  filter(!idx %in% idx_con_mimic,
         complete.cases(select(., all_of(causas_mincer))))

cat("N sin indicador a imputar (perfil completo):", nrow(data_sin_ind), "\n")

data_sin_ind <- data_sin_ind %>%
  mutate(
    pred_robusta = as.numeric(predict(rr, newdata = .)),
    tramo_pred   = cut(pred_robusta, breaks = cortes_pred,
                       include.lowest = TRUE, labels = FALSE),
    # Los que caen fuera del rango de cortes (pred extrema) se asignan al
    # tramo extremo más cercano para no quedar sin estrato
    tramo_pred   = case_when(
      is.na(tramo_pred) & pred_robusta <= cortes_pred[1] ~ 1L,
      is.na(tramo_pred) & pred_robusta >= cortes_pred[length(cortes_pred)] ~ as.integer(n_tramos_estrato),
      TRUE ~ as.integer(tramo_pred)
    ),
    estrato = interaction(Ocupacion, Jornada_Completa, tramo_pred, drop = TRUE)
  )

data_sin_ind$ln_sal_imputado <- vapply(seq_len(nrow(data_sin_ind)), function(i) {
  imputar_una(data_sin_ind[i, ])
}, numeric(1))

cat("Sin indicador imputados:", nrow(data_sin_ind), "\n\n")


# * 5. ENSAMBLADO DEL DATASET FINAL (27.363, CON FLAGS) ----
# Se reconstruye el dataset completo partiendo de data_original (intacto).
# Cada observación recibe:
#   Salario_final       — la predicción MIMIC si es válida y no atípica,
#                          o el valor imputado en caso contrario
#   Imputado            — 1 si el salario procede de imputación, 0 si es
#                          predicción MIMIC directa
#   Metodo_imputacion   — "MIMIC" / "Robusta+ruido (atipico)" /
#                          "Robusta+ruido (sin indicador)" / "No estimable"

cat("── 5. ENSAMBLADO DEL DATASET FINAL ───────────────────────\n")

tabla_mimic <- data_modelo_ss %>%
  transmute(
    idx,
    Salario_mimic      = Pred_Bartlett,
    ln_sal_imp_atip    = ln_sal_imputado,
    Flag_atipico_final,
    Flag_Robusta_bin,
    Flag_IQR_bin
  )

tabla_sinind <- data_sin_ind %>%
  transmute(idx, ln_sal_imp_sinind = ln_sal_imputado)

dataset_final <- data_original %>%
  select(idx, any_of(c("IDENCPERHOG", "NVIVI", "NPERS"))) %>%
  left_join(tabla_mimic,  by = "idx") %>%
  left_join(tabla_sinind, by = "idx") %>%
  mutate(
    # Salario_mimic_previo: predicción MIMIC ANTES de imputar (incluye el
    # valor original de los atípicos, para poder comparar antes/después).
    # Es NA para las observaciones sin indicador (no tenían predicción).
    Salario_mimic_previo = Salario_mimic,
    Salario_final = case_when(
      !is.na(ln_sal_imp_atip) & coalesce(Flag_atipico_final, 0L) == 1 ~ exp(ln_sal_imp_atip),
      !is.na(ln_sal_imp_sinind)                                        ~ exp(ln_sal_imp_sinind),
      !is.na(Salario_mimic)                                            ~ Salario_mimic,
      TRUE ~ NA_real_
    ),
    Imputado = case_when(
      !is.na(ln_sal_imp_atip) & coalesce(Flag_atipico_final, 0L) == 1 ~ 1L,
      !is.na(ln_sal_imp_sinind)                                        ~ 1L,
      !is.na(Salario_mimic)                                            ~ 0L,
      TRUE ~ NA_integer_
    ),
    Metodo_imputacion = case_when(
      !is.na(ln_sal_imp_atip) & coalesce(Flag_atipico_final, 0L) == 1 ~ "Robusta+ruido (atipico)",
      !is.na(ln_sal_imp_sinind)                                        ~ "Robusta+ruido (sin indicador)",
      !is.na(Salario_mimic)                                            ~ "MIMIC",
      TRUE ~ "No estimable"
    ),
    # Banderas de atípico (0/1), NA si la observación no pasó por detección
    Flag_atipico_final = coalesce(Flag_atipico_final, NA_integer_),
    Flag_atipico_robusta = coalesce(Flag_Robusta_bin, NA_integer_),
    Flag_atipico_iqr     = coalesce(Flag_IQR_bin, NA_integer_)
  )

cat("Dataset final:", nrow(dataset_final), "observaciones\n\n")
cat("Reparto por método:\n")
dataset_final %>% count(Metodo_imputacion) %>%
  mutate(Pct = round(n / sum(n) * 100, 1)) %>% print()
cat("\n")


# * 6. CHECKS DE VALIDEZ ----
# Verificaciones antes de dar la imputación por buena:
#   - Se conserva el total de observaciones (27.363)
#   - Ningún salario imputado negativo o nulo
#   - Ningún imputado por debajo del suelo (SMI prorrateado por jornada)
#   - Cuántas observaciones quedan como "No estimable" (perfil incompleto)

cat("── 6. CHECKS DE VALIDEZ ──────────────────────────────────\n")

n_total_ok  <- nrow(dataset_final) == nrow(data_original)
n_negativos <- sum(dataset_final$Salario_final <= 0, na.rm = TRUE)
n_no_estim  <- sum(dataset_final$Metodo_imputacion == "No estimable")

cat("  Total observaciones conservado:", n_total_ok,
    "(", nrow(dataset_final), "de", nrow(data_original), ")\n")
cat("  Salarios <= 0:", n_negativos,
    ifelse(n_negativos == 0, "(OK)", "(*** REVISAR ***)"), "\n")
cat("  Observaciones 'No estimable' (perfil incompleto):", n_no_estim, "\n")

# Check del suelo, solo sobre imputados: ningún valor imputado debe caer
# por debajo del SMI prorrateado según su jornada. No se comprueba techo
# (no existe tope de salario real).
chequeo_suelo <- dataset_final %>%
  filter(Imputado == 1) %>%
  left_join(data_original %>%
              mutate(Jornada_Completa = ifelse(as.numeric(PARCO1) == 1, 1, 0)) %>%
              select(idx, Jornada_Completa), by = "idx") %>%
  rowwise() %>%
  mutate(
    suelo          = smi_prorrateado(cur_data()),
    por_debajo_smi = Salario_final < suelo
  ) %>%
  ungroup()

cat("  Imputados por debajo del SMI prorrateado:",
    sum(chequeo_suelo$por_debajo_smi, na.rm = TRUE), "\n\n")


# * 7. SUMMARY Y GRÁFICOS ----
# Resumen numérico del alcance de la imputación y gráficos de diagnóstico.
# Los gráficos "antes vs después" usan como "antes" la predicción MIMIC
# (Salario_mimic, que para los atípicos es su predicción original antes de
# reimputar) y como "después" el Salario_final. Permiten comprobar
# visualmente que la imputación no distorsiona la distribución global.

cat("── 7. SUMMARY Y GRÁFICOS ─────────────────────────────────\n")

# --- Summary numérico del alcance de la imputación ---
n_imputados <- sum(dataset_final$Imputado, na.rm = TRUE)
n_atipicos  <- sum(dataset_final$Metodo_imputacion == "Robusta+ruido (atipico)")
n_sinind    <- sum(dataset_final$Metodo_imputacion == "Robusta+ruido (sin indicador)")
n_mimic     <- sum(dataset_final$Metodo_imputacion == "MIMIC")

cat("\nResumen de la imputación:\n")
cat("  N total dataset original:      ", nrow(dataset_final), "\n")
cat("  N predicción MIMIC directa:    ", n_mimic,
    sprintf(" (%.1f%%)\n", n_mimic / nrow(dataset_final) * 100))
cat("  N imputados (total):           ", n_imputados,
    sprintf(" (%.1f%%)\n", n_imputados / nrow(dataset_final) * 100))
cat("    - Atípicos reimputados:      ", n_atipicos,
    sprintf(" (%.1f%%)\n", n_atipicos / nrow(dataset_final) * 100))
cat("    - Sin indicador (perfil):    ", n_sinind,
    sprintf(" (%.1f%%)\n", n_sinind / nrow(dataset_final) * 100))

cat("\nEstadísticos del salario final por método:\n")
dataset_final %>%
  filter(!is.na(Salario_final)) %>%
  group_by(Metodo_imputacion) %>%
  summarise(
    N       = n(),
    Media   = round(mean(Salario_final), 0),
    Mediana = round(median(Salario_final), 0),
    SD      = round(sd(Salario_final), 0),
    Min     = round(min(Salario_final), 0),
    Max     = round(max(Salario_final), 0),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("\n")

# Percentil 99 conjunto para recortar colas en los gráficos (visualización)
p99_sal <- quantile(c(dataset_final$Salario_mimic, dataset_final$Salario_final),
                    0.99, na.rm = TRUE)

# --- IMPUTA_G1: distribución del salario antes vs después (TODOS los datos) ---
# "Antes" = predicción MIMIC (solo disponible para quien tiene predicción);
# "Después" = Salario_final (todas las observaciones, imputadas incluidas).
datos_g1 <- bind_rows(
  dataset_final %>%
    filter(!is.na(Salario_mimic)) %>%
    transmute(Salario = Salario_mimic, Momento = "Antes (MIMIC)"),
  dataset_final %>%
    filter(!is.na(Salario_final)) %>%
    transmute(Salario = Salario_final, Momento = "Después (imputado)")
) %>%
  filter(Salario < p99_sal)

g1_imp <- datos_g1 %>%
  ggplot(aes(x = Salario, color = Momento, fill = Momento)) +
  geom_density(alpha = 0.2, linewidth = 0.9) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = c("Antes (MIMIC)" = "#534AB7",
                                "Después (imputado)" = "#D85A30")) +
  scale_fill_manual(values  = c("Antes (MIMIC)" = "#534AB7",
                                "Después (imputado)" = "#D85A30")) +
  labs(x = "Salario mensual (euros/mes)", y = "Densidad", color = NULL, fill = NULL,
       title = "Distribución del salario antes vs después de la imputación",
       subtitle = "Todos los datos") +
  theme_minimal() +
  theme(legend.position = "bottom")
guardar_grafico(g1_imp, ruta_out_imputa, "IMPUTA_G1_antes_despues_todos")

# --- IMPUTA_G2: distribución antes vs después SOLO para los atípicos ---
# "Antes" = predicción MIMIC del atípico (el valor que lo hizo atípico);
# "Después" = su valor reimputado.
idx_atip_final <- dataset_final %>%
  filter(Metodo_imputacion == "Robusta+ruido (atipico)") %>% pull(idx)

datos_g2 <- bind_rows(
  dataset_final %>%
    filter(idx %in% idx_atip_final, !is.na(Salario_mimic)) %>%
    transmute(Salario = Salario_mimic, Momento = "Antes (MIMIC atípico)"),
  dataset_final %>%
    filter(idx %in% idx_atip_final, !is.na(Salario_final)) %>%
    transmute(Salario = Salario_final, Momento = "Después (imputado)")
) %>%
  filter(Salario < quantile(Salario, 0.99, na.rm = TRUE))

g2_imp <- datos_g2 %>%
  ggplot(aes(x = Salario, color = Momento, fill = Momento)) +
  geom_density(alpha = 0.2, linewidth = 0.9) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = c("Antes (MIMIC atípico)" = "#534AB7",
                                "Después (imputado)" = "#D85A30")) +
  scale_fill_manual(values  = c("Antes (MIMIC atípico)" = "#534AB7",
                                "Después (imputado)" = "#D85A30")) +
  labs(x = "Salario mensual (euros/mes)", y = "Densidad", color = NULL, fill = NULL,
       title = "Distribución del salario antes vs después de la imputación",
       subtitle = "Solo los atípicos detectados") +
  theme_minimal() +
  theme(legend.position = "bottom")
guardar_grafico(g2_imp, ruta_out_imputa, "IMPUTA_G2_antes_despues_atipicos")

# --- IMPUTA_G3: distribución del salario final por método de obtención ---
# Muestra si los valores imputados (atípicos y sin indicador) son coherentes
# con los estimados directamente por el MIMIC o están desplazados.
g3_imp <- dataset_final %>%
  filter(!is.na(Salario_final), Metodo_imputacion != "No estimable",
         Salario_final < p99_sal) %>%
  ggplot(aes(x = Salario_final, color = Metodo_imputacion, fill = Metodo_imputacion)) +
  geom_density(alpha = 0.15, linewidth = 0.9) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_viridis_d(end = 0.85) +
  scale_fill_viridis_d(end = 0.85) +
  labs(x = "Salario final (euros/mes)", y = "Densidad", color = NULL, fill = NULL,
       title = "Distribución del salario final por método de obtención",
       subtitle = "¿Son los valores imputados coherentes con los del MIMIC?") +
  theme_minimal() +
  theme(legend.position = "bottom")
guardar_grafico(g3_imp, ruta_out_imputa, "IMPUTA_G3_distribucion_por_metodo")

# --- IMPUTA_G4: desplazamiento de los atípicos (antes → después) ---
# Diferencia entre el valor reimputado y la predicción MIMIC original de
# cada atípico. Muestra en qué dirección y magnitud corrige la imputación.
g4_imp <- dataset_final %>%
  filter(Metodo_imputacion == "Robusta+ruido (atipico)",
         !is.na(Salario_mimic), !is.na(Salario_final)) %>%
  mutate(Desplazamiento = Salario_final - Salario_mimic) %>%
  filter(abs(Desplazamiento) < quantile(abs(Desplazamiento), 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Desplazamiento)) +
  geom_histogram(bins = 60, fill = "#1D9E75", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#D85A30", linewidth = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Salario imputado − MIMIC original (euros/mes)", y = "Frecuencia",
       title = "Desplazamiento de los atípicos tras la imputación",
       subtitle = "Negativo: la imputación corrige a la baja; positivo: al alza") +
  theme_minimal()
guardar_grafico(g4_imp, ruta_out_imputa, "IMPUTA_G4_desplazamiento_atipicos")

cat("\n")


# * 8. GUARDADO ----

cat("── 8. GUARDADO ───────────────────────────────────────────\n")

dataset_export <- dataset_final %>%
  select(idx, any_of(c("IDENCPERHOG", "NVIVI", "NPERS")),
         Salario_mimic_previo, Salario_final, Imputado, Metodo_imputacion,
         Flag_atipico_final, Flag_atipico_robusta, Flag_atipico_iqr)

saveRDS(dataset_export,
        file.path(ruta_out_imputa, paste0("dataset_imputado_", tag_fecha, ".rds")))
# write_csv2: separador ';' y coma decimal, para que Excel con configuración
# regional española abra el fichero correctamente (write_csv usa punto
# decimal, que Excel-ES malinterpreta como separador de miles).
write_csv2(dataset_export,
           file.path(ruta_out_imputa, paste0("dataset_imputado_", tag_fecha, ".csv")))

cat("Guardado:", nrow(dataset_export), "filas en", ruta_out_imputa, "\n")
cat("Columnas: idx, [ident.], Salario_mimic_previo (MIMIC antes de imputar),\n")
cat("          Salario_final, Imputado, Metodo_imputacion,\n")
cat("          Flag_atipico_final, Flag_atipico_robusta, Flag_atipico_iqr\n\n")

cat("=============================================================\n")
cat(" IMPUTACIÓN COMPLETADA\n")
cat(" N total:", nrow(dataset_final),
    "| Imputados:", sum(dataset_final$Imputado, na.rm = TRUE),
    "| No estimables:", n_no_estim, "\n")
cat("=============================================================\n")

sink()
