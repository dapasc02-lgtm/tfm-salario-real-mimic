# =============================================================================
# DEPURACION.R — DETECCIÓN DE ATÍPICOS SOBRE LA PREDICCIÓN DEL MIMIC
# =============================================================================
#
# Se ejecuta DESPUÉS de estimacion.R (necesita outputs/modelos_ajustados.rds).
#
# MÉTODOS:
#   1. Regresión robusta (residuo estudentizado heterocedástico) — PRINCIPAL.
#      Cada observación se compara contra su propia varianza esperada,
#      según cuántos indicadores (AT/SS/MDR) tiene disponibles. Detalle
#      matemático completo en funciones_depuracion.R.
#   2. IQR / Eurostat (Tukey) — REFUERZO. No decide por sí solo bajo
#      criterio_atipico_final = "solo_robusta".
#
# DESCARTADOS del criterio final (por bajo rendimiento bajo FDR riguroso):
#   - Isolation Forest (IF-1/IF-2/IF-3)
#   - Detector de discrepancias múltiples (error_AT/SS/MDR)
#
# CRITERIO FINAL: Flag_atipico_final, configurable vía
# criterio_atipico_final en config.R ("solo_robusta" / "robusta_o_iqr" /
# "robusta_y_iqr").
#
# PENDIENTE: imputación no implementada — este script es solo DETECCIÓN.
#
# Estructura: 1.Carga+preprocesamiento 
#             2.Robusta 
#             3.IQR 
#             4.Flag_atipico_final
#             5.Contingencia 
#             6.Guardado 
#             7.Gráficos
# =============================================================================

library(tidyverse)
library(readxl)
library(lavaan)
library(MASS)
library(moments)
library(here)

select <- dplyr::select

source(here("config.R"))
source(here("funciones_modelo.R"))
source(here("funciones_depuracion.R"))

crear_carpetas_salida()
crear_carpetas_depuracion()

fecha_ejecucion <- Sys.time()
tag_fecha       <- format(fecha_ejecucion, "%Y%m%d_%H%M")
sink(file.path(ruta_log, paste0("log_depuracion_", tag_fecha, ".txt")),
     split = TRUE, append = FALSE)

cat("=============================================================\n")
cat(" DEPURACIÓN — DETECCIÓN DE ATÍPICOS\n")
cat(" Versión:", version_pipeline, "| Fecha:", format(fecha_ejecucion),
    "| Criterio:", criterio_atipico_final, "\n")
cat("=============================================================\n\n")


# * 1. CARGA Y PREPROCESAMIENTO ----
# Se carga el modelo MIMIC ya ajustado por estimacion.R (no se reajusta
# aquí) y el Excel original. El preprocesamiento es idéntico al de
# estimacion.R (misma función preprocesar_datos), garantizando que
# trabajamos sobre exactamente el mismo conjunto de observaciones y
# variables. Se añaden las discrepancias entre fuentes (Dif_AT_SS, etc.)
# y se calcula la predicción Bartlett que servirá de base a los dos
# métodos de detección.

cat("── 1. CARGA DEL MODELO Y PREPROCESAMIENTO ─────────────────\n")

if (!file.exists(ruta_modelos_ajust)) {
  stop("No se encuentra el modelo MIMIC ajustado en:\n  ", ruta_modelos_ajust,
       "\nEjecutar primero estimacion.R.")
}

modelos_guardados <- readRDS(ruta_modelos_ajust)
ajuste_ind        <- modelos_guardados$ajuste_ind

cat("Modelo MIMIC ajustado el:", format(modelos_guardados$fecha_ajuste), "\n")

if (!file.exists(ruta_entrada)) stop("No se encuentra el fichero de entrada en: ", ruta_entrada)

data_salarios <- read_excel(ruta_entrada)
validar_columnas(data_salarios, columnas_esperadas)

datos          <- preprocesar_datos(data_salarios, contratos_parciales)
data_modelo_ss <- datos$ss %>% calcular_discrepancias()

cat("N Excel original:", nrow(data_salarios),
    "| N tras preprocesamiento:", nrow(data_modelo_ss),
    sprintf("| N descartado: %d (%.1f%%)\n",
            nrow(data_salarios) - nrow(data_modelo_ss),
            (nrow(data_salarios) - nrow(data_modelo_ss)) / nrow(data_salarios) * 100))

data_modelo_ss <- data_modelo_ss %>%
  mutate(
    Pred_Bartlett = predecir_bartlett(ajuste_ind, data_modelo_ss),
    ln_Pred       = log(Pred_Bartlett)
  )

cat("Pred_Bartlett — Media:", round(mean(data_modelo_ss$Pred_Bartlett, na.rm=TRUE), 0),
    "| Mediana:", round(median(data_modelo_ss$Pred_Bartlett, na.rm=TRUE), 0),
    "| SD:", round(sd(data_modelo_ss$Pred_Bartlett, na.rm=TRUE), 0), "euros/mes\n\n")


# * 2. REGRESIÓN ROBUSTA (PRINCIPAL) ----
# Método principal de detección. Se ajusta una regresión robusta de Huber
# (rlm, M-estimación con k=1.345) de ln_Pred sobre las causas estructurales
# del MIMIC. El residuo de cada observación se estudentiza de forma
# heterocedástica: la varianza esperada depende del número de indicadores
# disponibles (AT, SS, MDR), porque cuantos más indicadores tiene una
# observación, más precisa es su predicción Bartlett y menor es su
# varianza de medición. El p-valor empírico de cada residuo se obtiene
# via ECDF (sin asumir normalidad) y se corrige por FDR
# (Benjamini-Hochberg) para controlar la tasa de falsos positivos al
# evaluar simultáneamente ~25.000 contrastes. El criterio final es
# Flag_Robusta_bin = 1 si p_valor_bh < alpha (configurable en config.R).

cat("── 2. REGRESIÓN ROBUSTA (MÉTODO PRINCIPAL) ────────────────\n")

df_robusta <- data_modelo_ss %>%
  select(idx, ln_AT_mens, ln_SS_mens, ln_MDR, all_of(causas_base),
         Pred_Bartlett, Edad, Sexo_Orig, Meses_Emp, Ocupacion, Sector,
         Jornada_Completa, Educ_Superior, Contrato_Parc) %>%
  filter(complete.cases(select(., all_of(causas_base))))

cat("N disponible:", nrow(df_robusta), "\n")

resultado_het <- calcular_residuo_heterocedastico(
  ajuste = ajuste_ind, datos = df_robusta, causas = causas_base,
  indicadores = c("ln_AT_mens", "ln_SS_mens", "ln_MDR"),
  factor_latente = "Salario_Real", metodo_rlm = "M",
  k_huber = huber_k, alpha = huber_alpha_bonferroni
)

df_robusta <- bind_cols(df_robusta, resultado_het) %>%
  mutate(
    Flag_Robusta     = ifelse(flag_bin == 1, "Atípico", "Normal"),
    Flag_Robusta_bin = flag_bin
  )

cat("psi_hat:", round(resultado_het$psi_hat[1], 5),
    "| Media 1/I_i:", round(mean(df_robusta$var_medicion, na.rm=TRUE), 5), "\n")

cat("\nVarianza de medición por n.º de indicadores disponibles:\n")
df_robusta %>%
  group_by(n_indicadores_validos) %>%
  summarise(N = n(), Media_var_medicion = round(mean(var_medicion, na.rm=TRUE), 5),
            .groups = "drop") %>%
  print()

cat("\nEstadísticos de z (residuo estudentizado):\n")
df_robusta %>%
  summarise(Min = round(min(z, na.rm=TRUE), 2), P5 = round(quantile(z, .05, na.rm=TRUE), 2),
            P50 = round(median(z, na.rm=TRUE), 2), P95 = round(quantile(z, .95, na.rm=TRUE), 2),
            Max = round(max(z, na.rm=TRUE), 2)) %>%
  print()

cat("\nAtípicos (FDR, alpha=", huber_alpha_bonferroni, "):",
    sum(df_robusta$Flag_Robusta_bin),
    sprintf(" (%.1f%%)\n", mean(df_robusta$Flag_Robusta_bin) * 100))

cat("\nReparto por n.º de indicadores entre los atípicos detectados:\n")
df_robusta %>%
  filter(Flag_Robusta_bin == 1) %>%
  count(n_indicadores_validos) %>%
  mutate(Pct = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\nDe los atípicos con 1 solo indicador, cuál es ese indicador:\n")
df_robusta %>%
  filter(Flag_Robusta_bin == 1, n_indicadores_validos == 1) %>%
  mutate(Fuente_unica = case_when(
    !is.na(ln_AT_mens) ~ "Solo AT", !is.na(ln_SS_mens) ~ "Solo SS",
    !is.na(ln_MDR) ~ "Solo MDR", TRUE ~ "Error")) %>%
  count(Fuente_unica) %>%
  mutate(Pct = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\nTop 15 casos con mayor |z| (Pred en miles de euros/mes):\n")
df_robusta %>%
  arrange(desc(abs(z))) %>%
  head(15) %>%
  mutate(Pred_miles = round(Pred_Bartlett / 1000, 2)) %>%
  select(idx, z, n_indicadores_validos, Pred_miles, Edad, Meses_Emp,
         Ocupacion, Sector, Jornada_Completa, Contrato_Parc) %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("\n")


# * 3. IQR / EUROSTAT (REFUERZO) ----
# Método de refuerzo basado en el criterio de Tukey (IQR) aplicado sobre
# el residuo de una regresión robusta auxiliar de ln_Pred sobre las causas.
# Se calculan dos umbrales: normal (Q1 - 1.5·IQR, Q3 + 1.5·IQR) y extremo
# (Q1 - 3·IQR, Q3 + 3·IQR), factores configurables en config.R.
# Este método no decide de forma autónoma bajo el criterio por defecto
# ("solo_robusta"): su función es reforzar la confianza cuando coincide
# con la robusta (nivel de severidad "Robusta + IQR") o señalar casos
# que la robusta no alcanza cuando se usa "robusta_o_iqr".

cat("── 3. IQR / EUROSTAT (MÉTODO DE REFUERZO) ─────────────────\n")

df_iqr <- data_modelo_ss %>%
  select(idx, ln_Pred, all_of(causas_base), Pred_Bartlett) %>%
  filter(complete.cases(select(., ln_Pred, all_of(causas_base))))

formula_aux <- as.formula(paste("ln_Pred ~", paste(causas_base, collapse = " + ")))
ajuste_aux  <- rlm(formula_aux, data = df_iqr, psi = psi.huber, k = huber_k, maxit = 100)
df_iqr <- df_iqr %>% mutate(Residuo = residuals(ajuste_aux))

calcular_limites_iqr <- function(x, factor) {
  q1 <- quantile(x, 0.25, na.rm = TRUE); q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  c(q1 - factor * iqr, q3 + factor * iqr)
}
lim_normal  <- calcular_limites_iqr(df_iqr$Residuo, iqr_factor_normal)
lim_extremo <- calcular_limites_iqr(df_iqr$Residuo, iqr_factor_extremo)

cat("N disponible:", nrow(df_iqr), "\n")
cat("Límites Tukey — Normal: [", round(lim_normal[1],3), ",", round(lim_normal[2],3),
    "] | Extremo: [", round(lim_extremo[1],3), ",", round(lim_extremo[2],3), "]\n")

df_iqr <- df_iqr %>%
  mutate(
    Flag_IQR = case_when(
      Residuo < lim_extremo[1] | Residuo > lim_extremo[2] ~ "Atípico extremo",
      Residuo < lim_normal[1]  | Residuo > lim_normal[2]  ~ "Atípico",
      TRUE ~ "Normal"
    ),
    Flag_IQR_bin = flag_a_binario(Flag_IQR)
  )

cat("\nDistribución:\n")
df_iqr %>% count(Flag_IQR) %>% mutate(Pct = round(n/sum(n)*100, 2)) %>% print()
cat("\nTotal atípicos (Tukey):", sum(df_iqr$Flag_IQR_bin),
    sprintf(" (%.1f%%)\n\n", mean(df_iqr$Flag_IQR_bin) * 100))


# * 4. FLAG_ATIPICO_FINAL ----
# Se combinan los resultados de los dos métodos anteriores en una única
# bandera de atípico final, siguiendo el criterio configurado en config.R.
# Además de la bandera binaria, se calcula Nivel_severidad para distinguir
# cuándo ambos métodos coinciden (mayor confianza) de cuándo solo uno de
# ellos detecta la anomalía. Se imprime también un perfil comparado entre
# normales y atípicos (predicción media, jornada, educación, ocupación)
# para facilitar la interpretación cualitativa de los casos detectados.

cat("── 4. FLAG_ATIPICO_FINAL (criterio:", criterio_atipico_final, ") ──\n")
cat("Opciones: solo_robusta | robusta_o_iqr | robusta_y_iqr\n\n")

opciones_validas <- c("solo_robusta", "robusta_o_iqr", "robusta_y_iqr")
if (!criterio_atipico_final %in% opciones_validas) {
  stop("criterio_atipico_final no reconocido: '", criterio_atipico_final,
       "'. Debe ser uno de: ", paste(opciones_validas, collapse = ", "))
}

data_final <- data_modelo_ss %>%
  left_join(df_robusta %>% select(idx, z, p_valor_bh, n_indicadores_validos,
                                  Flag_Robusta, Flag_Robusta_bin), by = "idx") %>%
  left_join(df_iqr %>% select(idx, Residuo, Flag_IQR, Flag_IQR_bin), by = "idx")

data_final <- data_final %>%
  mutate(
    Flag_atipico_final = case_when(
      criterio_atipico_final == "solo_robusta" ~ Flag_Robusta_bin,
      criterio_atipico_final == "robusta_o_iqr" ~
        as.integer(coalesce(Flag_Robusta_bin, 0L) == 1 | coalesce(Flag_IQR_bin, 0L) == 1),
      criterio_atipico_final == "robusta_y_iqr" ~
        as.integer(coalesce(Flag_Robusta_bin, 0L) == 1 & coalesce(Flag_IQR_bin, 0L) == 1)
    ),
    Nivel_severidad = case_when(
      is.na(Flag_atipico_final) | Flag_atipico_final == 0 ~ "Normal",
      coalesce(Flag_Robusta_bin, 0L) == 1 & coalesce(Flag_IQR_bin, 0L) == 1 ~ "Robusta + IQR",
      coalesce(Flag_Robusta_bin, 0L) == 1 ~ "Solo Robusta",
      coalesce(Flag_IQR_bin, 0L) == 1 ~ "Solo IQR",
      TRUE ~ NA_character_
    )
  )

cat("Distribución de Flag_atipico_final:\n")
data_final %>% count(Flag_atipico_final) %>% mutate(Pct = round(n/sum(n)*100, 2)) %>% print()

cat("\nNiveles de severidad:\n")
data_final %>% filter(Flag_atipico_final == 1) %>%
  count(Nivel_severidad) %>% mutate(Pct = round(n/sum(n)*100, 1)) %>% print()

cat("\nPerfil comparado (Normal vs Atípico):\n")
data_final %>%
  group_by(Flag_atipico_final) %>%
  summarise(N = n(), Media_pred = round(mean(Pred_Bartlett, na.rm=TRUE), 0),
            SD_pred = round(sd(Pred_Bartlett, na.rm=TRUE), 0),
            Pct_jornada_completa = round(mean(Jornada_Completa, na.rm=TRUE)*100, 1),
            Pct_educ_superior = round(mean(Educ_Superior, na.rm=TRUE)*100, 1),
            Pct_directivos = round(mean(Ocup_Directivos, na.rm=TRUE)*100, 1),
            .groups = "drop") %>%
  print()
cat("\n")


# * 5. TABLA DE CONTINGENCIA (ROBUSTA x IQR) ----
# Cuantifica el solapamiento entre los dos métodos: cuántos casos detecta
# solo la robusta, solo el IQR, ambos a la vez, o ninguno. Permite
# evaluar en qué medida los métodos son complementarios o redundantes.
# La variable Combinacion se guarda en el dataset final para que sea
# posible filtrar subgrupos en análisis posteriores.

cat("── 5. TABLA DE CONTINGENCIA (ROBUSTA x IQR) ────────────────\n")

niveles_combinacion <- c("Ninguno", "Solo Robusta", "Solo IQR", "Robusta + IQR")

data_final <- data_final %>%
  mutate(
    Combinacion = case_when(
      coalesce(Flag_Robusta_bin, 0L) == 0 & coalesce(Flag_IQR_bin, 0L) == 0 ~ "Ninguno",
      coalesce(Flag_Robusta_bin, 0L) == 1 & coalesce(Flag_IQR_bin, 0L) == 0 ~ "Solo Robusta",
      coalesce(Flag_Robusta_bin, 0L) == 0 & coalesce(Flag_IQR_bin, 0L) == 1 ~ "Solo IQR",
      coalesce(Flag_Robusta_bin, 0L) == 1 & coalesce(Flag_IQR_bin, 0L) == 1 ~ "Robusta + IQR",
      TRUE ~ NA_character_
    ),
    Combinacion = factor(Combinacion, levels = niveles_combinacion)
  )

data_final %>%
  filter(!is.na(Combinacion)) %>%
  count(Combinacion, .drop = FALSE) %>%
  mutate(Pct = round(n / sum(n) * 100, 3)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("\n")


# * 6. GUARDAR RESULTADOS ----
# Se exporta un dataset reducido con las columnas de detección (flags,
# residuos, nivel de severidad y combinación de métodos) junto al idx
# y la predicción Bartlett, en formato .rds (para uso en R) y .csv
# (para consulta externa). El nombre del fichero incluye la fecha y hora
# de ejecución para mantener un histórico de cada pasada.

cat("── 6. GUARDAR RESULTADOS ──────────────────────────────────\n")

dataset_atipicos <- data_final %>%
  select(idx, Pred_Bartlett, z, p_valor_bh, n_indicadores_validos,
         Flag_Robusta, Flag_Robusta_bin, Residuo, Flag_IQR, Flag_IQR_bin,
         Flag_atipico_final, Nivel_severidad, Combinacion)

saveRDS(dataset_atipicos,
        file.path(ruta_out_depura, paste0("dataset_atipicos_", tag_fecha, ".rds")))
write_csv(dataset_atipicos,
          file.path(ruta_out_depura, paste0("dataset_atipicos_", tag_fecha, ".csv")))

cat("Guardado:", nrow(dataset_atipicos), "filas en", ruta_out_depura, "\n\n")


# * 7. GRÁFICOS ----
# Tres gráficos de diagnóstico:
#   G1: distribución del residuo estudentizado z — permite ver si la masa
#       central es simétrica y si hay colas pesadas que justifiquen la
#       detección de atípicos.
#   G2: composición de los atípicos detectados por combinación de métodos
#       (Solo Robusta / Solo IQR / Robusta+IQR) — visualiza el solapamiento
#       entre métodos de la sección 5.
#   G3: scatter de z vs predicción Bartlett, coloreado por Flag_atipico_final
#       — permite ver si los atípicos se concentran en algún tramo salarial
#       o si están distribuidos a lo largo de toda la escala.

cat("── 7. GRÁFICOS ────────────────────────────────────────────\n")

g1 <- data_final %>%
  filter(!is.na(z)) %>%
  ggplot(aes(x = z)) +
  geom_histogram(bins = 80, fill = "#534AB7", color = "white", alpha = 0.8) +
  labs(x = "z (residuo estudentizado heterocedástico)", y = "Frecuencia",
       title = "Distribución del residuo de la regresión robusta") +
  theme_minimal()
guardar_grafico_depuracion(g1, ruta_out_depura, "DEPURA_distribucion_z")

g2 <- data_final %>%
  filter(!is.na(Combinacion), Combinacion != "Ninguno") %>%
  count(Combinacion, .drop = FALSE) %>%
  ggplot(aes(x = Combinacion, y = n, fill = Combinacion)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  scale_fill_viridis_d(guide = "none") +
  labs(x = NULL, y = "N.º de trabajadores",
       title = "Composición de los atípicos detectados (Robusta vs IQR)") +
  theme_minimal()
guardar_grafico_depuracion(g2, ruta_out_depura, "DEPURA_composicion")

g3 <- data_final %>%
  filter(Pred_Bartlett < quantile(Pred_Bartlett, 0.99, na.rm = TRUE)) %>%
  mutate(Flag_lab = factor(ifelse(Flag_atipico_final == 1, "Atípico", "Normal"),
                           levels = c("Normal", "Atípico"))) %>%
  arrange(Flag_lab) %>%
  ggplot(aes(x = Pred_Bartlett, y = z, color = Flag_lab)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_manual(values = c("Normal" = "#2C7FB8", "Atípico" = "#D7191C")) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Predicción Bartlett (€/mes)", y = "z (residuo studentizado)",
       color = "Flag_atipico_final",
       title = "Residuo studentizado vs predicción") +
  theme_minimal() + theme(legend.position = "bottom")
guardar_grafico_depuracion(g3, ruta_out_depura, "DEPURA_z_vs_pred")

cat("\nGráficos guardados en", ruta_out_depura, "\n\n")

cat("=============================================================\n")
cat(" DEPURACIÓN COMPLETADA\n")
cat(" N evaluado:", nrow(data_final), "| N atípicos:",
    sum(dataset_atipicos$Flag_atipico_final, na.rm = TRUE),
    "| Criterio:", criterio_atipico_final, "\n")
cat(" PENDIENTE: imputación no implementada en este script.\n")
cat("=============================================================\n")

sink()
