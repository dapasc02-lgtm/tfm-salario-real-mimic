# =============================================================================
# FUNCIONES_DEPURACION.R — FUNCIONES DE DETECCIÓN DE ATÍPICOS
# =============================================================================
#
# Funciones puras de detección, usadas por depuracion.R. NO incluye
# funciones de preprocesamiento ni de ajuste del MIMIC — esas ya existen
# en funciones_modelo.R (preprocesar_datos, ajustar_mimic) y se
# reutilizan directamente, evitando duplicar lógica entre el pipeline
# de atípicos original y el de producción.
#
# Funciones incluidas:
#   calcular_discrepancias()          — diferencias entre fuentes (IQR)
#   flag_a_binario()                  — conversión de flag 3 niveles a 0/1
#   calcular_bonferroni_ecdf()        — corrección Bonferroni vía ECDF
#   calcular_fdr_ecdf()               — corrección FDR (Benjamini-Hochberg)
#   calcular_residuo_heterocedastico() — residuo studentizado de la robusta
#   guardar_grafico_depuracion()      — guardar gráfico ggplot a disco
# =============================================================================


# ── Discrepancias entre fuentes (usado por el IQR/Eurostat) ────────────────
calcular_discrepancias <- function(data_modelo_ss) {
  data_modelo_ss %>%
    mutate(
      Dif_AT_SS  = ln_AT_mens - ln_SS_mens,
      Dif_AT_MDR = ln_AT_mens - ln_MDR,
      Dif_SS_MDR = ln_SS_mens - ln_MDR
    )
}


# ── Conversión de flag de 3 niveles a bandera binaria ───────────────────────
# 1 si el caso es "Atípico" o "Atípico extremo", 0 si es "Normal".
flag_a_binario <- function(flag_3niveles) {
  as.integer(flag_3niveles %in% c("Atípico", "Atípico extremo"))
}


# ── Corrección de Bonferroni vía ECDF ───────────────────────────────────────
# Calcula un p-valor empírico a partir de la función de distribución
# empírica (ECDF) de la propia variable, en lugar de asumir una
# distribución teórica (normal) que puede no ser correcta.
#
# x_poblacion (opcional): si se proporciona, la ECDF se construye sobre
# esa población de referencia completa, y x solo determina cuántos
# tests se corrigen (n en el denominador alpha/n). Esto evita que el
# p-valor de cada caso quede artificialmente inflado al compararse solo
# contra un subconjunto ya preseleccionado, en lugar de contra la
# población completa.
calcular_bonferroni_ecdf <- function(x, alpha = 0.05, colas = "ambas",
                                     x_poblacion = NULL) {

  n <- length(x)
  alpha_corregido <- alpha / n

  if (is.null(x_poblacion)) x_poblacion <- x

  if (colas == "ambas") {
    x_ref    <- abs(x)
    x_ref_pob <- abs(x_poblacion)
  } else if (colas == "superior") {
    x_ref    <- x
    x_ref_pob <- x_poblacion
  } else {
    stop("colas debe ser 'ambas' o 'superior'")
  }

  ecdf_pob <- ecdf(x_ref_pob)
  p_valor  <- 1 - ecdf_pob(x_ref)

  umbral <- quantile(x_ref_pob, probs = 1 - alpha_corregido, na.rm = TRUE)

  list(
    p_valor         = p_valor,
    umbral          = umbral,
    alpha_corregido = alpha_corregido,
    flag_bin        = as.integer(p_valor < alpha_corregido)
  )
}


# ── Corrección de FDR (Benjamini-Hochberg) vía ECDF ─────────────────────────
# Alternativa a Bonferroni para el mismo problema (corrección por
# comparaciones múltiples), usando un umbral CRECIENTE con el rango del
# p-valor en lugar de un umbral fijo para todos:
#
#   Bonferroni: umbral fijo para todos        = alpha / n
#   FDR (BH):   umbral creciente con el rango = (i / n) * alpha
#
# Referencia: Benjamini, Y. y Hochberg, Y. (1995). Controlling the false
# discovery rate. Journal of the Royal Statistical Society B, 57(1).
#
# Este es el CRITERIO FINAL usado para clasificar atípicos en
# depuracion.R (la robusta lo aplica internamente sobre el p-valor de z;
# el IQR lo aplica sobre el residuo de la regresión auxiliar).
calcular_fdr_ecdf <- function(x, alpha = 0.05, colas = "ambas",
                              x_poblacion = NULL) {

  n <- length(x)

  if (is.null(x_poblacion)) x_poblacion <- x

  if (colas == "ambas") {
    x_ref     <- abs(x)
    x_ref_pob <- abs(x_poblacion)
  } else if (colas == "superior") {
    x_ref     <- x
    x_ref_pob <- x_poblacion
  } else {
    stop("colas debe ser 'ambas' o 'superior'")
  }

  ecdf_pob <- ecdf(x_ref_pob)
  p_valor  <- 1 - ecdf_pob(x_ref)

  p_ajustado <- p.adjust(p_valor, method = "BH")

  list(
    p_valor    = p_valor,
    p_ajustado = p_ajustado,
    flag_bin   = as.integer(p_ajustado < alpha)
  )
}


# ── Residuo estudentizado heterocedástico (regresión robusta) ──────────────
# Corrige el error conceptual de tratar ln_Pred (una ESTIMACIÓN de
# Bartlett) como si fuera el dato real al ajustar la regresión robusta.
#
# DESCOMPOSICIÓN: e_i = u_i + zeta_i
#   u_i    : error de medición de Bartlett, Var(u_i) = 1/I_i (distinta
#            por observación, según indicadores AT/SS/MDR disponibles)
#   zeta_i : dispersión estructural genuina, Var(zeta_i) = psi (igual
#            para todos)
#   Var(e_i) = psi + 1/I_i
#
# z_i = e_i / sqrt(psi_hat + 1/I_i)  — cada observación se estandariza
# contra SU PROPIA varianza esperada, no contra una escala fija para
# toda la muestra.
#
# NOTA sobre la extracción de cargas: no se puede usar
# lavInspect(ajuste,"est")$lambda con indexado de columna fijo, porque
# en un MIMIC con regresiones causales sobre los propios indicadores
# (p.ej. ln_SS_mens ~ Contrato_Parc + ...), lavaan amplía esa matriz con
# una columna por cada variable endógena, y las cargas reales quedan
# repartidas en columnas con el nombre del propio indicador. Se extraen
# por tanto con parameterEstimates(), igual que predecir_bartlett() en
# funciones_modelo.R.
calcular_residuo_heterocedastico <- function(ajuste, datos, causas,
                                             indicadores = c("ln_AT_mens",
                                                            "ln_SS_mens",
                                                            "ln_MDR"),
                                             factor_latente = "Salario_Real",
                                             metodo_rlm = "M",
                                             k_huber = 1.345,
                                             alpha = 0.05) {

  params <- lavaan::parameterEstimates(ajuste)

  lam <- params %>%
    dplyr::filter(op == "=~", lhs == factor_latente, rhs %in% indicadores) %>%
    dplyr::select(rhs, est) %>%
    tibble::deframe()
  lam <- lam[indicadores]

  th <- params %>%
    dplyr::filter(op == "~~", lhs == rhs, lhs %in% indicadores) %>%
    dplyr::select(lhs, est) %>%
    tibble::deframe()
  th <- th[indicadores]

  present <- !is.na(as.matrix(datos[, indicadores]))
  Info    <- as.numeric(present %*% (lam^2 / th))
  var_medicion <- ifelse(Info > 0, 1 / Info, NA_real_)

  eta_hat <- as.numeric(lavaan::lavPredict(ajuste, method = "Bartlett")[, factor_latente])

  ok <- is.finite(eta_hat) & is.finite(var_medicion)

  dd <- data.frame(eta_hat = eta_hat, datos[, causas, drop = FALSE])

  formula_rob <- stats::reformulate(causas, "eta_hat")

  ajuste_rob <- MASS::rlm(
    formula_rob, data = dd[ok, ],
    method = metodo_rlm,
    k      = if (metodo_rlm == "M") k_huber else NULL,
    maxit  = 200
  )

  fitted_rob <- rep(NA_real_, length(eta_hat))
  fitted_rob[ok] <- fitted(ajuste_rob)

  s2_rob  <- ajuste_rob$s^2
  psi_hat <- max(s2_rob - mean(var_medicion[ok], na.rm = TRUE), 0)

  e <- eta_hat - fitted_rob
  z <- e / sqrt(psi_hat + var_medicion)
  p <- 2 * (1 - pnorm(abs(z)))
  p_bh <- p.adjust(p, method = "BH")

  data.frame(
    eta_hat, var_medicion, fitted_rob, residuo = e, z, p_valor = p,
    p_valor_bh = p_bh,
    flag_bin   = as.integer(!is.na(p_bh) & p_bh < alpha),
    n_indicadores_validos = rowSums(present),
    psi_hat = psi_hat, s2_rob = s2_rob
  )
}


# ── Guardar gráfico ggplot ───────────────────────────────────────────────
guardar_grafico_depuracion <- function(p, carpeta, nombre, ancho = 10, alto = 6) {
  ruta <- file.path(carpeta, paste0(nombre, ".png"))
  ggsave(ruta, plot = p, width = ancho, height = alto, dpi = 300)
  cat("  Guardado:", nombre, ".png\n")
}
