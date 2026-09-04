# =============================================================================
# FUNCIONES_MODELO.R — PREPROCESAMIENTO Y AJUSTE COMPARTIDOS
# =============================================================================
#
# Funciones reutilizadas tanto por el script de producción como por el de
# evaluación, para garantizar que ambos aplican exactamente el mismo
# preprocesamiento y la misma especificación de modelo. Evita divergencias
# entre "lo que se evaluó" y "lo que se pone en producción".
#
# Requiere que config.R haya sido cargado previamente (usa sus constantes).
# =============================================================================

library(tidyverse)
library(lavaan)


# ── Validación de columnas de entrada ────────────────────────────────────────
validar_columnas <- function(data, columnas_esperadas) {
  faltantes <- setdiff(columnas_esperadas, names(data))
  if (length(faltantes) > 0) {
    stop("Faltan columnas en el fichero de entrada: ",
         paste(faltantes, collapse = ", "))
  }
  invisible(TRUE)
}



# ── Diagnóstico de pérdida de observaciones ──────────────────────────────────
# Calcula, sobre el dataset original completo (antes de cualquier filtro),
# cuántas observaciones se perderán en el preprocesamiento y por qué.
# Se llama desde estimacion.R justo después de la carga del Excel,
# antes de invocar preprocesar_datos().
#
# Grupos de pérdida:
#   A: Sin ninguna fuente salarial bruta  → solo imputación por perfil
#   B: Fuente bruta existe pero todos los ln_* caen (topes + mensualización)
#      B1: sin FACTOR_final
#      B2: solo SS anual topada
#      B3: SS + MDR ambas topadas (sin AT)
#      B4: otros casos residuales
#   C: Tienen algún ln_* válido pero HORASH = 9900
#      imputar mediana de grupo
#   OK: entran en el pipeline con normalidad
diagnostico_perdida <- function(data_salarios) {
  
  d <- data_salarios %>%
    mutate(
      Salario_AT  = ifelse(as.numeric(AT_SALARIO) / 100 == 0,
                           NA, as.numeric(AT_SALARIO) / 100),
      Salario_SS  = ifelse(as.numeric(SALARIO_SS_ANO) == 0,
                           NA, as.numeric(SALARIO_SS_ANO)),
      Salario_MDR = ifelse(as.numeric(SALARIO_SS_MDR) == 0,
                           NA, as.numeric(SALARIO_SS_MDR)),
      Tope_ANO    = as.character(TOPE_SS_ANO),
      Tope_MDR    = as.character(TOPE_SS_MDR),
      HORASH_real = case_when(
        as.numeric(HORASH) >= 9900 ~ NA_real_,
        TRUE ~ as.numeric(HORASH) / 100
      ),
      FACTOR_SS_num  = as.numeric(FACTOR_SS),
      FACTOR_EPA_num = as.numeric(FACTOR_EPA),
      FACTOR_final   = case_when(
        !is.na(FACTOR_SS_num) ~ FACTOR_SS_num,
        !is.na(FACTOR_EPA_num) &
          abs(FACTOR_EPA_num - 0.0833) < 0.001 ~ FACTOR_EPA_num,
        TRUE ~ NA_real_
      ),
      Salario_SS_Corr = ifelse(
        Tope_ANO %in% c("1. BMAX SALARIO_SS_ANO", "2. BMIN SALARIO_SS_ANO"),
        NA_real_, Salario_SS),
      Salario_MDR_Corr = ifelse(
        Tope_MDR %in% c("1. BMAX SALARIO_SS_MDR", "2. BMIN SALARIO_SS_MDR"),
        NA_real_, Salario_MDR)
    ) %>%
    mutate(
      ln_AT_raw  = ifelse(!is.na(Salario_AT) & !is.na(FACTOR_final) &
                            FACTOR_final > 0,
                          log(Salario_AT * FACTOR_final), NA_real_),
      ln_AT_mens = ifelse(!is.na(ln_AT_raw) & ln_AT_raw < 0,
                          NA_real_, ln_AT_raw),
      ln_SS_mens = ifelse(!is.na(Salario_SS_Corr) & !is.na(FACTOR_final) &
                            FACTOR_final > 0,
                          log(Salario_SS_Corr * FACTOR_final), NA_real_),
      ln_MDR     = ifelse(!is.na(Salario_MDR_Corr) & Salario_MDR_Corr > 0,
                          log(Salario_MDR_Corr), NA_real_)
    ) %>%
    mutate(
      tiene_fuente_bruta = !(is.na(Salario_AT) & is.na(Salario_SS) &
                               is.na(Salario_MDR)),
      tiene_ln           = !(is.na(ln_AT_mens) & is.na(ln_SS_mens) &
                               is.na(ln_MDR)),
      causa_perdida = case_when(
        # A: ninguna fuente salarial bruta
        !tiene_fuente_bruta ~ "A: Sin ninguna fuente salarial bruta",
        # B: fuente bruta existe pero todos los ln_* caen
        tiene_fuente_bruta & !tiene_ln & is.na(FACTOR_final)
        ~ "B1: Sin FACTOR_final (no se puede mensualizar)",
        tiene_fuente_bruta & !tiene_ln &
          !is.na(Salario_SS) & is.na(Salario_AT) & is.na(Salario_MDR) &
          Tope_ANO %in% c("1. BMAX SALARIO_SS_ANO", "2. BMIN SALARIO_SS_ANO")
        ~ "B2: Solo SS anual y esta topada",
        tiene_fuente_bruta & !tiene_ln &
          !is.na(Salario_SS) & is.na(Salario_AT) & !is.na(Salario_MDR) &
          Tope_ANO %in% c("1. BMAX SALARIO_SS_ANO", "2. BMIN SALARIO_SS_ANO") &
          Tope_MDR %in% c("1. BMAX SALARIO_SS_MDR", "2. BMIN SALARIO_SS_MDR")
        ~ "B3: SS + MDR ambas topadas (sin AT)",
        tiene_fuente_bruta & !tiene_ln
        ~ "B4: Otros (fuente bruta sin ln valido)",
        # C: tiene ln pero HORASH no disponible
        tiene_ln & is.na(HORASH_real)
        ~ "C: HORASH = 9900 (no procede en EPA)",
        # OK
        TRUE ~ "OK: Entra en pipeline"
      )
    )
  
  tabla <- d %>%
    count(causa_perdida) %>%
    mutate(Pct = round(n / nrow(d) * 100, 1)) %>%
    arrange(causa_perdida)
  
  n_ok      <- tabla$n[tabla$causa_perdida == "OK: Entra en pipeline"]
  n_perdido <- nrow(d) - n_ok
  
  cat("Total observaciones en el dataset original:", nrow(d), "\n\n")
  print(as.data.frame(tabla), row.names = FALSE)
  cat(sprintf(
    "\nResumen: %d entran en el pipeline (%.1f%%) | %d se descartan (%.1f%%)\n",
    n_ok,      n_ok      / nrow(d) * 100,
    n_perdido, n_perdido / nrow(d) * 100
  ))
  cat("\nNota: las observaciones del grupo C (HORASH = 9900) tienen fuente\n")
  cat("salarial valida. Su tratamiento (imputacion de HORASH por mediana de\n")
  cat("grupo Ocupacion x Jornada, o MIMIC alternativo sin Horas_Z)\n")
  
  invisible(tabla)
}


# ── Preprocesamiento completo ────────────────────────────────────────────────
# Devuelve una lista con dos datasets: 'base' (todas las observaciones con
# al menos un indicador) y 'ss' (subconjunto sin NA en Contrato_Parc,
# necesario para el modelo con causas en indicadores SS).
preprocesar_datos <- function(data_salarios, contratos_parciales) {
  
  data_salarios_modelo <- data_salarios %>%
    mutate(idx = row_number()) %>%
    mutate(
      Salario_AT  = as.numeric(AT_SALARIO) / 100,
      Salario_SS  = as.numeric(SALARIO_SS_ANO),
      Salario_MDR = as.numeric(SALARIO_SS_MDR),
      Edad        = as.numeric(EDAD1),
      Sexo_Orig   = as.numeric(SEXO1),
      Meses_Emp   = as.numeric(DCOM),
      Jornada     = as.numeric(PARCO1),
      Tope_ANO    = as.character(TOPE_SS_ANO),
      Tope_MDR    = as.character(TOPE_SS_MDR),
      HORASH_real = case_when(
        as.numeric(HORASH) >= 9900 ~ NA_real_,
        TRUE ~ as.numeric(HORASH) / 100
      )
    ) %>%
    mutate(
      Salario_AT  = ifelse(Salario_AT  == 0, NA, Salario_AT),
      Salario_SS  = ifelse(Salario_SS  == 0, NA, Salario_SS),
      Salario_MDR = ifelse(Salario_MDR == 0, NA, Salario_MDR)
    ) %>%
    mutate(
      Salario_SS_Corr = ifelse(
        Tope_ANO %in% c("1. BMAX SALARIO_SS_ANO",
                        "2. BMIN SALARIO_SS_ANO"),
        NA_real_, Salario_SS),
      Salario_MDR_Corr = ifelse(
        Tope_MDR %in% c("1. BMAX SALARIO_SS_MDR",
                        "2. BMIN SALARIO_SS_MDR"),
        NA_real_, Salario_MDR)
    ) %>%
    filter(!(is.na(Salario_AT) & is.na(Salario_SS_Corr) &
               is.na(Salario_MDR_Corr))) %>%
    filter(!is.na(Edad), !is.na(Meses_Emp),
           !is.na(Jornada), !is.na(Sexo_Orig)) %>%
    mutate(
      FACTOR_SS_num  = as.numeric(FACTOR_SS),
      FACTOR_EPA_num = as.numeric(FACTOR_EPA),
      FACTOR_final   = case_when(
        !is.na(FACTOR_SS_num)  ~ FACTOR_SS_num,
        !is.na(FACTOR_EPA_num) &
          abs(FACTOR_EPA_num - 0.0833) < 0.001 ~ FACTOR_EPA_num,
        TRUE ~ NA_real_
      )
    ) %>%
    mutate(
      ln_AT_mens = ifelse(!is.na(Salario_AT) & !is.na(FACTOR_final) &
                            FACTOR_final > 0,
                          log(Salario_AT * FACTOR_final), NA_real_),
      ln_SS_mens = ifelse(!is.na(Salario_SS_Corr) & !is.na(FACTOR_final) &
                            FACTOR_final > 0,
                          log(Salario_SS_Corr * FACTOR_final), NA_real_),
      ln_MDR     = ifelse(!is.na(Salario_MDR_Corr) & Salario_MDR_Corr > 0,
                          log(Salario_MDR_Corr), NA_real_)
    ) %>%
    mutate(
      ln_AT_mens = ifelse(!is.na(ln_AT_mens) &
                            ln_AT_mens < 0, NA_real_, ln_AT_mens)
    ) %>%
    filter(!(is.na(ln_AT_mens) & is.na(ln_SS_mens) & is.na(ln_MDR))) %>%
    mutate(
      Mujer_Dummy      = ifelse(Sexo_Orig == 6, 1, 0),
      Jornada_Completa = ifelse(Jornada == 1, 1, 0),
      Edad_Z           = as.numeric(scale(Edad)),
      Edad_Z2          = as.numeric(scale(Edad^2)),
      Exp_Z            = as.numeric(scale(Meses_Emp)),
      Exp_Z2           = Exp_Z^2,
      Extranjero       = ifelse(as.numeric(NAC1) == 3, 1, 0),
      Horas_Z          = as.numeric(scale(HORASH_real)),
      
      NFORMA_num   = as.numeric(NFORMA),
      NFORMA_grupo = floor(NFORMA_num / 10),
      Educacion    = case_when(
        NFORMA_grupo <= 2          ~ "Basica",
        NFORMA_grupo == 3          ~ "Secundaria",
        NFORMA_grupo %in% c(4, 5)  ~ "FP_Bach",
        NFORMA_grupo >= 6          ~ "Superior"
      ),
      Educ_Secundaria = ifelse(Educacion == "Secundaria", 1, 0),
      Educ_FP_Bach    = ifelse(Educacion == "FP_Bach",    1, 0),
      Educ_Superior   = ifelse(Educacion == "Superior",   1, 0),
      
      OCUP_cod   = as.character(OCUP),
      OCUP_grupo = as.numeric(substr(OCUP_cod, 1, 1)),
      Ocupacion  = case_when(
        OCUP_grupo %in% c(0, 1) ~ "Directivos",
        OCUP_grupo == 2         ~ "Tecnicos_Sup",
        OCUP_grupo == 3         ~ "Tecnicos_Ap",
        OCUP_grupo == 4         ~ "Administrativos",
        OCUP_grupo == 5         ~ "Servicios",
        OCUP_grupo == 6         ~ "Agricultura",
        OCUP_grupo == 7         ~ "Artesanos",
        OCUP_grupo == 8         ~ "Operadores",
        OCUP_grupo == 9         ~ "Elementales"
      ),
      Ocup_Directivos  = ifelse(Ocupacion == "Directivos",      1, 0),
      Ocup_TecSup      = ifelse(Ocupacion == "Tecnicos_Sup",    1, 0),
      Ocup_TecAp       = ifelse(Ocupacion == "Tecnicos_Ap",     1, 0),
      Ocup_Admin       = ifelse(Ocupacion == "Administrativos", 1, 0),
      Ocup_Servicios   = ifelse(Ocupacion == "Servicios",       1, 0),
      Ocup_Agricultura = ifelse(Ocupacion == "Agricultura",     1, 0),
      Ocup_Artesanos   = ifelse(Ocupacion == "Artesanos",       1, 0),
      Ocup_Operadores  = ifelse(Ocupacion == "Operadores",      1, 0),
      
      ACT_2dig = as.numeric(substr(as.character(ACT), 1, 2)),
      Sector   = case_when(
        ACT_2dig >= 1  & ACT_2dig <= 3  ~ "Agricultura",
        ACT_2dig >= 5  & ACT_2dig <= 39 ~ "Industria",
        ACT_2dig >= 41 & ACT_2dig <= 43 ~ "Construccion",
        ACT_2dig >= 45 & ACT_2dig <= 99 ~ "Servicios",
        TRUE ~ NA_character_
      ),
      Sector_Agri  = ifelse(Sector == "Agricultura",  1, 0),
      Sector_Ind   = ifelse(Sector == "Industria",    1, 0),
      Sector_Const = ifelse(Sector == "Construccion", 1, 0)
    ) %>%
    mutate(
      Grupo_cot     = as.numeric(NGRUP_MDR),
      Grupo_cot_rec = case_when(
        is.na(Grupo_cot)                   ~ NA_real_,
        Grupo_cot %in% c(1, 2, 3)         ~ 1,
        Grupo_cot %in% c(4, 5, 6, 7)      ~ 2,
        Grupo_cot %in% c(8, 9, 10, 11, 0) ~ 3
      ),
      GC_alta = ifelse(!is.na(Grupo_cot_rec) & Grupo_cot_rec == 1, 1, 0),
      GC_baja = ifelse(!is.na(Grupo_cot_rec) & Grupo_cot_rec == 3, 1, 0)
    ) %>%
    left_join(
      data_salarios %>%
        mutate(
          idx           = row_number(),
          CNTRATO_c     = as.character(cntrato),
          Contrato_Parc = ifelse(CNTRATO_c %in% contratos_parciales, 1, 0)
        ) %>%
        select(idx, Contrato_Parc),
      by = "idx"
    ) %>%
    filter(!is.na(HORASH_real))
  
  data_salarios_modelo_ss <- data_salarios_modelo %>%
    filter(!is.na(Contrato_Parc))
  
  list(base = data_salarios_modelo, ss = data_salarios_modelo_ss)
}




# ── Ajuste del modelo MIMIC ───────────────────────────────────────────────
ajustar_mimic <- function(data, medida, causas, causas_ind_ss = NULL) {
  
  formula <- if (is.null(causas_ind_ss)) {
    paste0(medida, causas)
  } else {
    paste0(medida, causas, causas_ind_ss)
  }
  
  ajuste <- tryCatch(
    sem(
      model         = formula,
      data          = data,
      missing       = "fiml",
      estimator     = "MLR",
      bounds        = TRUE,
      meanstructure = TRUE
    ),
    error = function(e) {
      stop("Error al ajustar el modelo SEM: ", conditionMessage(e))
    }
  )
  
  if (!lavInspect(ajuste, "converged")) {
    stop("El modelo MIMIC no convergi\u00f3. Revisar los datos de entrada ",
         "(posibles cambios en la estructura de missings o en la ",
         "distribuci\u00f3n de alguna variable causal).")
  }
  
  ajuste
}


# ── Predicción puntual (estimador de Bartlett) ───────────────────────────────
predecir_bartlett <- function(ajuste, df, devolver_score = FALSE) {
  intercepto <- parameterEstimates(ajuste) %>%
    filter(lhs == "ln_AT_mens", op == "~1") %>%
    pull(est)
  scores <- as.data.frame(lavPredict(ajuste, type = "lv",
                                     method = "Bartlett"))
  score_log <- intercepto + scores[, 1]
  
  if (devolver_score) {
    # Devuelve el score en escala log (mu_i), SIN exponenciar. Permite
    # aplicar después la corrección de sigma: exp(mu_i + sigma^2/2) para
    # la media de la lognormal, en lugar de exp(mu_i) que es la mediana.
    return(score_log)
  }
  
  exp(score_log)
}


# ── Parámetros de medición ───────────────────────────────────────────────
# Cargas factoriales, varianzas residuales, covarianzas entre indicadores,
# diagnóstico de casos Heywood e índices de ajuste global. Solo para uso
# en evaluación/auditoría del modelo (no depende de STP3).
imprimir_medicion <- function(ajuste, nombre) {
  cat(">>>", nombre, "<<<\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  params <- parameterEstimates(ajuste, standardized = TRUE) %>%
    filter(
      (op == "=~") |
        (op == "~~" & lhs == rhs &
           lhs %in% c("ln_AT_mens", "ln_SS_mens",
                      "ln_MDR", "Salario_Real")) |
        (op == "~~" & lhs != rhs &
           lhs %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR"))
    ) %>%
    mutate(
      Parametro = case_when(
        op == "=~"               ~ paste0("Carga (\u03bb) ", rhs),
        op == "~~" & lhs == rhs  ~ paste0("Varianza residual (\u03b8) ", lhs),
        op == "~~" & lhs != rhs  ~ paste0("Covarianza ", lhs, "-", rhs)
      ),
      across(where(is.numeric), ~round(., 4))
    ) %>%
    select(Parametro, est, se, pvalue, std.all)
  
  print(params, row.names = FALSE)
  
  vars_neg <- parameterEstimates(ajuste) %>%
    filter(op == "~~", lhs == rhs,
           lhs %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR"),
           est < -0.001)
  cat("Heywood:", ifelse(nrow(vars_neg) > 0, "S\u00cd", "NO"), "\n")
  
  cat("\n\u00cdndices de ajuste:\n")
  cat("  GL:        ", fitMeasures(ajuste, "df"), "\n")
  cat("  AIC:       ", round(AIC(ajuste), 2), "\n")
  cat("  BIC:       ", round(BIC(ajuste), 2), "\n")
  cat("  RMSEA rob: ", round(fitMeasures(ajuste, "rmsea.robust"), 4), "\n")
  cat("  CFI rob:   ", round(fitMeasures(ajuste, "cfi.robust"), 4), "\n")
  cat("  SRMR:      ", round(fitMeasures(ajuste, "srmr"), 4), "\n")
  cat("\n")
}


# ── Parámetros estructurales ─────────────────────────────────────────────
# Efecto de cada causa sobre el factor latente Salario_Real. Solo para uso
# en evaluación/auditoría del modelo.
imprimir_estructural <- function(ajuste, nombre) {
  cat(">>>", nombre, "<<<\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  params <- parameterEstimates(ajuste, standardized = TRUE) %>%
    filter(op == "~", lhs == "Salario_Real") %>%
    mutate(
      Variable = case_when(
        rhs == "Jornada_Completa" ~ "Jornada completa",
        rhs == "Mujer_Dummy"      ~ "Mujer",
        rhs == "Edad_Z"           ~ "Edad (Z)",
        rhs == "Edad_Z2"          ~ "Edad\u00b2 (Z)",
        rhs == "Exp_Z"            ~ "Experiencia (Z)",
        rhs == "Exp_Z2"           ~ "Experiencia\u00b2 (Z)",
        rhs == "Extranjero"       ~ "Extranjero",
        rhs == "Horas_Z"          ~ "Horas habituales (Z)",
        rhs == "Educ_Secundaria"  ~ "Educ: Secundaria",
        rhs == "Educ_FP_Bach"     ~ "Educ: FP/Bach",
        rhs == "Educ_Superior"    ~ "Educ: Superior",
        rhs == "Ocup_Directivos"  ~ "Ocup: Directivos",
        rhs == "Ocup_TecSup"      ~ "Ocup: T\u00e9c. Sup.",
        rhs == "Ocup_TecAp"       ~ "Ocup: T\u00e9c. Apoyo",
        rhs == "Ocup_Admin"       ~ "Ocup: Admin.",
        rhs == "Ocup_Servicios"   ~ "Ocup: Servicios",
        rhs == "Ocup_Agricultura" ~ "Ocup: Agricultura",
        rhs == "Ocup_Artesanos"   ~ "Ocup: Artesanos",
        rhs == "Ocup_Operadores"  ~ "Ocup: Operadores",
        rhs == "Sector_Agri"      ~ "Sector: Agricultura",
        rhs == "Sector_Ind"       ~ "Sector: Industria",
        rhs == "Sector_Const"     ~ "Sector: Construcci\u00f3n",
        TRUE                      ~ rhs
      ),
      Sig = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01  ~ "**",
        pvalue < 0.05  ~ "*",
        pvalue < 0.10  ~ ".",
        TRUE           ~ ""
      ),
      across(where(is.numeric), ~round(., 4))
    ) %>%
    select(Variable, est, se, z, pvalue, std.all, Sig) %>%
    arrange(desc(abs(std.all)))
  
  print(params, row.names = FALSE)
  cat("\n")
}


# ── Causas en indicadores SS (solo M1_ind) ───────────────────────────────
imprimir_causas_ind_ss <- function(ajuste) {
  cat(">>> M1_ind — Causas en indicadores SS <<<\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  parameterEstimates(ajuste, standardized = TRUE) %>%
    filter(op == "~",
           lhs %in% c("ln_SS_mens", "ln_MDR"),
           rhs %in% c("Contrato_Parc", "GC_alta", "GC_baja")) %>%
    mutate(
      Indicador = case_when(
        lhs == "ln_SS_mens" ~ "SS anual",
        lhs == "ln_MDR"     ~ "MDR mensual"
      ),
      Variable  = case_when(
        rhs == "Contrato_Parc" ~ "Contrato parcial",
        rhs == "GC_alta"       ~ "Grupo cot. alto (1-3)",
        rhs == "GC_baja"       ~ "Grupo cot. bajo (8-11)"
      ),
      Sig = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01  ~ "**",
        pvalue < 0.05  ~ "*",
        pvalue < 0.10  ~ ".",
        TRUE           ~ ""
      ),
      across(where(is.numeric), ~round(., 4))
    ) %>%
    select(Indicador, Variable, est, se, z, pvalue, std.all, Sig) %>%
    arrange(Indicador, desc(abs(std.all))) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  
  cat("\n")
}


# ── Auditoría: residuos e índices de modificación ────────────────────────
# Solo para uso en evaluación/auditoría del modelo.
imprimir_auditoria <- function(ajuste, nombre) {
  cat(">>>", nombre, "<<<\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  cat("\nRes\u00edduos |z| > 1.96 (indicadores ~ causas):\n")
  tryCatch({
    res_z <- lavResiduals(ajuste, type = "cor.bentler")$cov.z
    res_df <- as.data.frame(res_z) %>%
      rownames_to_column("var1") %>%
      pivot_longer(-var1, names_to = "var2", values_to = "z") %>%
      filter(
        var1 %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR"),
        !var2 %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR"),
        abs(z) > 1.96
      ) %>%
      arrange(desc(abs(z)))
    
    if (nrow(res_df) > 0) {
      print(head(res_df, 20), row.names = FALSE)
      cat("Total residuos problem\u00e1ticos:", nrow(res_df), "\n")
    } else {
      cat("Ning\u00fan residuo supera 1.96\n")
    }
  }, error = function(e) cat("No disponible:", e$message, "\n"))
  
  cat("\nRes\u00edduos entre indicadores (cor.bentler):\n")
  tryCatch({
    lavResiduals(ajuste, type = "cor.bentler")$cov.z[
      c("ln_AT_mens", "ln_SS_mens", "ln_MDR"),
      c("ln_AT_mens", "ln_SS_mens", "ln_MDR")
    ] %>% round(3) %>% print()
  }, error = function(e) cat("No disponible\n"))
  
  cat("\n\u00cdndices de modificaci\u00f3n — covarianzas entre indicadores:\n")
  tryCatch({
    mi_cov <- modindices(ajuste, sort. = TRUE, minimum.value = 0) %>%
      filter(op == "~~", lhs != rhs,
             lhs %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR") |
               rhs %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR")) %>%
      select(lhs, op, rhs, mi, epc) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
    if (nrow(mi_cov) > 0) print(mi_cov, row.names = FALSE)
    else cat("Sin \u00edndices disponibles\n")
  }, error = function(e) cat("No disponible\n"))
  
  cat("\n\u00cdndices de modificaci\u00f3n — regresiones (top 10):\n")
  tryCatch({
    mi_reg <- modindices(ajuste, sort. = TRUE, minimum.value = 0) %>%
      filter(op == "~",
             lhs %in% c("ln_AT_mens", "ln_SS_mens", "ln_MDR")) %>%
      head(10) %>%
      select(lhs, op, rhs, mi, epc) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
    if (nrow(mi_reg) > 0) print(mi_reg, row.names = FALSE)
    else cat("Sin \u00edndices disponibles\n")
  }, error = function(e) cat("No disponible\n"))
  
  cat("\n")
}


# ── Guardar gráfico ggplot ────────────────────────────────────────────────
# Construye la ruta completa a partir de la carpeta y el nombre, y guarda
# el gráfico a 300 ppp. Usada tanto en evaluación como en cualquier futuro
# script (Isolation Forest, etc.) que necesite guardar gráficos igual.
guardar_grafico <- function(p, carpeta, nombre, ancho = 10, alto = 6) {
  ruta <- file.path(carpeta, paste0(nombre, ".png"))
  ggsave(ruta, plot = p, width = ancho, height = alto, dpi = 300)
  cat("  Guardado:", nombre, ".png\n")
}


# ── Preprocesamiento SOLO de causas (sin filtrar por indicadores) ────────────
# Variante de preprocesar_datos() usada por imputacion.R para recuperar el
# PERFIL sociodemográfico de TODAS las observaciones del dataset original,
# incluidas las que no tienen ningún indicador administrativo (las ~1.629
# que preprocesar_datos() descarta). Construye exactamente las mismas
# variables de causa que preprocesar_datos() (Edad_Z, Educ_*, Ocup_*,
# Sector_*, etc.), pero SIN aplicar los filtros que eliminan por falta de
# indicador o de HORASH, ya que esas observaciones no van a alimentar el
# MIMIC sino un modelo de imputación basado solo en el perfil.
#
# No calcula los indicadores logarítmicos (ln_AT_mens, etc.) porque no se
# usan aquí. Devuelve un único data.frame con idx y todas las causas de
# perfil, para el conjunto completo de observaciones que tengan el perfil
# demográfico básico disponible.
preprocesar_datos_causas <- function(data_salarios, contratos_parciales) {
  
  data_salarios %>%
    mutate(idx = row_number()) %>%
    mutate(
      Edad        = as.numeric(EDAD1),
      Sexo_Orig   = as.numeric(SEXO1),
      Meses_Emp   = as.numeric(DCOM),
      Jornada     = as.numeric(PARCO1),
      HORASH_real = case_when(
        as.numeric(HORASH) >= 9900 ~ NA_real_,
        TRUE ~ as.numeric(HORASH) / 100
      )
    ) %>%
    mutate(
      Mujer_Dummy      = ifelse(Sexo_Orig == 6, 1, 0),
      Jornada_Completa = ifelse(Jornada == 1, 1, 0),
      Edad_Z           = as.numeric(scale(Edad)),
      Edad_Z2          = as.numeric(scale(Edad^2)),
      Exp_Z            = as.numeric(scale(Meses_Emp)),
      Exp_Z2           = Exp_Z^2,
      Extranjero       = ifelse(as.numeric(NAC1) == 3, 1, 0),
      Horas_Z          = as.numeric(scale(HORASH_real)),
      
      NFORMA_num   = as.numeric(NFORMA),
      NFORMA_grupo = floor(NFORMA_num / 10),
      Educacion    = case_when(
        NFORMA_grupo <= 2          ~ "Basica",
        NFORMA_grupo == 3          ~ "Secundaria",
        NFORMA_grupo %in% c(4, 5)  ~ "FP_Bach",
        NFORMA_grupo >= 6          ~ "Superior"
      ),
      Educ_Secundaria = ifelse(Educacion == "Secundaria", 1, 0),
      Educ_FP_Bach    = ifelse(Educacion == "FP_Bach",    1, 0),
      Educ_Superior   = ifelse(Educacion == "Superior",   1, 0),
      
      OCUP_cod   = as.character(OCUP),
      OCUP_grupo = as.numeric(substr(OCUP_cod, 1, 1)),
      Ocupacion  = case_when(
        OCUP_grupo %in% c(0, 1) ~ "Directivos",
        OCUP_grupo == 2         ~ "Tecnicos_Sup",
        OCUP_grupo == 3         ~ "Tecnicos_Ap",
        OCUP_grupo == 4         ~ "Administrativos",
        OCUP_grupo == 5         ~ "Servicios",
        OCUP_grupo == 6         ~ "Agricultura",
        OCUP_grupo == 7         ~ "Artesanos",
        OCUP_grupo == 8         ~ "Operadores",
        OCUP_grupo == 9         ~ "Elementales"
      ),
      Ocup_Directivos  = ifelse(Ocupacion == "Directivos",      1, 0),
      Ocup_TecSup      = ifelse(Ocupacion == "Tecnicos_Sup",    1, 0),
      Ocup_TecAp       = ifelse(Ocupacion == "Tecnicos_Ap",     1, 0),
      Ocup_Admin       = ifelse(Ocupacion == "Administrativos", 1, 0),
      Ocup_Servicios   = ifelse(Ocupacion == "Servicios",       1, 0),
      Ocup_Agricultura = ifelse(Ocupacion == "Agricultura",     1, 0),
      Ocup_Artesanos   = ifelse(Ocupacion == "Artesanos",       1, 0),
      Ocup_Operadores  = ifelse(Ocupacion == "Operadores",      1, 0),
      
      ACT_2dig = as.numeric(substr(as.character(ACT), 1, 2)),
      Sector   = case_when(
        ACT_2dig >= 1  & ACT_2dig <= 3  ~ "Agricultura",
        ACT_2dig >= 5  & ACT_2dig <= 39 ~ "Industria",
        ACT_2dig >= 41 & ACT_2dig <= 43 ~ "Construccion",
        ACT_2dig >= 45 & ACT_2dig <= 99 ~ "Servicios",
        TRUE ~ NA_character_
      ),
      Sector_Agri  = ifelse(Sector == "Agricultura",  1, 0),
      Sector_Ind   = ifelse(Sector == "Industria",    1, 0),
      Sector_Const = ifelse(Sector == "Construccion", 1, 0)
    ) %>%
    mutate(
      Grupo_cot     = as.numeric(NGRUP_MDR),
      Grupo_cot_rec = case_when(
        is.na(Grupo_cot)                   ~ NA_real_,
        Grupo_cot %in% c(1, 2, 3)         ~ 1,
        Grupo_cot %in% c(4, 5, 6, 7)      ~ 2,
        Grupo_cot %in% c(8, 9, 10, 11, 0) ~ 3
      ),
      GC_alta = ifelse(!is.na(Grupo_cot_rec) & Grupo_cot_rec == 1, 1, 0),
      GC_baja = ifelse(!is.na(Grupo_cot_rec) & Grupo_cot_rec == 3, 1, 0)
    ) %>%
    left_join(
      data_salarios %>%
        mutate(
          idx           = row_number(),
          CNTRATO_c     = as.character(cntrato),
          Contrato_Parc = ifelse(CNTRATO_c %in% contratos_parciales, 1, 0)
        ) %>%
        select(idx, Contrato_Parc),
      by = "idx"
    )
}


# ── Imputación de HORASH por mediana de grupo ────────────────────────────────
# Recibe el dataset de entrada (con las columnas HORASH, OCUP, PARCO1) y
# devuelve el MISMO dataset con la columna HORASH sobrescrita: los valores
# 9900 ("no procede" en la EPA) se sustituyen por la mediana de horas de su
# grupo Ocupación × Jornada (en la escala original × 100). El resto de
# valores de HORASH quedan intactos.
#
# Se usa de forma idéntica en estimacion.R, depuracion.R e imputacion.R para
# garantizar que los tres scripts trabajan sobre exactamente el mismo
# conjunto de observaciones (las que recuperan las ~477 con HORASH = 9900).
#
# Imputación DETERMINISTA (mediana de grupo): más transparente y auditable
# que un método estocástico para una variable auxiliar de baja dispersión
# intragrupo. Si algún grupo no tuviera mediana (todos NA), se usa la
# mediana global como fallback.
imputar_horash <- function(data_salarios) {
  
  data_horas <- data_salarios %>%
    mutate(
      idx = row_number(),
      HORASH_real = case_when(
        as.numeric(HORASH) >= 9900 ~ NA_real_,
        TRUE ~ as.numeric(HORASH) / 100
      ),
      OCUP_grupo_h    = as.numeric(substr(as.character(OCUP), 1, 1)),
      Jornada_Comp_h  = ifelse(as.numeric(PARCO1) == 1, 1, 0),
      HORASH_faltante = is.na(HORASH_real)
    ) %>%
    group_by(OCUP_grupo_h, Jornada_Comp_h) %>%
    mutate(HORASH_imputada = ifelse(HORASH_faltante,
                                    median(HORASH_real, na.rm = TRUE),
                                    HORASH_real)) %>%
    ungroup()
  
  mediana_global_horas <- median(data_horas$HORASH_real, na.rm = TRUE)
  data_horas <- data_horas %>%
    mutate(HORASH_imputada = ifelse(is.na(HORASH_imputada),
                                    mediana_global_horas, HORASH_imputada))
  
  n_imputadas <- sum(data_horas$HORASH_faltante)
  cat("HORASH imputadas (código 9900):", n_imputadas,
      "| mediana de grupo Ocupación × Jornada",
      sprintf("(fallback global: %.1f h)\n", mediana_global_horas))
  
  # Sobrescribir HORASH (escala original × 100) solo en los faltantes.
  # Se usa una clave temporal (.idx_tmp) para no interferir con una posible
  # columna idx preexistente en el dataset de entrada.
  data_salarios %>%
    mutate(.idx_tmp = row_number()) %>%
    left_join(data_horas %>% select(idx, HORASH_imputada, HORASH_faltante) %>%
                rename(.idx_tmp = idx),
              by = ".idx_tmp") %>%
    mutate(HORASH = ifelse(HORASH_faltante,
                           round(HORASH_imputada * 100), HORASH)) %>%
    select(-.idx_tmp, -HORASH_imputada, -HORASH_faltante)
}