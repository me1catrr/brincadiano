# Métricas de luz con LightLogR sobre la serie por minuto que exporta Brincadiano.
#
# Lo llama scripts/94_validar_luz_lightlogr.jl. LightLogR no lee ActiGraph: la serie
# ya viene reconstruida (script 09). No lanzarlo a mano sin la autorización del 94.
#
#   Rscript lightlogr_metricas.R <entrada.csv> <salida_dias.csv> <salida_global.csv>
#   (directorio de trabajo: scripts/R, para el renv)
#
# Qué hace
#   Por día: minutos con luz ≥ 100 y ≥ 1000, y primera, última y hora media con luz
#   ≥ 100, en horas. En el registro: IS e IV de luz_vigilia. El 94 borra los dos CSV
#   de salida en cuanto los ha leído; la comparación queda en
#   outputs/validation/tables/luz_lightlogr.{csv,tex}.
#
# Depende de
#   El renv de scripts/R/ (LightLogR 0.10.6) y del CSV del 94.
#
# Lo usa
#   Solo el 94. No tiene relación con ggir_benchmark.R ni con ggir_sueno.R.
#
# Entrada
#   CSV con Datetime (hora local, ISO; se lee como reloj, sin convertir zona),
#   luz_valida (bruta; NA si el minuto no es «válido») y luz_vigilia (bruta; NA si
#   el minuto es sueño, retirada o sin dato).
#
# Salida
#   <salida_dias.csv>     dia_desde_inicio, min_100, min_1000, primera/ultima/media_100_h
#   <salida_global.csv>   IS, IV
#   El día de calendario no se escribe: solo el entero desde el primer día del archivo.
#
# Consola: «LightLogR: listo».

main <- function() {
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(LightLogR); library(dplyr)})
d <- read.csv(args[1])
d$Datetime <- as.POSIXct(d$Datetime, tz = "UTC", format = "%Y-%m-%dT%H:%M:%S")  # reloj local tal cual
d$dia <- as.Date(d$Datetime, tz = "UTC")
h <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else as.numeric(x) / 3600
res <- d %>% group_by(dia) %>% group_modify(function(x, k) {
  tt <- hms::as_hms(x$Datetime)
  tm <- timing_above_threshold(x$luz_valida, tt, threshold = 100, na.rm = TRUE)
  tibble::tibble(
    min_100 = as.numeric(duration_above_threshold(x$luz_valida, x$Datetime, threshold = 100, na.rm = TRUE)) / 60,
    min_1000 = as.numeric(duration_above_threshold(x$luz_valida, x$Datetime, threshold = 1000, na.rm = TRUE)) / 60,
    primera_100_h = h(tm$first), ultima_100_h = h(tm$last), media_100_h = h(tm$mean))
}) %>% ungroup()
res$dia_desde_inicio <- as.integer(res$dia - min(d$dia)); res$dia <- NULL
write.csv(res, args[2], row.names = FALSE)
g <- data.frame(
  IS = suppressWarnings(interdaily_stability(d$luz_vigilia, d$Datetime, na.rm = TRUE)),
  IV = suppressWarnings(intradaily_variability(d$luz_vigilia, d$Datetime, na.rm = TRUE)))
write.csv(g, args[3], row.names = FALSE)
cat("LightLogR: listo\n")
}
main()
