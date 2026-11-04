# Benchmark de GGIR para sueño (partes 1–4) y, opcionalmente, intensidad (parte 5).
#
# Lo llaman dos scripts, con salidas distintas:
#   93  Rscript ggir_sueno.R <dir> <archivos.gt3x...>
#       y otra vez con --eje=2. Directorios: outputs/validation/ggir_sueno/ y ggir_sueno_ejeY/.
#   96  Rscript ggir_sueno.R --eje=2 --parte5 <dir> <archivos.gt3x...>
#       Directorio: outputs/validation/ggir_parte5/.
# El directorio de trabajo tiene que ser scripts/R (renv). No lanzarlo a mano: GGIR lee
# el info.txt completo. La autorización la exigen el 93 y el 96 (--confirmo).
#
#   Rscript ggir_sueno.R [--eje=1|2|3] [--parte5] <directorio_salida> <archivo1.gt3x> [...]
#
# Qué hace
#   Misma base que ggir_benchmark.R, más sensor.location = "hip" (GGIR 3.3.9 activa los
#   tres ángulos, HASPT.algo = "HorAngle" y sleepwindowType = "TimeInBed"). Sin --eje,
#   GGIR elige el eje longitudinal. Con --parte5 ejecuta también la parte 5: ventanas de
#   medianoche a medianoche y umbrales 63,3 / 142,6 / 464,6 mg. El temporal se borra al
#   salir. La fecha de cada noche pasa a «días desde el primer día del registro» antes
#   de escribirse.
#
# Depende de
#   El renv de scripts/R/ (GGIR 3.3.9). No lee Arrow ni tablas de Julia.
#
# Lo usan
#   93 (sin --parte5, dos ejes) y 96 (--eje=2 --parte5). No comparten directorio de salida.
#
# Entrada
#   Uno o más .gt3x. Opciones antes de las rutas: --eje=1|2|3 y --parte5.
#
# Salida (no se versiona; la escribe quien llama, en el directorio que le pase)
#   ggir_<ID>_angulos5s.csv        timestamp, anglex, angley, anglez
#   ggir_part4_noches.csv          inicio, fin, SPT, sueño, WASO, SRI, eficiencia, eje
#   ggir_part5_dias.csv            solo con --parte5: minutos IN/LIG/MOD/VIG y tiempo despierto
#
# Consola: «GGIR sueño: listo», el formato de calendar_date con los dígitos enmascarados
# y los nombres de columna de la parte 4 (y de la 5 si se pidió). No imprime valores.

main <- function() {
args <- commandArgs(trailingOnly = TRUE)
eje <- NULL; parte5 <- FALSE
while (length(args) > 0 && grepl("^--", args[1])) {
  if (grepl("^--eje=", args[1])) eje <- as.integer(sub("^--eje=", "", args[1]))
  if (args[1] == "--parte5") parte5 <- TRUE
  args <- args[-1]
}
outdir <- args[1]
files <- normalizePath(args[-1])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tmp <- file.path(tempdir(), paste0("ggir_sueno_", Sys.getpid()))
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

suppressMessages(library(GGIR))
invisible(capture.output(suppressWarnings(suppressMessages(GGIR(
  mode = if (parte5) 1:5 else 1:4, datadir = files, outputdir = tmp, studyname = "brincadiano",
  overwrite = TRUE, do.report = if (parte5) c(2, 4, 5) else c(2, 4), idloc = 2, print.filename = FALSE,
  storefolderstructure = FALSE, do.parallel = FALSE, verbose = FALSE,
  desiredtz = "Europe/Madrid", windowsizes = c(5, 900, 3600),
  do.enmo = TRUE, do.mad = TRUE, do.cal = TRUE, nonwear_approach = "2023",
  sensor.location = "hip", do.anglex = TRUE, do.angley = TRUE, do.anglez = TRUE,
  includedaycrit = 10, includenightcrit = 10, epochvalues2csv = FALSE,
  visualreport = FALSE, outliers.only = FALSE, do.visual = FALSE, longitudinal_axis = eje,
  threshold.lig = 63.3, threshold.mod = 142.6, threshold.vig = 464.6, timewindow = "MM"
)))))

base <- list.files(tmp, pattern = "^output_", full.names = TRUE)[1]
id_de <- function(f) sub("(\\.gt3x)?(\\.RData)?$", "", sub("^meta_", "", basename(f)))

# 1. Ángulos por epoch de 5 s (parte 1) y primer día de cada registro.
primer_dia <- list()
for (f in list.files(file.path(base, "meta", "basic"), full.names = TRUE)) {
  e <- new.env(); load(f, envir = e); M <- e$M
  id <- id_de(f)
  ms <- M$metashort[, intersect(c("timestamp", "anglex", "angley", "anglez"), names(M$metashort)), drop = FALSE]
  write.csv(ms, file.path(outdir, paste0("ggir_", id, "_angulos5s.csv")), row.names = FALSE)
  primer_dia[[id]] <- as.Date(substr(as.character(ms$timestamp[1]), 1, 10))
}

# 2. Resumen por noche (parte 4, sin limpiar): columnas de la lista blanca.
p <- file.path(base, "results", "QC", "part4_nightsummary_sleep_full.csv")
if (!file.exists(p)) p <- file.path(base, "results", "part4_nightsummary_sleep_cleaned.csv")
d <- read.csv(p, check.names = FALSE)
cat("columnas de la parte 4:", paste(names(d), collapse = ", "), "\n")
fecha <- as.character(d$calendar_date)
cat("formato de calendar_date:", gsub("[0-9]", "9", fecha[1]), "\n")
fx <- as.Date(fecha, tryFormats = c("%Y-%m-%d", "%d/%m/%Y"), optional = TRUE)
d$id <- id_de(d$filename)
d$dia_desde_inicio <- as.integer(fx - do.call(c, primer_dia[d$id]))
keep <- c("id", "night", "dia_desde_inicio", "weekday", "sleeponset", "wakeup", "SptDuration",
          "SleepDurationInSpt", "WASO", "guider", "guider_inbedStart", "guider_inbedEnd", "guider_inbedDuration",
          "SleepRegularityIndex1", "SriFractionValid", "sleepefficiency", "sleeplatency", "number_of_awakenings", "daysleeper", "fraction_night_invalid",
          "cleaningcode", "sleeplog_used", "acc_available", "longitudinal_axis", "number_sib_sleepperiod")
d <- d[, intersect(keep, names(d)), drop = FALSE]
write.csv(d, file.path(outdir, "ggir_part4_noches.csv"), row.names = FALSE)

# 3. Parte 5 (opcional): minutos por intensidad en el tiempo despierto de cada día.
if (parte5) {
  p5 <- list.files(file.path(base, "results"), pattern = "^part5_daysummary_MM.*\\.csv$", full.names = TRUE)[1]
  d5 <- read.csv(p5, check.names = FALSE)
  f5 <- as.character(d5$calendar_date)
  cat("formato de calendar_date (parte 5):", gsub("[0-9]", "9", f5[1]), "\n")
  fx5 <- as.Date(f5, tryFormats = c("%Y-%m-%d", "%d/%m/%Y"), optional = TRUE)
  d5$id <- id_de(d5$filename)
  d5$dia_desde_inicio <- as.integer(fx5 - do.call(c, primer_dia[d5$id]))
  patron5 <- "^(id|dia_desde_inicio|weekday|window_number|dur_spt_min|dur_day_min|nonwear_perc_day|nonwear_perc_spt|dur_day_total_(IN|LIG|MOD|VIG)_min)$"
  d5 <- d5[, grepl(patron5, names(d5)), drop = FALSE]
  write.csv(d5, file.path(outdir, "ggir_part5_dias.csv"), row.names = FALSE)
}
cat("GGIR sueño: listo\n")
}
main()
