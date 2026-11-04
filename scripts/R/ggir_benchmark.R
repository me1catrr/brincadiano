# Benchmark de GGIR para calibración, ENMO, no uso y ritmo (partes 1 y 2).
#
# Lo llama scripts/92_validar_ggir.jl. No está pensado para ejecutarse a mano: GGIR
# lee el info.txt completo del .gt3x. Sin la autorización que exige el script 92, no lanzarlo.
#
#   Rscript ggir_benchmark.R <directorio_salida> <archivo1.gt3x> [archivo2.gt3x ...]
#   (el directorio de trabajo tiene que ser scripts/R, para que renv active la librería)
#
# Qué hace
#   GGIR 3.3.9, mode 1:2, epochs de 5 s, ventanas de 15 min, ENMO, MAD, calibración,
#   no uso «2023», L5/M10 e IS/IV (cosinor = TRUE). Trabaja en un directorio temporal
#   que se borra al salir, haya error o no. De la salida solo se copian columnas de
#   una lista blanca; el ID sale del nombre del archivo (idloc = 2).
#
# Depende de
#   El renv de este directorio (scripts/R/install.R, GGIR 3.3.9). No lee salidas de Julia.
#
# Lo usa
#   Solo el 92. El sueño y la parte 5 van por ggir_sueno.R (scripts 93 y 96).
#
# Entrada
#   Uno o más .gt3x. El primer argumento es el directorio de salida
#   (el 92 le pasa outputs/validation/ggir/).
#
# Salida (la escribe el 92 bajo outputs/validation/ggir/; no se versiona)
#   ggir_calibracion.csv                  error y escala, sin filename
#   ggir_<ID>_epochs5s.csv                timestamp, ENMO, MAD
#   ggir_<ID>_ventanas15min.csv           timestamp, nonwearscore, clippingscore
#   ggir_part2_summary.csv                IS, IV y medias de ENMO de la lista blanca
#   ggir_part2_daysummary.csv             L5/M10 por día, de la misma lista
#
# Consola: «GGIR: listo» y, una vez, los nombres de columna de ritmo de part2_summary.
# No imprime valores de la salida de GGIR.

main <- function() {
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
files <- normalizePath(args[-1])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tmp <- file.path(tempdir(), paste0("ggir_", Sys.getpid()))
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

suppressMessages(library(GGIR))
invisible(capture.output(suppressMessages(GGIR(
  mode = 1:2, datadir = files, outputdir = tmp, studyname = "brincadiano",
  overwrite = TRUE, do.report = 2, idloc = 2, print.filename = FALSE,
  storefolderstructure = FALSE, do.parallel = FALSE, verbose = FALSE,
  desiredtz = "Europe/Madrid", windowsizes = c(5, 900, 3600),
  do.enmo = TRUE, do.mad = TRUE, do.cal = TRUE, nonwear_approach = "2023",
  winhr = c(5, 10), includedaycrit = 10, epochvalues2csv = FALSE,
  visualreport = FALSE, cosinor = TRUE      # cosinor = TRUE: activa IS/IV (g.IVIS) en la parte 2
))))

base <- list.files(tmp, pattern = "^output_", full.names = TRUE)[1]
id_de <- function(f) sub("(\\.gt3x)?(\\.RData)?$", "", sub("^meta_", "", basename(f)))

# 1. Calibración (informe de calidad de GGIR), columnas de la lista blanca.
qc <- read.csv(file.path(base, "results", "QC", "data_quality_report.csv"), check.names = FALSE)
keep <- intersect(c("filename", "cal.error.start", "cal.error.end", "scale.x", "scale.y", "scale.z",
                    "offset.x", "offset.y", "offset.z", "n.hours.considered", "QCmessage"), names(qc))
qc <- qc[, keep, drop = FALSE]
qc$id <- id_de(qc$filename); qc$filename <- NULL
write.csv(qc, file.path(outdir, "ggir_calibracion.csv"), row.names = FALSE)

# 2. Series de la parte 1 (sin imputar): epochs de 5 s y ventanas de 15 min.
for (f in list.files(file.path(base, "meta", "basic"), full.names = TRUE)) {
  e <- new.env(); load(f, envir = e); M <- e$M
  id <- id_de(f)
  ms <- M$metashort[, intersect(c("timestamp", "ENMO", "MAD"), names(M$metashort)), drop = FALSE]
  ml <- M$metalong[, intersect(c("timestamp", "nonwearscore", "clippingscore"), names(M$metalong)), drop = FALSE]
  write.csv(ms, file.path(outdir, paste0("ggir_", id, "_epochs5s.csv")), row.names = FALSE)
  write.csv(ml, file.path(outdir, paste0("ggir_", id, "_ventanas15min.csv")), row.names = FALSE)
}

# 3. Resúmenes de la parte 2: solo columnas de ritmo y volumen (lista blanca por patrón).
patron <- "^(ID|filename|calendar_date|weekday)$|(^|_)IS(_|$)|(^|_)IV(_|$)|interdaily|intradaily|^L5|^M5|^L10|^M10|AD_mean_ENMO|mean_ENMO_mg_24hr|ENMO_fullRecordingMean|wear_dur_def_proto_day"
for (nm in c("part2_summary.csv", "part2_daysummary.csv")) {
  p <- file.path(base, "results", nm)
  if (!file.exists(p)) next
  d <- read.csv(p, check.names = FALSE)
  # Solo nombres de columnas candidatas (nunca valores), para documentar qué se extrae.
  if (nm == "part2_summary.csv") cat("columnas de ritmo en part2_summary:",
      paste(grep("IS|IV|interd|intrad", names(d), value = TRUE), collapse = ", "), "\n")
  d <- d[, grepl(patron, names(d)), drop = FALSE]
  if ("filename" %in% names(d)) { d$id <- id_de(d$filename); d$filename <- NULL }
  write.csv(d, file.path(outdir, paste0("ggir_", nm)), row.names = FALSE)
}
cat("GGIR: listo\n")
}
main()
