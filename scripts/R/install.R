# Entorno R de los benchmarks. No forma parte del pipeline Julia y no lee .gt3x.
#
# Lo usan, como procesos lanzados desde Julia (directorio de trabajo scripts/R, para
# que renv active la librería):
#   90  read.gt3x, embebido en el propio script
#   92  ggir_benchmark.R          GGIR, partes 1–2
#   93 y 96  ggir_sueno.R         GGIR, partes 1–4 y, el 96, también la 5
#   94  lightlogr_metricas.R      LightLogR
#   95  read.gt3x con imputeZeroes = TRUE, embebido en el propio script
#
# Reproducir en otra máquina:
#   brew install r cmake pkgconf udunits gdal geos proj abseil
#   cd scripts/R && Rscript -e 'renv::restore()'
#
# Creación original (2026-10-02, R 4.6.1 de Homebrew). LightLogR depende de suntools y
# de sf, que necesitan la pila geoespacial (GDAL/GEOS/PROJ, udunits, Abseil/cmake para
# s2). Como renv revierte toda la instalación si un paquete falla, se instaló en dos
# pasos. Versiones fijadas en renv.lock: GGIR 3.3.9, read.gt3x 1.2.0, LightLogR 0.10.6.
# La librería (scripts/R/renv/library/) no se versiona; sí renv.lock.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = parallel::detectCores())
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!file.exists("renv.lock")) renv::init(bare = TRUE, restart = FALSE)
renv::install(c("GGIR", "read.gt3x"), prompt = FALSE)
renv::install("LightLogR", prompt = FALSE)
renv::snapshot(packages = c("GGIR", "read.gt3x", "LightLogR"), prompt = FALSE)
# Versiones fijadas: GGIR 3.3.9, read.gt3x 1.2.0, LightLogR 0.10.6 (ver renv.lock).
