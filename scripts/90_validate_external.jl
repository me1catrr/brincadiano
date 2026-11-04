# 90 · Validación del lector frente a read.gt3x (R) y pygt3x (Python).
#
# Qué hace
#   Vuelve a leer cada .gt3x con Brincadiano.GT3X y lo compara, muestra a muestra, con
#   read.gt3x (imputeZeroes = FALSE) y con pygt3x. Mide la diferencia absoluta máxima
#   sin redondear y tras redondear a 3 decimales con empates lejos de cero (la
#   convención de ActiLife) y cuenta las muestras que el lector externo añade, porque
#   imputa huecos, o que pierde. No forma parte del pipeline 01–10: no escribe
#   data/interim/ y no lee lo que haya escrito el 01.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/90_validate_external.jl --confirmo [--r] [--python] [ID ...]
#   Sin --confirmo solo describe lo que haría y termina. Hace falta autorización
#   explícita: lanza R y Python sobre data/raw/. Sin --r ni --python ejecuta los dos.
#
# Depende de
#   data/raw/*.gt3x.
#   R 4.6 + renv en scripts/R/ (read.gt3x 1.2.0; lo prepara scripts/R/install.R).
#   Python 3.12 en .venv/ (pygt3x 0.7.1; scripts/python/requirements.txt).
#   Los dos se lanzan como procesos (Rscript y .venv/bin/python), sin RCall ni
#   PythonCall. El paquete Julia no depende de ellos.
#
# Lo usan
#   Ningún script. El 95 repite la idea de la imputación (ceros de read.gt3x y última
#   muestra de pygt3x/ActiLife) sobre el pipeline completo, no sobre este CSV. La
#   coincidencia con ActiLife, que repite la última muestra, está en test/actilife_csv.jl
#   (maxdiff_filled = 0) y no pasa por aquí. Capítulo 14, vía `make sync-assets`.
#
# Entrada
#   data/raw/*.gt3x, o los ID nombrados en la línea de comandos.
#
# Salida (no se versiona)
#   outputs/validation/external_readers.csv
#       una fila por registro y lector: muestras, diferencia máxima y muestras añadidas
#   Los CSV intermedios de R y Python son temporales, fuera del repositorio, y se borran
#   al terminar. Contienen solo tiempo y aceleración.
#
# Convenciones de los lectores externos
#   read.gt3x con imputeZeroes = FALSE devuelve solo las muestras registradas; el tiempo
#   POSIXct es hora local tratada como UTC.
#   pygt3x rellena el idle sleep con la última muestra y lo marca en IdleSleepMode; el
#   Timestamp va en segundos Unix de la hora local.
# Consola: solo el ID y las diferencias. info.txt no se lee fuera de Brincadiano.GT3X.

using Brincadiano.GT3X
using Brincadiano.Rutas: Rutas, raw_files, find_raw
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const OUT = joinpath(ROOT, "outputs", "validation")
const RENV = joinpath(ROOT, "scripts", "R")                  # proyecto renv
const PYTHON = joinpath(ROOT, ".venv", "bin", "python")

const R_SCRIPT = raw"""
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(read.gt3x))
x <- read.gt3x(args[1], imputeZeroes = FALSE, asDataFrame = TRUE, verbose = FALSE)
out <- data.frame(t = format(as.numeric(x$time), digits = 15), x = x$X, y = x$Y, z = x$Z)
write.csv(out, args[2], row.names = FALSE, quote = FALSE)
"""

const PY_SCRIPT = raw"""
import sys
from pygt3x.reader import FileReader
with FileReader(sys.argv[1]) as reader:
    df = reader.to_pandas()
df = df.reset_index()
ism = df["IdleSleepMode"].astype(int) if "IdleSleepMode" in df.columns else 0
out = df.assign(ism=ism)[["Timestamp", "X", "Y", "Z", "ism"]]
out.columns = ["t", "x", "y", "z", "ism"]
out.to_csv(sys.argv[2], index=False, float_format="%.6f")
"""

round3(v) = round(Float64(v), RoundNearestTiesAway; digits = 3)

"""Índice de muestra (1-based) desde info.start, a partir de segundos Unix locales."""
sample_index(t_unix::AbstractVector{<:Real}, start::DateTime, fs) =
    round.(Int, (t_unix .- datetime2unix(start)) .* fs) .+ 1

function ours_index(r)
    fs = r.info.sample_rate
    ms = Dates.value.(r.acc.time .- r.info.start)
    return (ms .÷ 1000) .* fs .+ (mod.(ms, 1000) .* fs) .÷ 1000 .+ 1
end

function compare(r, ext::DataFrame; drop_flagged = false)
    fs = r.info.sample_rate
    drop_flagged && hasproperty(ext, :ism) && (ext = ext[ext.ism .== 0, :])
    ie = sample_index(ext.t, r.info.start, fs)
    io = ours_index(r)
    pos = Dict(i => k for (k, i) in enumerate(io))
    common = [(k, pos[i]) for (k, i) in enumerate(ie) if haskey(pos, i)]
    ke = first.(common); ko = last.(common)
    raw = maximum(c -> maximum(abs.(Float64.(ext[ke, c]) .- Float64.(r.acc[ko, c]))), (:x, :y, :z); init = 0.0)
    rnd = maximum(c -> maximum(abs.(round3.(ext[ke, c]) .- round3.(r.acc[ko, c]))), (:x, :y, :z); init = 0.0)
    return (n_ours = nrow(r.acc), n_ext = length(ie), n_common = length(common),
            n_ext_only = length(ie) - length(common), n_ours_only = nrow(r.acc) - length(common),
            maxdiff_raw = raw, maxdiff_round3 = rnd)
end

function run_external(tool, gt3x, csvout)
    script = tempname() * (tool == "read.gt3x" ? ".R" : ".py")
    write(script, tool == "read.gt3x" ? R_SCRIPT : PY_SCRIPT)
    # Rscript se ejecuta desde scripts/R para que renv active su librería fijada.
    cmd = tool == "read.gt3x" ? Cmd(`Rscript $script $gt3x $csvout`; dir = RENV) :
                                `$PYTHON $script $gt3x $csvout`
    try
        run(pipeline(cmd; stdout = devnull))
    finally
        rm(script; force = true)
    end
    return CSV.read(csvout, DataFrame)
end

function main(args)
    if !("--confirmo" in args)
        println("""
        Validación externa NO ejecutada (falta --confirmo).
        Este script lanzaría read.gt3x (R) y/o pygt3x (Python) sobre data/raw/*.gt3x como
        procesos locales y compararía muestra a muestra con Brincadiano.GT3X.
        Ejecutar solo con autorización explícita del investigador.""")
        return
    end
    tools = String[]
    "--r" in args && push!(tools, "read.gt3x")
    "--python" in args && push!(tools, "pygt3x")
    isempty(tools) && (tools = ["read.gt3x", "pygt3x"])
    "read.gt3x" in tools && Sys.which("Rscript") === nothing && error("Rscript no está disponible")
    "pygt3x" in tools && !isfile(PYTHON) && error("Falta .venv/ (ver scripts/python/requirements.txt)")
    ids = filter(a -> !startswith(a, "--"), args)
    files = isempty(ids) ? raw_files("gt3x") : find_raw.(ids .* ".gt3x")
    rows = DataFrame()
    for f in files
        id = splitext(basename(f))[1]
        r = read_gt3x(f)
        for tool in tools
            tmp = tempname() * ".csv"
            try
                ext = run_external(tool, f, tmp)
                c = compare(r, ext; drop_flagged = tool == "pygt3x")
                push!(rows, (; id, tool, c...); cols = :union)
                @printf("%-12s %-10s comunes=%d solo_ext=%d solo_nuestras=%d maxdiff=%.3g (redondeo 3 dec: %.3g)\n",
                        id, tool, c.n_common, c.n_ext_only, c.n_ours_only, c.maxdiff_raw, c.maxdiff_round3)
            finally
                rm(tmp; force = true)
            end
        end
    end
    mkpath(OUT)
    CSV.write(joinpath(OUT, "external_readers.csv"), rows)
    println("Resumen en outputs/validation/external_readers.csv")
end

main(ARGS)
