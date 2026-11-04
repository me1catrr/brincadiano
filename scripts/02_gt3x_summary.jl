# 02 · Resumen por registro del lector GT3X (cobertura, huecos, validación).
#
# Qué hace
#   Recorre data/interim/ y resume, en una fila por registro, lo que escribió el
#   script 01: duración, porcentaje de muestras presentes, huecos de idle sleep frente
#   al resto, errores de checksum y número de registros de luz. Repite la comprobación
#   de contabilidad (registradas + ausentes = esperadas) y lo dice en consola.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/02_gt3x_summary.jl
#
# Depende de
#   01. Exige meta.arrow en cada carpeta de data/interim/.
#
# Lo usan
#   Ningún script posterior. La tabla entra en el capítulo 5 del informe después de
#   `make sync-assets` (report/), que la copia a
#   report/assets/04_formato_gt3x/brincadiano/.
#
# Entrada
#   data/interim/<ID>/{meta,acc,gaps,lux,diagnostics}.arrow
#
# Salida (una fila por ID; no se versiona)
#   outputs/gt3x/tables/resumen_registros.csv
#   outputs/gt3x/tables/resumen_registros.tex
#
# Consola: número de registros y si todas las contabilidades cuadran. Solo el ID.

using Brincadiano.GT3X
using Brincadiano.GT3X: expected_samples
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT = joinpath(ROOT, "outputs", "gt3x", "tables")

ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "meta.arrow")), readdir(INTERIM)))
isempty(ids) && error("No hay registros en data/interim/; ejecuta antes scripts/01_read_gt3x.jl")
rows = map(ids) do id
    r = read_interim(joinpath(INTERIM, id))
    exp = expected_samples(r)
    (id = id, formato = String(r.format), horas = exp / r.info.sample_rate / 3600,
     fs = r.info.sample_rate, escala = r.scale,
     registrado_pct = 100 * nrow(r.acc) / exp,
     huecos_idle = count(==("idle_sleep"), r.gaps.kind),
     huecos_otros = count(!=("idle_sleep"), r.gaps.kind),
     checksum_err = get(r.diagnostics, "checksum_errors", 0),
     registros_lux = nrow(r.lux),
     cuadra = nrow(r.acc) + sum(r.gaps.n_missing) == exp)
end
df = DataFrame(rows)
mkpath(OUT)
CSV.write(joinpath(OUT, "resumen_registros.csv"), df)

open(joinpath(OUT, "resumen_registros.tex"), "w") do io
    println(io, "% Generado por scripts/02_gt3x_summary.jl — no editar a mano; no versionar.")
    println(io, "\\begin{tabular}{@{}lrrrrrr@{}}\n\\toprule")
    println(io, "ID & Horas & Registrado & Huecos \\emph{idle sleep} & Otros huecos & Errores \\emph{checksum} & Registros de luz \\\\\n\\midrule")
    for r in eachrow(df)
        @printf(io, "\\texttt{%s} & %.1f & \\SI{%.1f}{\\percent} & %d & %d & %d & %d \\\\\n",
                replace(r.id, "_" => "\\_"), r.horas, r.registrado_pct, r.huecos_idle, r.huecos_otros,
                r.checksum_err, r.registros_lux)
    end
    println(io, "\\bottomrule\n\\end{tabular}")
end
println("Resumen de ", nrow(df), " registros en outputs/gt3x/tables/ (todas las contabilidades cuadran: ", all(df.cuadra), ")")
