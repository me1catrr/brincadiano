# 95 · Efecto de imputar el idle sleep sobre el pipeline.
#
# Qué hace
#   Repite, en memoria, el tramo 04–06 del pipeline con tres variantes de la misma
#   señal y escribe la diferencia. No reescribe data/interim/ ni los epochs.
#   «huecos» es la serie de Brincadiano (el idle sleep no se rellena). «ceros» es
#   GT3X.impute_idle_sleep con method = :zeros, como read.gt3x(imputeZeroes = TRUE).
#   «última» repite la última muestra, como ActiLife y pygt3x. Para cada variante:
#   autocalibración, no uso, ENMO, días válidos (≥ 600 min), L5/M10, IS/IV, episodios
#   vigorosos (464,6 mg) y θ, R. Desde H28, θ y R usan el peso corregido por la
#   cobertura del reloj (coverage_weights), el mismo del script 06. En el registro corto (el que empieza por VTL) comprueba
#   además que «ceros» coincide muestra a muestra con read.gt3x. La variante «última»
#   frente a ActiLife ya está en test/actilife_csv.jl (maxdiff_filled = 0) y no se
#   relanza aquí.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/95_imputacion.jl --confirmo
#   Sin --confirmo termina tras un aviso. Hace falta autorización: el registro corto
#   lanza read.gt3x en R.
#
# Depende de
#   01  data/interim/<ID>/ y, para la comprobación, data/raw/<ID>.gt3x
#   05  solo como filtro: entra en la tabla quien tenga epochs_60s.arrow y no sea el
#       registro corto. Los números se recalculan; no se leen las tablas del 04 ni del 06.
#   scripts/R/ (read.gt3x 1.2.0), igual que el 90, que es la comparación del lector
#       sin imputar. Este script es la comparación del pipeline con imputación.
#
# Lo usan
#   El informe, tras `make sync-assets`. Ningún script posterior lee estas tablas.
#
# Entrada
#   data/interim/<ID>/{acc,gaps,lux,capsense,meta}.arrow
#   data/raw/<VTL>.gt3x    solo la comprobación contra read.gt3x
#
# Salida (individual; no se versiona)
#   outputs/validation/tables/imputacion.csv       una fila por registro y variante
#   outputs/validation/tables/imputacion.tex       las tres variantes en columnas
#   outputs/validation/tables/imputacion_readgt3x.csv
#       muestras, comunes, diferencia máxima y recuento de ceros frente a read.gt3x
#
# Consola: solo el ID, la variante y los estadísticos. El CSV temporal de R se borra.

using Brincadiano.GT3X, Brincadiano.QC, Brincadiano.Metrics, Brincadiano.Circadian
using Brincadiano.Rutas: Rutas, find_raw
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT = joinpath(ROOT, "outputs", "validation", "tables")
const RENV = joinpath(ROOT, "scripts", "R")
const MIN_USO = 600
const VARIANTES = (("huecos", nothing), ("ceros", :zeros), ("última", :last))

fmt(x; d = 1) = (x isa Real && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
hhmm(h) = isnan(h) ? "---" : Dates.format(Time(0) + Minute(round(Int, 60h) % 1440), "HH:MM")
dif_h(a, b) = (d = mod(a - b, 24); d > 12 ? d - 24 : d)

function ejecutar(r, acc, gaps)
    fs = r.info.sample_rate
    ac = autocalibrate(acc, fs)
    off, sc = ac.aplicar ? (ac.offset, ac.scale) : ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    q = qc_report(acc, gaps, r.lux, r.capsense, fs; clip = something(r.info.accel_max, 8.0) - 0.5)
    runs = q.tramos
    o = floor(first(r.acc.time), Day); stop = last(r.acc.time)
    e60 = epoch_metrics(acc, gaps, fs, Second(60); origin = o, stop = stop, offset = off, scale = sc); flag_nonwear!(e60, runs)
    e1 = epoch_metrics(acc, gaps, fs, Second(1); origin = o, stop = stop, offset = off, scale = sc); flag_nonwear!(e1, runs)
    dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    pm = profile24h(e60, validos; nonwear_as = :missing); pr = profile24h(e60, validos; nonwear_as = :rest)
    rm_, rr = rar_metrics(pm), rar_metrics(pr)
    ancla = rm_.estimable ? rm_ : rr
    iv = is_iv(e60, validos; nonwear_as = rm_.estimable ? :missing : :rest)
    ep = vigorous_episodes(e1, validos; threshold = 464.6)
    # Misma corrección que el θ principal (script 06): la dosis de cada episodio
    # se divide por los días en que se observó su hora del reloj.
    cob = clock_coverage(e1, validos; nonwear_as = rm_.estimable ? :missing : :rest)
    ep.peso = coverage_weights(ep, cob)
    tm = timing(ep, ancla.L5_medio_h; peso = :peso)
    enmo = sum(dd.enmo_medio_mg[dd.minutos_uso .>= MIN_USO]; init = 0.0) / max(length(validos), 1)
    return (cal_aplicada = ac.aplicar, cal_error_despues_mg = ac.error_despues_mg,
            nouso_h = sum(runs.horas[runs.etiqueta .== "no uso probable"]; init = 0.0),
            dias_validos = length(validos), enmo_medio_mg = enmo,
            L5_estimable = rm_.estimable, ancla_h = ancla.L5_medio_h, RA = ancla.RA, IS = iv.IS, IV = iv.IV,
            episodios = tm.n, θ_h = tm.θ_h, R = tm.R)
end

# Comprobación de la variante «ceros» frente a read.gt3x(imputeZeroes = TRUE).
const R_SCRIPT = raw"""
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(read.gt3x))
x <- read.gt3x(args[1], imputeZeroes = TRUE, asDataFrame = TRUE, verbose = FALSE)
out <- data.frame(t = format(as.numeric(x$time), digits = 15), x = x$X, y = x$Y, z = x$Z)
write.csv(out, args[2], row.names = FALSE, quote = FALSE)
"""

function comprobar_readgt3x(id)
    r = read_gt3x(find_raw(id * ".gt3x"))
    fs = r.info.sample_rate
    a0, _ = impute_idle_sleep(r.acc, r.gaps, fs; method = :zeros)
    script = tempname() * ".R"; csv = tempname() * ".csv"
    write(script, R_SCRIPT)
    try
        run(pipeline(Cmd(`Rscript $script $(find_raw(id * ".gt3x")) $csv`; dir = RENV); stdout = devnull, stderr = devnull))
        ext = CSV.read(csv, DataFrame)
        # Índice de muestra desde el inicio, como en scripts/90 (tiempos locales en segundos Unix).
        ie = round.(Int, (ext.t .- datetime2unix(r.info.start)) .* fs) .+ 1
        ms = Dates.value.(a0.time .- r.info.start)
        io = (ms .÷ 1000) .* fs .+ (mod.(ms, 1000) .* fs) .÷ 1000 .+ 1
        pos = Dict(i => k for (k, i) in enumerate(io))
        com = [(k, pos[i]) for (k, i) in enumerate(ie) if haskey(pos, i)]
        ke, ko = first.(com), last.(com)
        dif = maximum(c -> maximum(abs.(Float64.(ext[ke, c]) .- Float64.(a0[ko, c]))), (:x, :y, :z); init = 0.0)
        ceros_ext = count((ext.x .== 0) .& (ext.y .== 0) .& (ext.z .== 0))
        ceros_nuestros = count((a0.x .== 0) .& (a0.y .== 0) .& (a0.z .== 0))
        return (id = id, n_readgt3x = nrow(ext), n_nuestra = nrow(a0), n_comunes = length(com),
                max_dif_g = dif, ceros_readgt3x = ceros_ext, ceros_nuestros = ceros_nuestros)
    finally
        rm(script; force = true); rm(csv; force = true)
    end
end

function tabla_tex(df, ids)
    io = IOBuffer()
    vs = first.(VARIANTES)
    println(io, "% Generado por scripts/95_imputacion.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X *{$(3length(ids))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(["\\multicolumn{3}{c}{$(replace(id, "_" => "\\_"))}" for id in ids], " & "), " \\\\")
    println(io, join(["\\cmidrule(lr){$(2 + 3(k-1))-$(4 + 3(k-1))}" for k in eachindex(ids)], ""))
    println(io, " & ", join(repeat(["huecos", "ceros", "última"], length(ids)), " & "), " \\\\")
    println(io, "\\midrule")
    g(id, v) = df[(df.id .== id) .& (df.variante .== v), :][1, :]
    fila(nombre, f) = println(io, nombre, " & ", join([f(g(id, v)) for id in ids for v in vs], " & "), " \\\\")
    fila("Error de calibración tras ajustar (m\\BRg)", x -> fmt(x.cal_error_despues_mg) * (x.cal_aplicada ? "" : "\\textsuperscript{*}"))
    fila("No uso probable (h)", x -> fmt(x.nouso_h))
    fila("Días válidos", x -> string(x.dias_validos))
    fila("ENMO medio (m\\BRg)", x -> fmt(x.enmo_medio_mg))
    fila("Ancla (punto medio de L5)", x -> hhmm(x.ancla_h) * (x.L5_estimable ? "" : "\\textsuperscript{†}"))
    fila("IS / IV", x -> fmt(x.IS; d = 2) * " / " * fmt(x.IV; d = 2))
    fila("Episodios vigorosos", x -> string(x.episodios))
    fila("\$\\bar\\theta\$ (h tras el ancla) / \$R\$", x -> fmt(x.θ_h) * " / " * fmt(x.R; d = 2))
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

if !("--confirmo" in ARGS)
    println("Contraste R1 (imputación del idle sleep). Lanza read.gt3x en R: añade --confirmo.")
    exit(0)
end
mkpath(OUT)
ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "epochs_60s.arrow")) && find_raw(d * ".gt3x") !== nothing &&
                       !startswith(d, "VTL"), readdir(INTERIM)))
filas = NamedTuple[]
for id in ids
    r = read_interim(joinpath(INTERIM, id)); fs = r.info.sample_rate
    for (nombre, m) in VARIANTES
        acc, gaps = m === nothing ? (r.acc, r.gaps) : impute_idle_sleep(r.acc, r.gaps, fs; method = m)
        x = ejecutar(r, acc, gaps)
        push!(filas, (id = id, variante = nombre, x...))
        @printf("%-12s %-7s cal %.1f mg%s · no uso %.1f h · %d días · ENMO %.1f mg · ancla %s%s · IS %.2f IV %.2f · %d episodios · θ %.1f h R %.2f\n",
                id, nombre, x.cal_error_despues_mg, x.cal_aplicada ? "" : " (no aplicada)", x.nouso_h, x.dias_validos, x.enmo_medio_mg,
                hhmm(x.ancla_h), x.L5_estimable ? "" : " (L5*)", x.IS, x.IV, x.episodios, x.θ_h, x.R)
    end
end
df = DataFrame(filas)
CSV.write(joinpath(OUT, "imputacion.csv"), df)
write(joinpath(OUT, "imputacion.tex"), tabla_tex(df, ids))
vtl = filter(d -> startswith(d, "VTL") && find_raw(d * ".gt3x") !== nothing, readdir(INTERIM))
if !isempty(vtl)
    try
        c = comprobar_readgt3x(first(vtl))
        CSV.write(joinpath(OUT, "imputacion_readgt3x.csv"), DataFrame([c]))
        @printf("%-12s read.gt3x(imputeZeroes): %d muestras frente a %d nuestras; %d comunes; dif. máx. %.2g g; ceros %d frente a %d\n",
                c.id, c.n_readgt3x, c.n_nuestra, c.n_comunes, c.max_dif_g, c.ceros_readgt3x, c.ceros_nuestros)
    catch e
        @printf("%-12s read.gt3x(imputeZeroes) no se pudo contrastar (%s). La tabla de imputación ya está escrita.\n",
                first(vtl), sprint(showerror, e)[1:min(end, 180)])
    end
end
println("Tablas en ", relpath(OUT, ROOT))
