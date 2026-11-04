# 92 · Validación del pipeline frente a GGIR y de la media circular frente a numpy.
#
# Qué hace
#   1. Lanza GGIR 3.3.9 (scripts/R/ggir_benchmark.R, partes 1 y 2) sobre data/raw/*.gt3x.
#      GGIR trabaja en un directorio temporal; aquí solo llega una lista blanca de
#      columnas, a outputs/validation/ggir/.
#   2. Compara con lo ya calculado por Brincadiano: autocalibración (recalculada sobre
#      el interim), ENMO y MAD de 5 s en epochs con muestreo completo y en uso, no uso
#      en ventanas de 15 min, y L5/M10 e IS/IV del resumen del script 06.
#   3. Recalcula θ y R con numpy a partir de los episodios del 06, como implementación
#      independiente de la media circular ponderada.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/92_validar_ggir.jl --confirmo
#   julia --project=scripts scripts/92_validar_ggir.jl --confirmo --sin-ggir
#   Sin --confirmo solo describe lo que haría. Hace falta autorización: lanza R (y
#   Python para numpy) sobre los .gt3x. --sin-ggir reutiliza outputs/validation/ggir/
#   si ya existe ggir_calibracion.csv; sirve para repetir la comparación sin relanzar GGIR.
#
# Depende de
#   01  data/interim/<ID>/ (aceleración, para calibrar y para las ventanas de no uso)
#   05  epochs_5s.arrow
#   06  outputs/circadian/tables/ritmo_resumen.csv
#       outputs/timing/tables/episodios_<ID>.csv y timing_resumen.csv
#   scripts/R/ggir_benchmark.R y el renv de scripts/R/ (GGIR 3.3.9)
#   .venv/bin/python con numpy, para el paso 3
#   No usa las tablas del 04: el no uso de 15 min se vuelve a calcular aquí.
#
# Lo usan
#   10, si existen los cruces por epoch y por ventana de 15 min, para Bland-Altman y
#   la línea de tiempo de no uso (capítulo 14). `make sync-assets` copia la tabla .tex.
#
# Entrada
#   data/raw/*.gt3x, más los Arrow y CSV de «Depende de».
#
# Salida de la lista blanca de GGIR (no se versiona; el temporal de GGIR se borra)
#   outputs/validation/ggir/ggir_calibracion.csv
#   outputs/validation/ggir/ggir_<ID>_epochs5s.csv
#   outputs/validation/ggir/ggir_<ID>_ventanas15min.csv
#   outputs/validation/ggir/ggir_part2_summary.csv
#   outputs/validation/ggir/ggir_part2_daysummary.csv
# Cruces que lee el script 10 (no se versionan)
#   outputs/validation/ggir_<ID>_enmo5s_cruce.csv   ENMO y MAD (Bland–Altman, script 10)
#   outputs/validation/ggir_<ID>_nouso15min_cruce.csv
# Resumen de la comparación (no se versiona)
#   outputs/validation/tables/ggir_calibracion.csv
#   outputs/validation/tables/ggir_enmo_mad.csv
#   outputs/validation/tables/ggir_no_uso.csv
#   outputs/validation/tables/ggir_ritmo.csv
#   outputs/validation/tables/ggir_circular.csv
#   outputs/validation/tables/ggir_comparacion.tex
#
# Consola: solo el ID y los estadísticos (correlación, acuerdo, IS/IV, diferencia de θ).
# GGIR lee el info.txt completo en su temporal; ese directorio se borra al terminar y
# de su salida solo se copian las columnas de la lista blanca.

using Brincadiano.GT3X, Brincadiano.QC, Brincadiano.Metrics, Brincadiano.Circadian
using Brincadiano.Rutas: Rutas, raw_files

using CSV, DataFrames, Dates, Printf

# Notación científica para LaTeX con coma decimal: 9,8·10⁻¹².
sci(x) = (e = floor(Int, log10(abs(x))); m = abs(x) / 10.0^e;
          "\$" * replace(@sprintf("%.1f", m), "." => "{,}", "-" => "−") * "\\cdot10^{$e}\$")

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const INTERIM = joinpath(ROOT, "data", "interim")
const GG = joinpath(ROOT, "outputs", "validation", "ggir")
const VALID = joinpath(ROOT, "outputs", "validation")
const OUT = joinpath(ROOT, "outputs", "validation", "tables")
const RENV = joinpath(ROOT, "scripts", "R")
const PYTHON = joinpath(ROOT, ".venv", "bin", "python")

# Coma decimal (convención del informe).
fmt(x; d = 2) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
mean_(v) = isempty(v) ? NaN : sum(v) / length(v)
function pearson(a, b)
    ma, mb = mean_(a), mean_(b)
    return sum((a .- ma) .* (b .- mb)) / sqrt(sum((a .- ma) .^ 2) * sum((b .- mb) .^ 2))
end
pctl(v, p) = (s = sort(v); s[clamp(ceil(Int, p * length(s)), 1, length(s))])
ggtime(s) = DateTime(first(s, 19))                     # "2018-09-25T17:00:00+0200" → hora local

function correr_ggir()
    files = raw_files("gt3x")
    err = tempname()
    try
        run(pipeline(Cmd(`Rscript ggir_benchmark.R $GG $files`; dir = RENV); stdout = devnull, stderr = err))
    catch
        # Solo las líneas de error (rutas y mensajes de GGIR; sin datos del encabezado).
        foreach(l -> occursin(r"Error"i, l) && println("  ", l), eachline(err))
        rethrow()
    finally
        rm(err; force = true)
    end
end

function comparar(id, cal_gg, p2, d2)
    dir = joinpath(INTERIM, id)
    r = read_interim(dir); fs = r.info.sample_rate
    # 1. Calibración.
    ac = autocalibrate(r.acc, fs)
    g = cal_gg[cal_gg.id .== id, :]
    cal = (id = id, err_antes_nuestro = ac.error_antes_mg, err_despues_nuestro = ac.error_despues_mg,
           err_antes_ggir = 1000g.var"cal.error.start"[1], err_despues_ggir = 1000g.var"cal.error.end"[1],
           max_dif_offset_mg = maximum(abs.(1000 .* collect(ac.offset) .- 1000 .* [g.var"offset.x"[1], g.var"offset.y"[1], g.var"offset.z"[1]])),
           max_dif_escala = maximum(abs.(collect(ac.scale) .- [g.var"scale.x"[1], g.var"scale.y"[1], g.var"scale.z"[1]])))
    # 2. ENMO y MAD en 5 s.
    e5 = load_epochs(dir, 5)
    gg = CSV.read(joinpath(GG, "ggir_$(id)_epochs5s.csv"), DataFrame)
    gg.start = ggtime.(String.(gg.timestamp))
    j = innerjoin(e5, select(gg, :start, :ENMO => :enmo_gg, :MAD => :mad_gg), on = :start)
    ok = (j.cov .== 1) .& .!j.nonwear
    a, b = j.enmo[ok], 1000 .* j.enmo_gg[ok]
    ma, mb = j.mad[ok], 1000 .* j.mad_gg[ok]
    # Cruces por epoch y por ventana, para las figuras de validación (scripts/10_figuras.jl).
    CSV.write(joinpath(VALID, "ggir_$(id)_enmo5s_cruce.csv"),
              DataFrame(start = j.start[ok], enmo = a, enmo_ggir = b, mad = ma, mad_ggir = mb))
    idle = (j.n .== 0) .& (j.idle_frac .== 1) .& .!j.nonwear
    enmo = (id = id, epochs = count(ok), r_enmo = pearson(a, b), dif_media_mg = mean_(a .- b),
            mae_mg = mean_(abs.(a .- b)), p99_mg = pctl(abs.(a .- b), 0.99),
            pct_1mg = 100 * count(abs.(a .- b) .< 1) / length(a), r_mad = pearson(ma, mb),
            mae_mad_mg = mean_(abs.(ma .- mb)), epochs_idle = count(idle),
            ggir_idle_enmo_medio_mg = 1000mean_(j.enmo_gg[idle]))
    # 3. No uso en ventanas de 15 min.
    w = nonwear_windows(r.acc, r.gaps, fs)
    gl = CSV.read(joinpath(GG, "ggir_$(id)_ventanas15min.csv"), DataFrame)
    gl.start = ggtime.(String.(gl.timestamp))
    k = innerjoin(select(w, :start, :score, :nonwear), select(gl, :start, :nonwearscore), on = :start)
    k = k[.!ismissing.(k.score), :]
    gnw = k.nonwearscore .>= 2
    CSV.write(joinpath(VALID, "ggir_$(id)_nouso15min_cruce.csv"), DataFrame(start = k.start, nuestro = k.nonwear, ggir = gnw))
    po = count(gnw .== k.nonwear) / nrow(k)
    pe = mean_(gnw) * mean_(k.nonwear) + (1 - mean_(gnw)) * (1 - mean_(k.nonwear))
    nw = (id = id, ventanas = nrow(k), acuerdo_pct = 100po, kappa = pe < 1 ? (po - pe) / (1 - pe) : NaN,
          nouso_nuestro_h = count(k.nonwear) / 4, nouso_ggir_h = count(gnw) / 4,
          score_identico_pct = 100 * count(k.score .== k.nonwearscore) / nrow(k))
    # 4. Ritmo: L5/M10 (media de los días de GGIR) e IS/IV.
    rr = CSV.read(joinpath(ROOT, "outputs", "circadian", "tables", "ritmo_resumen.csv"), DataFrame)
    ro = rr[rr.id .== id, :]
    dd = d2[(d2.id .== id) .& .!ismissing.(d2.var"L5hr_ENMO_mg_0-24hr"), :]
    circ_h(hs) = (c = circular_mean(2π .* mod.(hs, 24) ./ 24); phase_to_hours(c.θ))
    s2 = p2[p2.id .== id, :]
    ritmo = (id = id, L5_medio_nuestro_h = ro.L5_estimable[1] ? ro.L5_medio_h[1] : ro.L5_reposo_medio_h[1],
             L5_medio_ggir_h = circ_h(Float64.(dd.var"L5hr_ENMO_mg_0-24hr") .+ 2.5),
             M10_medio_nuestro_h = ro.M10_medio_h[1], M10_medio_ggir_h = circ_h(Float64.(dd.var"M10hr_ENMO_mg_0-24hr") .+ 5),
             IS_nuestro = ro.IS[1], IS_reposo_nuestro = ro.IS_reposo[1], IS_ggir = s2.IS[1],
             IV_nuestro = ro.IV[1], IV_reposo_nuestro = ro.IV_reposo[1], IV_ggir = s2.IV[1])
    return cal, enmo, nw, ritmo
end

function circular_numpy(ids)
    # numpy si el Python 3.12 del .venv sigue instalado. Si el intérprete ya no está,
    # la misma fórmula con math: no se instala nada. El peso es el corregido por cobertura.
    py_numpy = raw"""
import sys, numpy as np, pandas as pd
for f in sys.argv[1:]:
    d = pd.read_csv(f)
    phi = 2 * np.pi * d["fase_h"].to_numpy() / 24
    w = (d["peso"] if "peso" in d else d["dosis"]).to_numpy()
    C = np.sum(w * np.cos(phi)) / np.sum(w); S = np.sum(w * np.sin(phi)) / np.sum(w)
    print(f"{np.mod(np.arctan2(S, C), 2*np.pi) * 24 / (2*np.pi):.10f} {np.hypot(C, S):.10f}")
"""
    py_stdlib = raw"""
import csv, math, sys
for f in sys.argv[1:]:
    with open(f, newline="") as fh:
        rows = list(csv.DictReader(fh))
    clave = "peso" if rows and "peso" in rows[0] else "dosis"
    phi = [2 * math.pi * float(r["fase_h"]) / 24 for r in rows]
    w = [float(r[clave]) for r in rows]
    sw = sum(w)
    C = sum(wi * math.cos(p) for wi, p in zip(w, phi)) / sw
    S = sum(wi * math.sin(p) for wi, p in zip(w, phi)) / sw
    th = math.atan2(S, C) % (2 * math.pi) * 24 / (2 * math.pi)
    print(f"{th:.10f} {math.hypot(C, S):.10f}")
"""
    files = [joinpath(ROOT, "outputs", "timing", "tables", "episodios_$(id).csv") for id in ids]
    if isfile(PYTHON)
        py, exe, motor = py_numpy, PYTHON, "numpy"
    else
        py, exe, motor = py_stdlib, something(Sys.which("python3"), "python3"), "math"
        println("numpy: falta ", PYTHON, "; contraste circular con math (misma fórmula, sin instalar)")
    end
    script = tempname() * ".py"; write(script, py)
    out = try readlines(`$exe $script $files`) finally rm(script; force = true) end
    tm = CSV.read(joinpath(ROOT, "outputs", "timing", "tables", "timing_resumen.csv"), DataFrame)
    return DataFrame(map(zip(ids, out)) do (id, l)
        θn, Rn = parse.(Float64, split(l))
        t = tm[tm.id .== id, :]
        (id = id, θ_julia_h = t.θ_h[1], θ_numpy_h = θn, R_julia = t.R[1], R_numpy = Rn,
         dif_θ_h = abs(mod(t.θ_h[1] - θn + 12, 24) - 12), dif_R = abs(t.R[1] - Rn), motor = motor)
    end)
end

function tablas_tex(cal, enmo, nw, ritmo, circ)
    open(joinpath(OUT, "ggir_comparacion.tex"), "w") do io
        println(io, "% Generado por scripts/92_validar_ggir.jl — tabla individual, no versionar.")
        println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X *{$(nrow(cal))}{r}@{}}")
        println(io, "\\toprule")
        println(io, " & ", join(replace.(cal.id, "_" => "\\_"), " & "), " \\\\")
        println(io, "\\midrule")
        fila(lab, v) = println(io, lab, " & ", join(v, " & "), " \\\\")
        sub(t) = println(io, "\\multicolumn{$(nrow(cal) + 1)}{@{}l}{\\textit{$t}} \\\\")
        sub("Calibración")
        fila("Error antes / después: \\Brincadiano{} (m\\BRg)", [fmt(a; d = 1) * " / " * fmt(b; d = 1) for (a, b) in zip(cal.err_antes_nuestro, cal.err_despues_nuestro)])
        fila("Error antes / después: GGIR (m\\BRg)", [fmt(a; d = 1) * " / " * fmt(b; d = 1) for (a, b) in zip(cal.err_antes_ggir, cal.err_despues_ggir)])
        fila("Diferencia máxima de desplazamiento (m\\BRg)", fmt.(cal.max_dif_offset_mg; d = 1))
        fila("Diferencia máxima de escala", fmt.(cal.max_dif_escala; d = 4))
        sub("ENMO y MAD en \\emph{epochs} de 5\\,s (muestreo completo, en uso)")
        fila("\\emph{Epochs} comparados", string.(enmo.epochs))
        fila("Correlación de ENMO", fmt.(enmo.r_enmo; d = 4))
        fila("Diferencia media / absoluta media de ENMO (m\\BRg)", [fmt(a) * " / " * fmt(b) for (a, b) in zip(enmo.dif_media_mg, enmo.mae_mg)])
        fila("\\emph{Epochs} con diferencia $(raw"$<$") 1\\,m\\BRg{} (\\%)", fmt.(enmo.pct_1mg; d = 1))
        fila("Correlación de MAD", fmt.(enmo.r_mad; d = 4))
        fila("ENMO de GGIR en \\emph{idle sleep} (m\\BRg)", fmt.(enmo.ggir_idle_enmo_medio_mg; d = 2))
        sub("No uso (ventanas de 15\\,min)")
        fila("Acuerdo (\\%) / $(raw"$\kappa$")", [fmt(a; d = 1) * " / " * fmt(b) for (a, b) in zip(nw.acuerdo_pct, nw.kappa)])
        fila("Horas de no uso: \\Brincadiano{} / GGIR", [fmt(a; d = 1) * " / " * fmt(b; d = 1) for (a, b) in zip(nw.nouso_nuestro_h, nw.nouso_ggir_h)])
        sub("Ritmo")
        hh(h) = isnan(h) ? "---" : Dates.format(Time(0) + Minute(round(Int, 60h) % 1440), "HH:MM")
        fila("Punto medio de L5: \\Brincadiano{} / GGIR", [hh(a) * " / " * hh(b) for (a, b) in zip(ritmo.L5_medio_nuestro_h, ritmo.L5_medio_ggir_h)])
        fila("Punto medio de M10: \\Brincadiano{} / GGIR", [hh(a) * " / " * hh(b) for (a, b) in zip(ritmo.M10_medio_nuestro_h, ritmo.M10_medio_ggir_h)])
        fila("IS: \\Brincadiano{} (ausente / reposo) / GGIR", [fmt(a) * " (" * fmt(c) * ") / " * fmt(b) for (a, c, b) in zip(ritmo.IS_nuestro, ritmo.IS_reposo_nuestro, ritmo.IS_ggir)])
        fila("IV: \\Brincadiano{} (ausente / reposo) / GGIR", [fmt(a) * " (" * fmt(c) * ") / " * fmt(b) for (a, c, b) in zip(ritmo.IV_nuestro, ritmo.IV_reposo_nuestro, ritmo.IV_ggir)])
        sub("Media circular: Julia frente a una implementación independiente (" * something(circ.motor[1], "numpy") * ")")
        fila("Diferencia en $(raw"$\bar\theta$") (h) / en $(raw"$R$")", [sci(a) * " / " * sci(b) for (a, b) in zip(circ.dif_θ_h, circ.dif_R)])
        println(io, "\\bottomrule")
        println(io, "\\end{tabularx}")
    end
end

if !("--confirmo" in ARGS)
    println("""Validación frente a GGIR NO ejecutada (falta --confirmo).
    Este script lanzaría GGIR (R) sobre data/raw/*.gt3x en un directorio temporal, extraería
    columnas de una lista blanca y compararía con Brincadiano. Solo con autorización.""")
else
    mkpath(GG); mkpath(OUT)
    # --sin-ggir: reutiliza la salida ya extraída en outputs/validation/ggir (depuración).
    if "--sin-ggir" in ARGS && isfile(joinpath(GG, "ggir_calibracion.csv"))
        println("GGIR: se reutiliza la salida extraída en outputs/validation/ggir")
    else
        println("GGIR 3.3.9: ejecutando (≈ 1 min)…"); correr_ggir()
    end
    cal_gg = CSV.read(joinpath(GG, "ggir_calibracion.csv"), DataFrame; missingstring = "NA")
    p2 = CSV.read(joinpath(GG, "ggir_part2_summary.csv"), DataFrame; missingstring = "NA")
    d2 = CSV.read(joinpath(GG, "ggir_part2_daysummary.csv"), DataFrame; missingstring = "NA")
    ids = sort(intersect(String.(cal_gg.id), String.(p2.id)))
    res = map(id -> comparar(id, cal_gg, p2, d2), ids)
    cal, enmo, nw, ritmo = (DataFrame(getindex.(res, k)) for k in 1:4)
    circ = circular_numpy(ids)
    for (nm, df) in (("calibracion", cal), ("enmo_mad", enmo), ("no_uso", nw), ("ritmo", ritmo), ("circular", circ))
        CSV.write(joinpath(OUT, "ggir_$(nm).csv"), df)
    end
    tablas_tex(cal, enmo, nw, ritmo, circ)
    for i in eachindex(ids)
        @printf("%-10s calibración: error %.1f→%.1f mg (GGIR %.1f→%.1f); ENMO 5 s: r = %.4f, |dif| media %.2f mg, %.1f %% < 1 mg\n",
                ids[i], cal.err_antes_nuestro[i], cal.err_despues_nuestro[i], cal.err_antes_ggir[i], cal.err_despues_ggir[i],
                enmo.r_enmo[i], enmo.mae_mg[i], enmo.pct_1mg[i])
        @printf("%-10s no uso: acuerdo %.1f %% (κ = %.2f), %.1f h frente a %.1f h; IS %.2f/%.2f (GGIR %.2f), IV %.2f/%.2f (GGIR %.2f); circular: |Δθ| %.1e h\n",
                ids[i], nw.acuerdo_pct[i], nw.kappa[i], nw.nouso_nuestro_h[i], nw.nouso_ggir_h[i],
                ritmo.IS_nuestro[i], ritmo.IS_reposo_nuestro[i], ritmo.IS_ggir[i], ritmo.IV_nuestro[i], ritmo.IV_reposo_nuestro[i],
                ritmo.IV_ggir[i], circ.dif_θ_h[i])
    end
    println("Tablas en ", relpath(OUT, ROOT))
end
