# 05 · Métricas ENMO y MAD en epochs de 1, 5 y 60 s.
#
# Qué hace
#   Agrega la aceleración en epochs alineados al reloj (1, 5 y 60 s) con ENMO y MAD,
#   cobertura y la marca de no uso heredada de los tramos del script 04. Antes estima
#   la autocalibración de van Hees (2014), la misma que usa GGIR y sin temperatura:
#   el offset y la escala se aplican a las métricas solo si la esfera de reposo está
#   bien cubierta y el error baja de 10 mg. El sesgo en reposo de la tabla se mide
#   antes de corregir, como diagnóstico. ENMO tiene dos columnas en el Arrow:
#   `enmo` usa solo las muestras presentes; `enmo_idle0` trata el idle sleep como 0
#   (la única imputación de esta serie, y va con nombre propio).
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/05_metricas.jl [ID ...]    # sin ID, todos los que tengan meta.arrow
#   julia --project=scripts scripts/05_metricas.jl --fixture   # ilustración versionable; no toca interim
#
# Depende de
#   01 (aceleración y huecos) y 04 (qc_tramos.arrow). Falla con un mensaje explícito
#   si falta el tramo de no uso.
#
# Lo usan
#   06 (epochs de 1 s para los episodios vigorosos y de 60 s para el perfil de 24 h),
#   08 (epochs de 60 s, para el L5 de cada noche), 09 (epochs de 60 s, para el estado
#   de la luz), 10 (epochs de 5 y 60 s), 92 (epochs de 5 s frente a GGIR) y 96
#   (epochs de 5 s, minutos por intensidad). El 95 no lee estos Arrow: repite el
#   cálculo, pero solo entra en la comparación quien ya tenga epochs_60s.arrow.
#   `make sync-assets` copia las tablas y la ventana al capítulo 8.
#
# Entrada
#   data/interim/<ID>/*.arrow del script 01
#   data/interim/<ID>/qc_tramos.arrow del script 04
#   Con --fixture: test/fixtures/legacy_sample.gt3x (MIT, no es de ningún participante)
#
# Salida (individual; no se versiona)
#   data/interim/<ID>/epochs_1s.arrow
#   data/interim/<ID>/epochs_5s.arrow
#   data/interim/<ID>/epochs_60s.arrow
#       columnas de tiempo, n, cobertura, idle_frac, enmo, enmo_idle0, mad, nonwear
#   outputs/metrics/tables/metricas_resumen.{csv,tex}     una fila por registro; se reescribe entera
#   outputs/metrics/tables/metricas_<ID>_dias.csv         uso y ENMO medio por día
#   outputs/metrics/figures/<ID>_epochs_ventana.{pdf,png} 5 min alrededor del minuto más intenso
# Salida versionada (--fixture)
#   report/assets/07_metricas/report/figures/metricas-ejemplo.pdf
#
# Un día válido, aquí y en el 06, son al menos 600 min de uso (las 10 h del 04).
# Consola: solo el ID, el sesgo, si la calibración se aplicó y el ENMO de los días válidos.

using Brincadiano.GT3X, Brincadiano.QC, Brincadiano.Metrics, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_T = joinpath(ROOT, "outputs", "metrics", "tables")
const OUT_F = joinpath(ROOT, "outputs", "metrics", "figures")
const EPOCHS = (1, 5, 60)
const MIN_USO = 600                   # minutos de uso para un día válido (10 h, como en H4)

# Coma decimal (convención del informe).
fmt(x; d = 1) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
pct(v, p) = (w = sort(filter(isfinite, v)); isempty(w) ? NaN : w[clamp(ceil(Int, p * length(w)), 1, length(w))])

function metricas(id)
    dir = joinpath(INTERIM, id)
    isfile(joinpath(dir, "qc_tramos.arrow")) || error("$id: falta qc_tramos.arrow (ejecuta scripts/04_qc.jl)")
    r = read_interim(dir); runs = load_runs(dir); fs = r.info.sample_rate
    # Autocalibración (van Hees 2014, como GGIR): solo se aplica si cumple sus criterios.
    ac = autocalibrate(r.acc, fs)
    off, sc = ac.aplicar ? (ac.offset, ac.scale) : ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    ep = Dict{Int,DataFrame}()
    for len in EPOCHS
        e = epoch_metrics(r.acc, r.gaps, fs, Second(len); offset = off, scale = sc)
        flag_nonwear!(e, runs)
        save_epochs(dir, e, len)
        ep[len] = e
    end
    cal = calibration_check(r.acc, fs)
    dias = daily_summary(ep[60]); dias.valido = dias.minutos_uso .>= MIN_USO
    mkpath(OUT_T); mkpath(OUT_F)
    CSV.write(joinpath(OUT_T, "metricas_$(id)_dias.csv"), dias)

    e1, e5, e60 = ep[1], ep[5], ep[60]
    uso5 = .!e5.nonwear .& (e5.n .> 0)
    validos = dias.dia[dias.valido]
    enmo_validos = isempty(validos) ? NaN : sum(dias.enmo_medio_mg[dias.valido]) / length(validos)
    fila = (id = id, sesgo_mg = cal.sesgo_mg, ventanas_reposo = cal.ventanas,
            cal_error_antes_mg = ac.error_antes_mg, cal_error_despues_mg = ac.error_despues_mg,
            cal_aplicada = ac.aplicar, offset_mg = round.(1000 .* ac.offset; digits = 1),
            escala = round.(ac.scale; digits = 4),
            enmo_medio_mg = enmo_validos, dias_validos = length(validos),
            mad5_mediana_mg = pct(e5.mad[uso5], 0.5), enmo1_p99_mg = pct(e1.enmo[e1.n .> 0], 0.99),
            enmo5_p99_mg = pct(e5.enmo[e5.n .> 0], 0.99), enmo60_p99_mg = pct(e60.enmo_idle0, 0.99),
            pct_epochs1_completos = 100 * count(e1.cov .== 1) / max(count(e1.n .> 0), 1))

    # Ventana de 5 min alrededor del minuto con más ENMO (con muestreo completo).
    full = (e60.cov .== 1) .& .!e60.nonwear
    if any(full)
        k = findall(full)[argmax(e60.enmo[full])]
        fig = plot_metrics_window(r, e60.start[k] - Minute(2), Minute(5); offset = off, scale = sc)
        guardar(joinpath(OUT_F, "$(id)_epochs_ventana"), fig;
                procedencia = (; tipo = :datos, registros = [id], funcion = "plot_metrics_window", entradas = ["data/interim/$id/"]))
    end
    @printf("%-12s sesgo en reposo %+.1f mg (%d ventanas); autocalibración %s: error %.1f → %.1f mg\n",
            id, cal.sesgo_mg, cal.ventanas, ac.aplicar ? "aplicada" : "NO aplicada", ac.error_antes_mg, ac.error_despues_mg)
    @printf("%-12s ENMO medio %.1f mg en %d días válidos; ENMO 1 s p99 %.0f mg\n",
            id, enmo_validos, length(validos), fila.enmo1_p99_mg)
    return fila
end

function tabla_tex(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/05_metricas.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}l *{8}{>{\\raggedleft\\arraybackslash}X}@{}}")
    println(io, "\\toprule")
    println(io, "ID & Sesgo en reposo (mg) & Error de calibración antes / después (mg) & ENMO medio (mg) & Días válidos & MAD 5\\,s mediana (mg) & ",
            "ENMO 1\\,s p99 (mg) & ENMO 5\\,s p99 (mg) & ENMO 60\\,s p99 (mg) \\\\")
    println(io, "\\midrule")
    for r in eachrow(df)
        println(io, replace(r.id, "_" => "\\_"), " & ", fmt(r.sesgo_mg), " & ",
                fmt(r.cal_error_antes_mg), " / ", fmt(r.cal_error_despues_mg), r.cal_aplicada ? "" : "\\textsuperscript{*}",
                " & ", fmt(r.enmo_medio_mg), " & ",
                r.dias_validos, " & ", fmt(r.mad5_mediana_mg), " & ", fmt(r.enmo1_p99_mg; d = 0), " & ",
                fmt(r.enmo5_p99_mg; d = 0), " & ", fmt(r.enmo60_p99_mg; d = 0), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

function fixture()
    r = read_gt3x(joinpath(ROOT, "test", "fixtures", "legacy_sample.gt3x"))
    fig = plot_metrics_window(r, first(r.acc.time), Second(80))
    dir = joinpath(ROOT, "report", "assets", "07_metricas", "report", "figures")
    mkpath(dir)
    guardar(joinpath(dir, "metricas-ejemplo"), fig; png = false,
            procedencia = (; tipo = :fixture, funcion = "plot_metrics_window", entradas = ["test/fixtures/legacy_sample.gt3x"]))
    println("  ", relpath(joinpath(dir, "metricas-ejemplo.pdf"), ROOT))
end

if "--fixture" in ARGS
    fixture()
else
    ids = filter(a -> !startswith(a, "--"), ARGS)
    isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "meta.arrow")), readdir(INTERIM))))
    res = DataFrame(map(metricas, ids))
    CSV.write(joinpath(OUT_T, "metricas_resumen.csv"), res)
    write(joinpath(OUT_T, "metricas_resumen.tex"), tabla_tex(res))
    println("Resumen en ", relpath(joinpath(OUT_T, "metricas_resumen.{csv,tex}"), ROOT))
end
