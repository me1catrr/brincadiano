# 08 · Periodo de reposo nocturno sin diario para un monitor en cadera.
#
# Qué hace
#   Estima, noche a noche (de 12:00 a 12:00), el periodo de reposo nocturno sin diario. Los
#   ángulos de 5 s siguen a GGIR 3.3.9 (10 Hz, mediana de 5 s). La ventana es la regla
#   HorAngle de GGIR para cadera: el eje longitudinal a menos de 60° de la horizontal
#   cuenta como tumbado, los tramos duran más de 30 min, los huecos de menos de 60 min
#   se unen y se queda el tramo más largo de la noche. El sueño dentro de la ventana
#   es el SIB de van Hees 2015 (eje z, 5°, 5 min). Adaptaciones de este proyecto: el
#   idle sleep hereda el último ángulo (queda escrito en angulos_5s.arrow); el no uso
#   nunca es «tumbado»; si hay al menos 2 h de no uso entre las 18:00 y las 10:00, la
#   ventana es SPT*, el tramo de no uso más largo (supuesto declarado, análogo al L5*).
#   Si ya existe el resumen del script 06, resta el ancla (L5 o L5*) del punto medio
#   de cada ventana. Sin ese CSV el sueño se calcula igual y la diferencia queda vacía.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/08_sueno.jl [ID ...]      # sin ID, quien tenga epochs_60s.arrow
#   julia --project=scripts scripts/08_sueno.jl --sintetico   # ilustración con datos simulados
#
# Depende de
#   01 (aceleración y huecos), 04 (qc_tramos.arrow, vía load_runs) y 05
#   (epochs_60s.arrow, para el L5 de cada noche y como señal de que el registro está
#   procesado). El 06 es opcional: outputs/circadian/tables/ritmo_resumen.csv.
#
# Lo usan
#   09, 10, 93, 94 y 96, que leen data/interim/<ID>/angulos_5s.arrow. El 10 y el
#   informe leen además sueno_<ID>_noches.csv. `make sync-assets` copia tablas y
#   figuras al capítulo 10. --sintetico no toca interim ni outputs/.
#
# Entrada
#   data/interim/<ID>/*.arrow, qc_tramos.arrow y epochs_60s.arrow
#   outputs/circadian/tables/ritmo_resumen.csv    opcional
#
# Salida (individual; no se versiona; las noches van numeradas, sin fecha de calendario)
#   data/interim/<ID>/angulos_5s.arrow            ángulo x/y/z, no uso e idle sleep heredado
#   outputs/sleep/tables/sueno_<ID>_noches.csv    una fila por noche
#   outputs/sleep/tables/sueno_resumen.{csv,tex}  una fila por registro; se reescribe entera
#   outputs/sleep/figures/<ID>_noches.{pdf,png}   de 12:00 a 12:00
# Salida versionada (--sintetico; dos noches inventadas, no de ningún participante)
#   report/assets/09_sueno/report/figures/sueno-ejemplo-sintetico.pdf
#
# Un registro sin noches completas se omite. Consola: solo el ID y los estadísticos
# de la ventana. El eje longitudinal es el de mayor mediana de |ángulo| en los epochs
# en uso (Sleep.longitudinal_axis). GGIR, en el script 93, lo elige de otra forma:
# correlación entre mitades de días; por eso aquel script lo fuerza al eje Y.

using Brincadiano.GT3X, Brincadiano.QC, Brincadiano.Metrics, Brincadiano.Sleep, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_T = joinpath(ROOT, "outputs", "sleep", "tables")
const OUT_F = joinpath(ROOT, "outputs", "sleep", "figures")
const RITMO = joinpath(ROOT, "outputs", "circadian", "tables", "ritmo_resumen.csv")
const EJES = ("X", "Y", "Z")

fmt(x; d = 1) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
hhmm(h) = isnan(h) ? "---" : Dates.format(Time(0) + Minute(round(Int, 60h) % 1440), "HH:MM")
# Diferencia circular en horas, en (−12, 12].
dif_h(a, b) = (d = mod(a - b, 24); d > 12 ? d - 24 : d)

function sueno(id)
    dir = joinpath(INTERIM, id)
    isfile(joinpath(dir, "epochs_60s.arrow")) || error("$id: faltan los epochs (ejecuta scripts/05_metricas.jl)")
    r = read_interim(dir); runs = load_runs(dir); fs = r.info.sample_rate
    ac = autocalibrate(r.acc, fs)
    off, sc = ac.aplicar ? (ac.offset, ac.scale) : ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    ep = epoch_angles(r.acc, r.gaps, fs, Second(5); offset = off, scale = sc)
    flag_nonwear!(ep, runs)
    epf = fill_idle(ep)
    save_angles(dir, epf)
    eje = longitudinal_axis(epf; mask = .!epf.nonwear)
    noches = sleep_nights(epf; axis = eje, nonwear = epf.nonwear,
                          rec_start = first(r.acc.time), rec_stop = last(r.acc.time))
    if nrow(noches) == 0
        println(rpad(id, 12), " sin noches completas: se omite")
        return nothing
    end
    e60 = load_epochs(dir, 60)
    noches.L5_noche_h = [night_l5(e60, n.noche; nonwear_as = n.metodo == "retirada" ? :rest : :missing)
                         for n in eachrow(noches)]
    noches.spt_menos_L5_h = [dif_h(a, b) for (a, b) in zip(noches.punto_medio_h, noches.L5_noche_h)]
    mkpath(OUT_T); mkpath(OUT_F)
    out = select(noches, Not([:inicio, :fin]))               # sin fechas de calendario
    out.noche = 1:nrow(out)
    CSV.write(joinpath(OUT_T, "sueno_$(id)_noches.csv"), out)

    s = sleep_summary(noches)
    ancla = NaN; nombre = ""
    if isfile(RITMO)
        rt = CSV.read(RITMO, DataFrame)
        k = findfirst(==(id), rt.id)
        if k !== nothing
            ancla = rt.L5_estimable[k] ? rt.L5_medio_h[k] : rt.L5_reposo_medio_h[k]
            nombre = rt.L5_estimable[k] ? "L5" : "L5*"
        end
    end
    ok = noches.metodo .!= "sin ventana"
    fila = (id = id, eje_longitudinal = EJES[eje], pct_idle_heredado = 100 * count(epf.imputado) / nrow(epf),
            s..., ancla = nombre, ancla_h = ancla,
            spt_menos_ancla_h = isnan(ancla) ? NaN : dif_h(mod(s.punto_medio_h, 24), ancla),
            spt_menos_L5_noche_h = any(ok) ? sum(noches.spt_menos_L5_h[ok]) / count(ok) : NaN)

    fig = plot_sleep_nights(epf, noches; axis = eje)
    guardar(joinpath(OUT_F, "$(id)_noches"), fig;
            procedencia = (; tipo = :datos, registros = [id], funcion = "plot_sleep_nights", entradas = ["data/interim/$id/angulos_5s.arrow"]))

    @printf("%-12s eje longitudinal %s; %d noches: %d por postura, %d por retirada (SPT*), %d sin periodo; ángulo heredado del idle sleep en el %.0f %% de los epochs\n",
            id, EJES[eje], s.noches, s.postura, s.retirada, s.sin_ventana, fila.pct_idle_heredado)
    @printf("%-12s periodo %s–%s (DE %.1f / %.1f h), punto medio %s (DE %.1f h), %.1f h; sueño %.1f h, eficiencia %.2f; jet lag social %.1f h\n",
            id, hhmm(s.inicio_h), hhmm(s.fin_h), s.inicio_de_h, s.fin_de_h, hhmm(s.punto_medio_h), s.punto_medio_de_h,
            s.duracion_h, s.sueno_h, s.eficiencia, s.jetlag_social_h)
    @printf("%-12s punto medio del periodo − %s del perfil: %+.1f h; − L5 de cada noche: %+.1f h (media)\n",
            id, nombre, fila.spt_menos_ancla_h, fila.spt_menos_L5_noche_h)
    return fila
end

function tabla_tex(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/08_sueno.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}l *{9}{>{\\raggedleft\\arraybackslash}X}@{}}")
    println(io, "\\toprule")
    println(io, "ID & Noches (postura / retirada) & Inicio (DE, h) & Fin (DE, h) & Punto medio (DE, h) & ",
            "Duración (h) & Sueño en el periodo (h) & Eficiencia & \\emph{Jet lag} social (h) & Punto medio \$-\$ L5 (h) \\\\")
    println(io, "\\midrule")
    for r in eachrow(df)
        println(io, replace(r.id, "_" => "\\_"), " & ", r.postura, " / ", r.retirada, " & ",
                hhmm(r.inicio_h), " (", fmt(r.inicio_de_h), ") & ", hhmm(r.fin_h), " (", fmt(r.fin_de_h), ") & ",
                hhmm(r.punto_medio_h), " (", fmt(r.punto_medio_de_h), ") & ", fmt(r.duracion_h), " & ",
                fmt(r.sueno_h), " & ", fmt(r.eficiencia; d = 2), " & ", fmt(r.jetlag_social_h), " & ",
                fmt(r.spt_menos_ancla_h), r.ancla == "L5*" ? "\\textsuperscript{*}" : "", " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

# Dos noches SIMULADAS (no proceden de ningún participante): la primera con el dispositivo
# puesto (tumbado de 21:40 a 07:50, con un despertar de 40 min a las 03:00 y vueltas en la
# cama), la segunda retirado de 21:15 a 08:05.
function sintetico()
    d0 = DateTime(2024, 3, 1)                               # viernes
    t = [d0 + Second(5k) for k in 0:(3 * 17280 - 1)]
    h = [Dates.value(x - d0) / 3_600_000 for x in t]
    rng = 12345
    ruido() = (rng = (1103515245 * rng + 12345) % 2^31; rng / 2^31 - 0.5)
    tumb1 = [(21 + 40 / 60 <= x < 31 + 50 / 60) && !(27 <= x < 27 + 40 / 60) for x in h]
    sentado = [(19 <= mod(x, 24) < 21) for x in h]          # sofá: tronco algo reclinado
    retir = [(24 + 21.25 <= x < 48 + 8 + 5 / 60) for x in h]
    angy = [tumb1[i] ? 8 + 6ruido() : sentado[i] ? 68 + 6ruido() : 84 + 4ruido() for i in eachindex(h)]
    angz = [tumb1[i] ? 40 * (floor(Int, h[i] / 0.75) % 3) + 2ruido() : 10ruido() for i in eachindex(h)]
    angy[retir] .= 2.0; angz[retir] .= 88.0                 # sobre la mesilla, inmóvil
    n = length(t)
    ep = DataFrame(start = t, n = fill(50, n), idle_frac = zeros(n), missing_frac = zeros(n),
                   anglex = zeros(n), angley = angy, anglez = angz, imputado = tumb1 .& iseven.(floor.(Int, 6 .* h)),
                   nonwear = retir)
    noches = sleep_nights(ep; axis = 2, nonwear = ep.nonwear)
    fig = plot_sleep_nights(ep, noches[1:2, :]; axis = 2)
    dir = joinpath(ROOT, "report", "assets", "09_sueno", "report", "figures")
    mkpath(dir)
    guardar(joinpath(dir, "sueno-ejemplo-sintetico"), fig; png = false,
            procedencia = (; tipo = :sintetico, funcion = "plot_sleep_nights"))
    println("  ", relpath(joinpath(dir, "sueno-ejemplo-sintetico.pdf"), ROOT))
    show(stdout, MIME("text/plain"), select(noches, :metodo, :inicio_h, :fin_h, :duracion_h, :sueno_h)); println()
end

if "--sintetico" in ARGS
    sintetico()
else
    ids = filter(a -> !startswith(a, "--"), ARGS)
    isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "epochs_60s.arrow")), readdir(INTERIM))))
    res = filter(!isnothing, map(sueno, ids))
    df = DataFrame(res)
    CSV.write(joinpath(OUT_T, "sueno_resumen.csv"), df)
    write(joinpath(OUT_T, "sueno_resumen.tex"), tabla_tex(df))
    println("Resumen en ", relpath(joinpath(OUT_T, "sueno_resumen.{csv,tex}"), ROOT))
end
