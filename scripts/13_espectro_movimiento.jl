# 13 · Espectro del movimiento (hito H21, escala de la señal a 30 Hz).
#
# Qué hace
#   En ventanas de 5 s alineadas al reloj, completas (sin huecos) y en uso, de los días
#   válidos (≥ 600 min), la densidad espectral de la norma del vector calibrada
#   (Spectral.window_spectra: Hann, sin rellenar nada) y sus rasgos: frecuencia dominante
#   (0,5–5 Hz: la cadencia si hay marcha o carrera), fracción de potencia en 0,25–5, 5–10 y
#   10–15 Hz, y entropía espectral. Con ello:
#     · PSD media por clase de intensidad (ENMO de 5 s con los umbrales de Hildebrand) y el
#       suelo de ruido (1 % de ventanas más quietas) frente a la cuantización (1/256 𝑔);
#     · cadencia por clase;
#     · QC espectral de la actividad vigorosa: la ventana de más potencia de cada episodio
#       del script 06 y la de cada segundo vigoroso aislado (≥ 464,6 mg fuera de episodios,
#       cap. 11); «vibración probable» si más de la mitad de la potencia está en 10–15 Hz.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/13_espectro_movimiento.jl [ID ...]
#
# Depende de
#   01 (acc.arrow), 05 (epochs_1s, epochs_5s, epochs_60s), 06 (episodios_<ID>.csv).
#
# Salida (individual; no se versiona)
#   outputs/metrics/tables/espectro_movimiento.{csv,tex}     rasgos por clase y registro
#   outputs/metrics/figures/psd_por_intensidad.{pdf,png}     PSD media por clase
#   outputs/vigorous/figures/cadencia.{pdf,png}              frecuencia dominante por clase
#   outputs/vigorous/tables/qc_espectral.{csv,tex}           episodios y segundos aislados
#   outputs/vigorous/figures/qc_espectral.{pdf,png}          episodios frente a segundos aislados
#
# Consola: solo el ID y las cifras resumen. Detalle en los capítulos 7 y 8.

using Brincadiano.GT3X, Brincadiano.Metrics, Brincadiano.Spectral, Brincadiano.Viz
using Brincadiano.Rutas: Rutas
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT = joinpath(ROOT, "outputs")
const MIN_USO = 600
const UMBRALES = (63.3, 142.6, 464.6)          # Hildebrand, cadera (mg)
const CLASES = ("sedentaria", "ligera", "moderada", "vigorosa")
const UMBRAL_VIBRACION = 0.5                   # fracción de potencia en 10–15 Hz
const Q_MG = 1000 / 256                        # cuantización (mg)

fmt(x; d = 2) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")
cuantil(v, p) = isempty(v) ? NaN : (s = sort(collect(v)); s[clamp(round(Int, p * (length(s) - 1)) + 1, 1, length(s))])
mediana(v) = cuantil(v, 0.5)
clase(e) = e < UMBRALES[1] ? CLASES[1] : e < UMBRALES[2] ? CLASES[2] : e < UMBRALES[3] ? CLASES[3] : CLASES[4]

function movimiento(id)
    dir = joinpath(INTERIM, id)
    fep = joinpath(OUT, "timing", "tables", "episodios_$(id).csv")
    (isfile(joinpath(dir, "epochs_60s.arrow")) && isfile(fep)) || return nothing
    e60 = load_epochs(dir, 60); dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    isempty(validos) && return nothing
    r = read_interim(dir); fs = r.info.sample_rate
    ac = autocalibrate(r.acc, fs)
    off, sc = ac.aplicar ? (ac.offset, ac.scale) : ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    ws = window_spectra(r.acc, fs, Second(5); offset = off, scale = sc, psd = true)
    fr = ws.frecuencias
    e5 = load_epochs(dir, 5)
    d = ws.rasgos; d.fila = 1:nrow(d)
    d = innerjoin(d, select(e5, :start, :enmo, :nonwear), on = :start)
    d = d[[Date(t) in validos for t in d.start] .& .!d.nonwear .& isfinite.(d.enmo), :]
    d.clase = clase.(d.enmo)
    posibles = 17280 * length(validos)

    # PSD media por clase y suelo de ruido.
    psd_clase = [(c, vec(sum(ws.psd[d.fila[d.clase .== c], :]; dims = 1)) ./ max(count(==(c), d.clase), 1))
                 for c in CLASES if any(==(c), d.clase)]
    quietas = d.fila[d.potencia .<= cuantil(d.potencia, 0.01)]
    suelo = vec(sum(ws.psd[quietas, :]; dims = 1)) ./ length(quietas)
    k = fr .>= 0.25
    suelo_nivel = sum(suelo[k]) / count(k)                     # m𝑔²/Hz, media en 0,25–15 Hz
    filas = NamedTuple[]
    for c in CLASES
        g = d[d.clase .== c, :]
        push!(filas, (id = id, clase = c, ventanas = nrow(g), pct = 100nrow(g) / nrow(d),
                      f_dom = mediana(g.f_dom), f_dom_q1 = cuantil(g.f_dom, 0.25), f_dom_q3 = cuantil(g.f_dom, 0.75),
                      frac_1 = mediana(g.frac_1), frac_2 = mediana(g.frac_2), frac_3 = mediana(g.frac_3),
                      entropia = mediana(g.entropia), potencia = mediana(g.potencia),
                      completas_pct = 100nrow(d) / posibles, suelo_mg2 = mediana(d.potencia[d.fila .∈ Ref(Set(quietas))]),
                      suelo_psd = suelo_nivel))
    end

    # QC espectral de la actividad vigorosa.
    ep = CSV.read(fep, DataFrame; types = Dict(:inicio => DateTime, :fin => DateTime, :medio => DateTime))
    idx = Dict(t => i for (i, t) in enumerate(d.start))
    ve = map(eachrow(ep)) do e
        ks = [idx[t] for t in floor(e.inicio, Second(5)):Second(5):e.fin - Millisecond(1) if haskey(idx, t)]
        isempty(ks) ? 0 : ks[argmax(d.potencia[ks])]
    end
    e1 = load_epochs(dir, 1)
    pos = Dict(t => i for (i, t) in enumerate(e1.start))
    en_ep = falses(nrow(e1))
    for e in eachrow(ep), t in e.inicio:Second(1):e.fin - Second(1)
        i = get(pos, t, 0); i > 0 && (en_ep[i] = true)
    end
    vig = [(Date(t) in validos) && !nw && isfinite(x) && x >= UMBRALES[3] for (t, nw, x) in zip(e1.start, e1.nonwear, e1.enmo)]
    ais = findall(vig .& .!en_ep)
    va = [get(idx, floor(e1.start[i], Second(5)), 0) for i in ais]
    epw = d[ve[ve .> 0], :]; aisw = d[va[va .> 0], :]
    h_ais = hour.(e1.start[ais])
    qc = (id = id, episodios = nrow(ep), con_ventana = count(>(0), ve),
          vibracion = count(epw.frac_3 .> UMBRAL_VIBRACION), frac_3_max = isempty(epw.frac_3) ? NaN : maximum(epw.frac_3),
          f_dom_ep = mediana(epw.f_dom), f_dom_ep_q1 = cuantil(epw.f_dom, 0.25), f_dom_ep_q3 = cuantil(epw.f_dom, 0.75),
          segundos_vig = count(vig), aislados = length(ais), aislados_pct = 100length(ais) / max(count(vig), 1),
          f_dom_ais = mediana(aisw.f_dom), frac_1_ais = mediana(aisw.frac_1), frac_1_ep = mediana(epw.frac_1),
          vibracion_ais = count(aisw.frac_3 .> UMBRAL_VIBRACION),
          enmo5_ais = mediana(aisw.enmo), enmo5_ep = mediana(epw.enmo),
          hora_pico_ais = argmax([count(==(h), h_ais) for h in 0:23]) - 1,
          pct_hora_pico_ais = 100maximum([count(==(h), h_ais) for h in 0:23]) / max(length(ais), 1))

    @printf("%-12s %.0f %% de ventanas completas; cadencia (mediana) ligera %.2f, moderada %.2f, vigorosa %.2f Hz; suelo %.1f mg² (cuantización %.1f)\n",
            id, filas[1].completas_pct, filas[2].f_dom, filas[3].f_dom, filas[4].f_dom, filas[1].suelo_mg2, Q_MG^2 / 12)
    @printf("%-12s episodios: %d, con vibración %d (máx. %.2f en 10–15 Hz), cadencia %.2f Hz; segundos aislados %d (%.0f %%), cadencia %.2f Hz, %d con vibración\n",
            id, qc.episodios, qc.vibracion, qc.frac_3_max, qc.f_dom_ep, qc.aislados, qc.aislados_pct, qc.f_dom_ais, qc.vibracion_ais)
    clases_f = [(c, d.f_dom[d.clase .== c]) for c in CLASES[2:4]]
    return DataFrame(filas), qc, (id, psd_clase), suelo, (id, clases_f), (id, epw, aisw), fr
end

function tex_movimiento(df)
    ids = unique(df.id)
    io = IOBuffer()
    println(io, "% Generado por scripts/13_espectro_movimiento.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}X *{$(4length(ids))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(["\\multicolumn{4}{c}{\\texttt{" * replace(id, "_" => "\\_") * "}}" for id in ids], " & "), " \\\\")
    println(io, join(["\\cmidrule(l){$(2 + 4(j - 1))-$(5 + 4(j - 1))}" for j in eachindex(ids)]))
    println(io, "Clase", repeat(" & Ventanas (\\%) & \$f_\\mathrm{dom}\$ (Hz) & 0,25--5\\,Hz (\\%) & Entropía", length(ids)), " \\\\")
    println(io, "\\midrule")
    for c in CLASES
        print(io, uppercasefirst(c))
        for id in ids
            r = df[(df.id .== id) .& (df.clase .== c), :][1, :]
            print(io, " & ", fmt(r.pct; d = 1), " & ", fmt(r.f_dom; d = 1), " (", fmt(r.f_dom_q1; d = 1), "--", fmt(r.f_dom_q3; d = 1), ")",
                  " & ", fmt(100r.frac_1; d = 0), " & ", fmt(r.entropia))
        end
        println(io, " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

function tex_qc(qc)
    io = IOBuffer()
    println(io, "% Generado por scripts/13_espectro_movimiento.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}X *{$(nrow(qc))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(["\\texttt{" * replace(id, "_" => "\\_") * "}" for id in qc.id], " & "), " \\\\")
    println(io, "\\midrule")
    filas = [
        ("Episodios vigorosos / con vibración probable", r -> "$(r.episodios) / $(r.vibracion)"),
        ("Máxima potencia en 10--15\\,Hz de un episodio", r -> fmt(r.frac_3_max)),
        ("Cadencia de los episodios (Hz, mediana y RIC)", r -> fmt(r.f_dom_ep; d = 1) * " (" * fmt(r.f_dom_ep_q1; d = 1) * "--" * fmt(r.f_dom_ep_q3; d = 1) * ")"),
        ("Segundos vigorosos / aislados (\\%)", r -> "$(r.segundos_vig) / " * fmt(r.aislados_pct; d = 0)),
        ("Cadencia de los aislados (Hz)", r -> fmt(r.f_dom_ais; d = 1)),
        ("Potencia en 0,25--5\\,Hz: episodios / aislados (\\%)", r -> fmt(100r.frac_1_ep; d = 0) * " / " * fmt(100r.frac_1_ais; d = 0)),
        ("ENMO de 5\\,s: episodios / aislados (m\\BRg)", r -> fmt(r.enmo5_ep; d = 0) * " / " * fmt(r.enmo5_ais; d = 0)),
        ("Aislados con vibración probable", r -> string(r.vibracion_ais)),
        ("Hora con más aislados (\\% de ellos)", r -> "$(r.hora_pico_ais):00 (" * fmt(r.pct_hora_pico_ais; d = 0) * ")"),
    ]
    for (nombre, f) in filas
        println(io, nombre, " & ", join([f(r) for r in eachrow(qc)], " & "), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

# Figuras y tablas .tex del informe: solo los registros de prueba (MELIPOP y pruebas). Las figuras
# superponen registros (≤ 4 colores) y no deben mezclar proyectos; los CSV llevan todos los IDs.
informe(id) = Rutas.clasificar(id).proyecto in ("MELIPOP", "pruebas")
ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "epochs_60s.arrow")), readdir(INTERIM))))
res = filter(!isnothing, map(movimiento, ids))
for d in ("metrics", "vigorous"), k in ("tables", "figures"); mkpath(joinpath(OUT, d, k)); end
fr = res[1][7]
df = reduce(vcat, [r[1] for r in res]); qc = DataFrame([r[2] for r in res])
CSV.write(joinpath(OUT, "metrics", "tables", "espectro_movimiento.csv"), df)
write(joinpath(OUT, "metrics", "tables", "espectro_movimiento.tex"), tex_movimiento(df[informe.(df.id), :]))
CSV.write(joinpath(OUT, "vigorous", "tables", "qc_espectral.csv"), qc)
write(joinpath(OUT, "vigorous", "tables", "qc_espectral.tex"), tex_qc(qc[informe.(qc.id), :]))
res = filter(r -> informe(r[3][1]), res)
regs = [r[3][1] for r in res]
guardar(joinpath(OUT, "metrics", "figures", "psd_por_intensidad"),
        plot_psd_classes([r[3] for r in res], fr; suelo = [r[4] for r in res], cuantizacion = Q_MG^2 / 12 / (fr[end]));
        procedencia = (; tipo = :datos, registros = regs, funcion = "plot_psd_classes", entradas = ["data/interim/<ID>/acc.arrow", "data/interim/<ID>/epochs_60s.arrow"]))
guardar(joinpath(OUT, "vigorous", "figures", "cadencia"), plot_cadence([r[5] for r in res], fr);
        procedencia = (; tipo = :datos, registros = regs, funcion = "plot_cadence", entradas = ["data/interim/<ID>/acc.arrow", "data/interim/<ID>/epochs_60s.arrow"]))
guardar(joinpath(OUT, "vigorous", "figures", "qc_espectral"), plot_episode_spectral_qc([r[6] for r in res]; umbral = UMBRAL_VIBRACION);
        procedencia = (; tipo = :datos, registros = regs, funcion = "plot_episode_spectral_qc", entradas = ["data/interim/<ID>/acc.arrow", "outputs/timing/tables/episodios_<ID>.csv"]))
println("Espectro del movimiento en outputs/metrics y outputs/vigorous")
