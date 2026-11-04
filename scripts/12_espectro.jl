# 12 · Análisis espectral del ritmo de reposo-actividad (hito H20, escala circadiana).
#
# Qué hace
#   Sobre los epochs de 60 s de los días válidos (≥ 600 min de uso), con el no uso como
#   ausente o como reposo según la misma regla que el ancla del script 06 (L5 estimable →
#   ausente; si no, reposo):
#     · periodograma de Lomb–Scargle generalizado de ENMO por minuto (0,2–12 ciclos/día),
#       sin rellenar huecos, y la ventana espectral del muestreo (siempre con el no uso
#       como ausente: muestra qué introduce el patrón de retiradas);
#     · periodograma χ² de Sokolove–Bushell con bins de 10 min (periodos de 2 a 36 h);
#     · IS repartido en los armónicos de 24 h (Parseval sobre el perfil horario).
#   Y sobre los epochs de 1 s, el timing respecto al ancla: los episodios del script 06
#   (≥ 464,6 mg) y los momentos armónicos de la dosis (Spectral.harmonic_timing: segundos
#   ≥ 464,6 y ≥ 142,6 mg sin definir episodios, y toda la actividad), con y sin corregir la
#   cobertura desigual del reloj (días incompletos, no uso; Spectral.clock_coverage).
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/12_espectro.jl [ID ...]
#
# Depende de
#   05: epochs_1s.arrow y epochs_60s.arrow. 06: ritmo_resumen.csv (ancla y regla del no
#   uso). Los episodios se recalculan con las mismas reglas del 06.
#
# Salida (individual; no se versiona)
#   outputs/circadian/tables/espectro_resumen.{csv,tex}   una fila por registro
#   outputs/circadian/tables/periodograma_<ID>.csv        Lomb–Scargle y ventana por frecuencia
#   outputs/circadian/figures/espectro_ritmo.{pdf,png}    periodogramas de los registros
#   outputs/circadian/figures/perfil_armonicos.{pdf,png}  perfil horario y su reconstrucción
#   outputs/timing/tables/timing_armonicos.{csv,tex}      θ y R: episodios, armónicos, cobertura
#
# Consola: solo el ID y las cifras resumen. Detalle en los capítulos 10 y 11.

using Brincadiano.Metrics, Brincadiano.Circadian, Brincadiano.Spectral, Brincadiano.Viz
using Brincadiano.Rutas: Rutas
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_C = joinpath(ROOT, "outputs", "circadian")
const OUT_T = joinpath(ROOT, "outputs", "timing")
const RITMO = joinpath(OUT_C, "tables", "ritmo_resumen.csv")
const MIN_USO = 600
const BIN_MIN = 10                     # bin del periodograma χ² (min)
const FMAX = 12.0                      # ciclos/día (periodo mínimo de 2 h)

fmt(x; d = 2) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
hhmm(h) = isnan(h) ? "---" : Dates.format(Time(0) + Minute(round(Int, 60h) % 1440), "HH:MM")
dif_h(a, b) = (d = mod(a - b, 24); d > 12 ? d - 24 : d)      # diferencia circular, (−12, 12]

# Serie por minuto de los días válidos: instantes (h desde la medianoche del primer día
# válido, centro del minuto) y ENMO con el no uso como NaN (:missing) o 0 (:rest).
function serie_minutos(e60, validos; nonwear_as)
    d0 = DateTime(minimum(validos))
    k = [Date(t) in validos for t in e60.start]
    t = [Dates.value(s - d0) / 3_600_000 + 1 / 120 for s in e60.start[k]]
    y = [nw ? (nonwear_as === :rest ? 0.0 : NaN) : v for (nw, v) in zip(e60.nonwear[k], e60.enmo_idle0[k])]
    return t, y
end

# Bins de BIN_MIN min por día válido (los días se encadenan: cada uno aporta 1440/BIN_MIN bins,
# así la posición en el ciclo de 24 h se conserva aunque falte algún día). NaN si el bin tiene
# menos de la mitad de minutos con dato.
function serie_bins(e60, validos; nonwear_as)
    dias = sort(collect(validos)); nb = 1440 ÷ BIN_MIN
    pos = Dict(d => i for (i, d) in enumerate(dias))
    s = zeros(nb * length(dias)); c = zeros(Int, length(s))
    for (t, nw, v) in zip(e60.start, e60.nonwear, e60.enmo_idle0)
        i = get(pos, Date(t), 0); i == 0 && continue
        b = (i - 1) * nb + (hour(t) * 60 + minute(t)) ÷ BIN_MIN + 1
        if nw
            nonwear_as === :rest && (c[b] += 1)
        elseif isfinite(v)
            s[b] += v; c[b] += 1
        end
    end
    return [c[b] >= BIN_MIN / 2 ? s[b] / c[b] : NaN for b in eachindex(s)]
end

function espectro(id, ritmo)
    dir = joinpath(INTERIM, id)
    k = findfirst(==(id), ritmo.id)
    (k === nothing || !isfile(joinpath(dir, "epochs_60s.arrow"))) && return nothing
    e1 = load_epochs(dir, 1); e60 = load_epochs(dir, 60)
    dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    isempty(validos) && return nothing
    estimable = ritmo.L5_estimable[k]
    nw_as = estimable ? :missing : :rest
    ancla = estimable ? ritmo.L5_medio_h[k] : ritmo.L5_reposo_medio_h[k]

    t, y = serie_minutos(e60, validos; nonwear_as = nw_as)
    T = maximum(t) - minimum(t)
    fr = freq_grid(T; fmin = 0.2, fmax = FMAX)
    g = gls_periodogram(t, y, fr)
    tm, ym = serie_minutos(e60, validos; nonwear_as = :missing)
    w = spectral_window(tm[isfinite.(ym)], fr)
    CSV.write(joinpath(OUT_C, "tables", "periodograma_$(id).csv"),
              DataFrame(frecuencia = g.frecuencia, potencia = g.potencia, amplitud = g.amplitud, ventana = w))
    gh = gls_periodogram(t, y, [1.0, 2.0, 3.0, 4.0])
    wh = spectral_window(tm[isfinite.(ym)], [1.0])
    q = chi2_periodogram(serie_bins(e60, validos; nonwear_as = nw_as), BIN_MIN / 60; periods_h = 2:BIN_MIN/60:36)
    qi = argmax(q.Qp .- q.Qp_crit)
    h = harmonic_is(e60, validos; nonwear_as = nw_as)
    a = h.armonicos

    fila = (id = id, no_uso = String(nw_as), dias = length(validos), minutos = count(isfinite, y),
            pico_cpd = g.frecuencia[argmax(g.potencia)], pico_periodo_h = g.periodo_h[argmax(g.potencia)],
            potencia_24h = gh.potencia[1], amplitud_24h = gh.amplitud[1], acrofase_24h = mod(gh.acrofase_h[1], 24),
            potencia_12h = gh.potencia[2], potencia_8h = gh.potencia[3], potencia_6h = gh.potencia[4],
            ventana_24h = wh[1], chi2_pico_h = q.periodo_h[qi], chi2_rel_24h = q.Qp_rel[q.periodo_h .≈ 24][1],
            IS = h.IS, IS_armonicos = h.IS_armonicos, IS_1 = a.IS_k[1], IS_2 = a.IS_k[2], IS_3 = a.IS_k[3],
            IS_4 = a.IS_k[4], IS_1a4 = sum(a.IS_k[1:4]), acrofase_1 = a.acrofase_h[1], amplitud_1 = a.amplitud[1],
            ancla_h = ancla)

    # Timing: episodios (cap. 11) frente a momentos armónicos, con y sin corregir la
    # cobertura desigual del reloj (días incompletos, no uso).
    cob = clock_coverage(e1, validos; nonwear_as = nw_as)
    ep = vigorous_episodes(e1, validos; threshold = 464.6)
    ep.peso = coverage_weights(ep, cob)
    filas_t = NamedTuple[]
    ref = timing(ep, ancla; peso = :peso)
    function fila_t(nombre, θ, R, R2)
        push!(filas_t, (id = id, variante = nombre, θ1_h = θ, R1 = R, R2 = R2, Δθ_h = dif_h(θ, ref.θ_h)))
    end
    fila_t("Episodios ≥ 464,6 mg (referencia)", ref.θ_h, ref.R, ref.R2)
    t0 = timing(ep, ancla); fila_t("Episodios, sin corregir la cobertura", t0.θ_h, t0.R, t0.R2)
    for (nombre, thr, corr) in (("Segundos ≥ 464,6 mg", 464.6, false), ("Segundos ≥ 464,6 mg, corregidos", 464.6, true),
                                ("Segundos ≥ 142,6 mg, corregidos", 142.6, true),
                                ("Toda la actividad", 0.0, false), ("Toda la actividad, corregida", 0.0, true))
        r = harmonic_timing(e1, validos, ancla; threshold = thr, K = 2, cobertura = corr, nonwear_as = nw_as).armonicos
        fila_t(nombre, r.θ_h[1], r.R[1], r.R[2])
    end
    ft = DataFrame(filas_t)
    cobertura_min, cobertura_max = extrema(cob)

    @printf("%-12s %d días, no uso %s: pico LS %.2f cpd (%.1f h), 24 h explica %.0f %% (ventana %.2f); χ² pico %.1f h; IS %.2f = %.0f %% en 24 h + %.0f %% en 12–6 h\n",
            id, length(validos), nw_as == :rest ? "= reposo" : "ausente", fila.pico_cpd, fila.pico_periodo_h,
            100fila.potencia_24h, fila.ventana_24h, fila.chi2_pico_h, h.IS, 100fila.IS_1 / h.IS,
            100(fila.IS_1a4 - fila.IS_1) / h.IS)
    @printf("%-12s cobertura del reloj de %d a %d días; θ episodios %.1f h (sin corregir %+.2f h); segundos corregidos %+.2f h; acrofase total corregida %s (Lomb–Scargle %s)\n",
            id, cobertura_min, cobertura_max, ref.θ_h, ft.Δθ_h[2], ft.Δθ_h[4], hhmm(mod(ft.θ1_h[7] + ancla, 24)), hhmm(fila.acrofase_24h))
    return fila, ft, (id, g, w, q), (id, h)
end

function tex_espectro(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/12_espectro.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}X *{$(nrow(df))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(["\\texttt{" * replace(r.id, "_" => "\\_") * "}" * (r.no_uso == "rest" ? "\\textsuperscript{*}" : "") for r in eachrow(df)], " & "), " \\\\")
    println(io, "\\midrule")
    filas = [
        ("Días válidos / minutos con dato", r -> "$(r.dias) / $(r.minutos)"),
        ("Pico de Lomb--Scargle (periodo, h)", r -> fmt(r.pico_periodo_h; d = 1)),
        ("Varianza explicada por 24 / 12 / 8 / 6\\,h (\\%)", r -> join(fmt.(100 .* [r.potencia_24h, r.potencia_12h, r.potencia_8h, r.potencia_6h]; d = 0), " / ")),
        ("Amplitud (m\\BRg) y acrofase de 24\\,h", r -> fmt(r.amplitud_24h; d = 1) * " · " * hhmm(r.acrofase_24h)),
        ("Ventana espectral en 1 ciclo/día", r -> fmt(r.ventana_24h)),
        ("Pico del periodograma \$\\chi^2\$ (h)", r -> fmt(r.chi2_pico_h; d = 1)),
        ("\$Q_p\$/umbral a 24\\,h", r -> fmt(r.chi2_rel_24h; d = 1)),
        ("IS", r -> fmt(r.IS)),
        ("IS en el armónico de 24\\,h (\\%)", r -> fmt(100r.IS_1 / r.IS; d = 0)),
        ("IS en los de 12, 8 y 6\\,h (\\%)", r -> fmt(100(r.IS_1a4 - r.IS_1) / r.IS; d = 0)),
    ]
    for (nombre, f) in filas
        println(io, nombre, " & ", join([f(r) for r in eachrow(df)], " & "), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

function tex_armonicos(df)
    ids = unique(df.id)
    io = IOBuffer()
    println(io, "% Generado por scripts/12_espectro.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}X *{$(4length(ids))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(["\\multicolumn{4}{c}{\\texttt{" * replace(id, "_" => "\\_") * "}}" for id in ids], " & "), " \\\\")
    println(io, join(["\\cmidrule(l){$(2 + 4(j - 1))-$(5 + 4(j - 1))}" for j in eachindex(ids)]))
    println(io, "Dosis", repeat(" & \$\\theta_1\$ (h) & \$\\Delta\\theta_1\$ (h) & \$R_1\$ & \$R_2\$", length(ids)), " \\\\")
    println(io, "\\midrule")
    for v in unique(df.variante)
        print(io, replace(v, r"(\d),(\d)" => s"\1{,}\2", " mg" => "\\,mg", "≥" => "\$\\geq\$"))
        for id in ids
            r = df[(df.id .== id) .& (df.variante .== v), :][1, :]
            ref = startswith(v, "Episodios ≥")
            print(io, " & ", fmt(r.θ1_h; d = 1), " & ", ref ? "---" : (r.Δθ_h > 0.05 ? "+" : "") * fmt(r.Δθ_h; d = 1),
                  " & ", fmt(r.R1), " & ", fmt(r.R2))
        end
        println(io, " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

# Figuras y tablas .tex del informe: solo los registros de prueba (MELIPOP y pruebas). Las figuras
# superponen registros (≤ 4 colores) y no deben mezclar proyectos; los CSV llevan todos los IDs.
informe(id) = Rutas.clasificar(id).proyecto in ("MELIPOP", "pruebas")

isfile(RITMO) || error("falta $(RITMO) (ejecuta scripts/06_circadiano.jl)")
ritmo = CSV.read(RITMO, DataFrame)
ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "epochs_60s.arrow")), readdir(INTERIM))))
for d in (OUT_C, OUT_T), k in ("tables", "figures"); mkpath(joinpath(d, k)); end
res = filter(!isnothing, map(id -> espectro(id, ritmo), ids))
df = DataFrame([r[1] for r in res]); dft = reduce(vcat, [r[2] for r in res])
CSV.write(joinpath(OUT_C, "tables", "espectro_resumen.csv"), df)
write(joinpath(OUT_C, "tables", "espectro_resumen.tex"), tex_espectro(df[informe.(df.id), :]))
CSV.write(joinpath(OUT_T, "tables", "timing_armonicos.csv"), dft)
write(joinpath(OUT_T, "tables", "timing_armonicos.tex"), tex_armonicos(dft[informe.(dft.id), :]))
res = filter(r -> informe(r[3][1]), res)
regs = [r[3][1] for r in res]
guardar(joinpath(OUT_C, "figures", "espectro_ritmo"), plot_rhythm_spectrum([r[3] for r in res]);
        procedencia = (; tipo = :datos, registros = regs, funcion = "plot_rhythm_spectrum", entradas = ["data/interim/<ID>/epochs_60s.arrow"]))
guardar(joinpath(OUT_C, "figures", "perfil_armonicos"), plot_harmonic_profile([r[4] for r in res]);
        procedencia = (; tipo = :datos, registros = regs, funcion = "plot_harmonic_profile", entradas = ["data/interim/<ID>/epochs_60s.arrow"]))
println("Espectro en outputs/circadian/{tables,figures} y outputs/timing/tables/timing_armonicos.*")
