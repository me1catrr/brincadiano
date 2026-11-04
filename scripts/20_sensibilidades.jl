# 20 · Sensibilidades pendientes (H36): sueño, episodios y submuestreo de días.
#
# Qué hace
#   Repite, cambiando una sola cosa cada vez, tres piezas que el informe declaraba pendientes:
#     A · Sueño: umbral de postura (45°, 60° referencia, 75°) y no uso que convierte la noche en
#         SPT* (1 h, 2 h referencia, 3 h). Por variante: noches por método, punto medio y
#         duración medios de la SPT.
#     B · Episodios vigorosos: tolerancia de pausa (0, 2 referencia, 5 s), duración mínima (1,
#         3 referencia, 5, 10 s) y epochs de 5 s en vez de 1 s. Por variante: n, θ, R, Δθ y |Δz|.
#     C · Submuestreo de días: todos los subconjuntos de k días válidos (k = 2 … n−1, ancla de
#         referencia fija); dispersión de θ, R, IS e IV frente al registro completo.
#   Control: la variante de referencia reproduce θ y R de 06 (mismo ancla, umbral y pesos).
#   No estima la fiabilidad entre niños (CCI): con tres registros no hay varianza entre niños;
#   `Circadian.timing_reliability` la calculará con la cohorte.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/20_sensibilidades.jl [ID ...]
#
# Depende de
#   05 (epochs 1, 5 y 60 s), 06 (ritmo_resumen.csv: ancla de referencia) y 08 (ángulos).
#
# Salida (individual; no se versiona)
#   outputs/sensibilidad/tables/sens_sueno.csv, sens_episodios.csv, sens_dias.csv
#   outputs/sensibilidad/tables/sens_dias_resumen.csv   por k: mediana y rango de |Δθ|, IS, IV
#
# Consola: solo el ID y cifras resumen.

using Brincadiano.Metrics, Brincadiano.Circadian, Brincadiano.Sleep
using CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT = joinpath(ROOT, "outputs", "sensibilidad", "tables")
const RITMO = joinpath(ROOT, "outputs", "circadian", "tables", "ritmo_resumen.csv")
const MIN_USO = 600
const UMBRAL_VIG = 464.6

dif_h(a, b) = (d = mod(a - b, 24); d > 12 ? d - 24 : d)
media_circular_h(h) = phase_to_hours(circular_mean(2π .* h ./ 24).θ)
combinaciones(v, k) = k == 0 ? [eltype(v)[]] :
    isempty(v) ? Vector{eltype(v)}[] :
    vcat([vcat([v[1]], c) for c in combinaciones(v[2:end], k - 1)], combinaciones(v[2:end], k))

function procesar(id, ritmo)
    dir = joinpath(INTERIM, id)
    k = findfirst(==(id), ritmo.id)
    (k === nothing || !isfile(joinpath(dir, "angulos_5s.arrow"))) && return nothing
    e1 = load_epochs(dir, 1); e5 = load_epochs(dir, 5); e60 = load_epochs(dir, 60)
    dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    isempty(validos) && return nothing
    estimable = ritmo.L5_estimable[k]
    ref_h = estimable ? ritmo.L5_medio_h[k] : ritmo.L5_reposo_medio_h[k]
    nw_as = estimable ? :missing : :rest

    # ---- A · sueño -------------------------------------------------------------------------
    ang = load_angles(dir)
    eje = longitudinal_axis(ang; mask = .!ang.nonwear)
    sueno = NamedTuple[]
    for (nombre, thr, nwh) in [("referencia (60°, 2 h)", 60, 2), ("umbral 45°", 45, 2), ("umbral 75°", 75, 2),
                               ("no uso ≥ 1 h", 60, 1), ("no uso ≥ 3 h", 60, 3)]
        n = sleep_nights(ang; axis = eje, nonwear = ang.nonwear, threshold = thr, max_nonwear = Hour(nwh))
        v = n[n.metodo .!= "sin ventana", :]
        push!(sueno, (id = id, variante = nombre, noches = nrow(n), postura = count(==("postura"), n.metodo),
            retirada = count(==("retirada"), n.metodo), sin_ventana = count(==("sin ventana"), n.metodo),
            punto_medio_h = isempty(v) ? NaN : media_circular_h(mod.(v.punto_medio_h, 24)),
            duracion_h = isempty(v) ? NaN : mean(v.duracion_h)))
    end
    sueno = DataFrame(sueno)
    sueno.dif_punto_medio_h = [dif_h(p, sueno.punto_medio_h[1]) for p in sueno.punto_medio_h]

    # ---- B · episodios ---------------------------------------------------------------------
    function con_peso(ep, dias = validos)
        e = copy(ep); e.peso = coverage_weights(e, clock_coverage(e1, dias; nonwear_as = nw_as)); e
    end
    variantes = [
        ("referencia (2 s, 3 s, 1 s)", e1, 2, 3), ("pausa 0 s", e1, 0, 3), ("pausa 5 s", e1, 5, 3),
        ("duración mín. 1 s", e1, 2, 1), ("duración mín. 5 s", e1, 2, 5), ("duración mín. 10 s", e1, 2, 10),
        ("epochs de 5 s", e5, 2, 5)]
    episodios = NamedTuple[]
    θref = NaN; zref = (NaN, NaN)
    for (nombre, e, tol, mins) in variantes
        ep = con_peso(vigorous_episodes(e, validos; threshold = UMBRAL_VIG, tolerance = tol, min_seconds = mins))
        tm = timing(ep, ref_h; peso = :peso)
        if startswith(nombre, "referencia"); θref = tm.θ_h; zref = (tm.R_cos, tm.R_sin); end
        push!(episodios, (id = id, variante = nombre, episodios = tm.n, θ_h = tm.θ_h, R = tm.R,
            Δθ_h = dif_h(tm.θ_h, θref), Δz = hypot(tm.R_cos - zref[1], tm.R_sin - zref[2])))
    end
    episodios = DataFrame(episodios)

    # ---- C · submuestreo de días -----------------------------------------------------------
    dias = sort(collect(validos)); nd = length(dias)
    ep_todos = vigorous_episodes(e1, validos; threshold = UMBRAL_VIG)
    function medidas(sub)
        s = Set(sub)
        ep = filter(r -> Date(r.medio) in s, ep_todos)
        ep.peso = coverage_weights(ep, clock_coverage(e1, s; nonwear_as = nw_as))
        tm = timing(ep, ref_h; peso = :peso)
        ii = is_iv(e60, s; nonwear_as = nw_as)
        (θ_h = tm.θ_h, R = tm.R, IS = ii.IS, IV = ii.IV, episodios = tm.n)
    end
    completo = medidas(dias)
    sub = NamedTuple[]
    for kd in 2:nd-1, c in combinaciones(dias, kd)
        m = medidas(c)
        push!(sub, (id = id, k = kd, θ_h = m.θ_h, Δθ_h = dif_h(m.θ_h, completo.θ_h), R = m.R, IS = m.IS, IV = m.IV, episodios = m.episodios))
    end
    sub = DataFrame(sub)
    @printf("%-12s %d días válidos | sueño: punto medio Δ %+.1f a %+.1f h | episodios: |Δθ| ≤ %.1f h (ref θ %.2f h, R %.2f) | días: control θ %.2f h\n",
            id, nd, minimum(sueno.dif_punto_medio_h), maximum(sueno.dif_punto_medio_h),
            maximum(abs, episodios.Δθ_h), θref, episodios.R[1], completo.θ_h)
    return sueno, episodios, sub
end

ritmo = CSV.read(RITMO, DataFrame)
ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "angulos_5s.arrow")), readdir(INTERIM))))
res = [(id, procesar(id, ritmo)) for id in ids]
res = [(id, r) for (id, r) in res if r !== nothing]
isempty(res) && error("Ningún registro procesable")
mkpath(OUT)
CSV.write(joinpath(OUT, "sens_sueno.csv"), reduce(vcat, [r[1] for (_, r) in res]))
CSV.write(joinpath(OUT, "sens_episodios.csv"), reduce(vcat, [r[2] for (_, r) in res]))
dias = reduce(vcat, [r[3] for (_, r) in res])
CSV.write(joinpath(OUT, "sens_dias.csv"), dias)
resumen = combine(groupby(dias, [:id, :k]), :id => length => :subconjuntos,
    :Δθ_h => (x -> maximum(abs, x)) => :max_abs_Δθ_h, :Δθ_h => (x -> median(abs.(x))) => :mediana_abs_Δθ_h,
    :R => minimum => :R_min, :R => maximum => :R_max, :IS => minimum => :IS_min, :IS => maximum => :IS_max,
    :IV => minimum => :IV_min, :IV => maximum => :IV_max)
CSV.write(joinpath(OUT, "sens_dias_resumen.csv"), resumen)
println("Salida en ", relpath(OUT, ROOT))
