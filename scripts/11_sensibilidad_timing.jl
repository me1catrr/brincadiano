# 11 · Sensibilidad del timing de la actividad vigorosa (hito H17).
#
# Qué hace
#   Repite el timing del script 06 (media circular ponderada por dosis, θ, R, bootstrap de
#   días) cambiando una sola cosa cada vez respecto a la referencia (ancla L5 o L5*, umbral
#   vigoroso de 464,6 mg, todos los días válidos):
#     · ancla: punto medio de la SPT (media de las noches), punto medio de la SPT de cada
#       noche, L5 de cada noche (cada episodio se refiere a la noche previa);
#     · umbral: ±20 % del vigoroso (371,7 y 557,5 mg) y el de AFMV (142,6 mg);
#     · días: solo laborables (lunes a viernes) o solo fin de semana.
#   Los nombres de las variantes son comunes a todos los registros: si L5 no es estimable,
#   la referencia es L5* (`ancla_ref`), y si el periodo de reposo nocturno es el del no uso, la SPT es
#   SPT* (`spt_retirada`); la tabla LaTeX lo marca con un asterisco en el ID.
#   Como en el 06, cada episodio pesa su dosis dividida por la cobertura del reloj (con los
#   días de la variante).
#   Para cada variante, Δθ frente a la referencia (circular, en (−12, 12] h) y la distancia
#   |Δz| en el plano (R cos θ, R sin θ), que es lo que ve el modelo metabólico.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/11_sensibilidad_timing.jl [ID ...]
#
# Depende de
#   05: epochs_1s.arrow y epochs_60s.arrow. 06: ritmo_resumen.csv (el ancla de referencia).
#   08: angulos_5s.arrow (periodos de reposo nocturno, recalculados aquí con las mismas reglas).
#
# Salida (individual; no se versiona)
#   outputs/timing/tables/sensibilidad_timing.{csv,tex}
#   outputs/timing/figures/sensibilidad_timing.{pdf,png}
#
# Consola: solo el ID y las cifras resumen. Detalle en el capítulo 12.

using Brincadiano.Metrics, Brincadiano.Circadian, Brincadiano.Sleep, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_T = joinpath(ROOT, "outputs", "timing")
const RITMO = joinpath(ROOT, "outputs", "circadian", "tables", "ritmo_resumen.csv")
const MIN_USO = 600
const UMBRAL_VIG = 464.6
const UMBRAL_AFMV = 142.6
const B = 2000

fmt(x; d = 1) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
fmt_signo(x; d = 1) = (isfinite(x) && round(x; digits = d) > 0 ? "+" : "") * fmt(x; d = d)
# Nombre de la variante en LaTeX: «%» escapado, coma decimal y espacio fino ante las unidades.
tex_nombre(v) = replace(v, r"(\d),(\d)" => s"\1{,}\2", " %" => "\\,\\%", " mg" => "\\,mg")
dif_h(a, b) = (d = mod(a - b, 24); d > 12 ? d - 24 : d)      # diferencia circular, (−12, 12]
media_circular_h(h) = phase_to_hours(circular_mean(2π .* h ./ 24).θ)
instante(d::Date, h::Real) = DateTime(d) + Millisecond(round(Int, h * 3_600_000))

function sensibilidad(id, ritmo)
    dir = joinpath(INTERIM, id)
    k = findfirst(==(id), ritmo.id)
    (k === nothing || !isfile(joinpath(dir, "angulos_5s.arrow"))) && return nothing
    e1 = load_epochs(dir, 1); e60 = load_epochs(dir, 60)
    dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    isempty(validos) && return nothing
    estimable = ritmo.L5_estimable[k]
    ref_h = estimable ? ritmo.L5_medio_h[k] : ritmo.L5_reposo_medio_h[k]
    nombre_ref = estimable ? "L5" : "L5*"

    # Periodos de reposo nocturno con las reglas del script 08.
    ang = load_angles(dir)
    noches = sleep_nights(ang; axis = longitudinal_axis(ang; mask = .!ang.nonwear), nonwear = ang.nonwear)
    noches = noches[noches.metodo .!= "sin ventana", :]
    spt_star = all(noches.metodo .== "retirada")
    l5n = [night_l5(e60, n.noche; nonwear_as = n.metodo == "retirada" ? :rest : :missing) for n in eachrow(noches)]
    spt_noche = [instante(n.noche, n.punto_medio_h) for n in eachrow(noches)]
    l5_noche = [instante(n.noche, h) for (n, h) in zip(eachrow(noches), l5n) if isfinite(h)]
    spt_media = media_circular_h(mod.(noches.punto_medio_h, 24))

    # Peso de cada episodio: dosis / cobertura del reloj, con la cobertura de los días que
    # entran en la variante (todos los válidos, o solo laborables o fin de semana).
    nw_as = estimable ? :missing : :rest
    function con_peso(ep, dias = validos)
        e = copy(ep); e.peso = coverage_weights(e, clock_coverage(e1, dias; nonwear_as = nw_as)); e
    end
    lab = Set(d for d in validos if dayofweek(d) <= 5); fds = Set(d for d in validos if dayofweek(d) >= 6)
    ep_vig = con_peso(vigorous_episodes(e1, validos; threshold = UMBRAL_VIG))
    variantes = [
        ("referencia", "Referencia (L5, 464,6 mg)", ep_vig, ref_h),
        ("ancla", "SPT, media de las noches", ep_vig, spt_media),
        ("ancla", "SPT de cada noche", ep_vig, spt_noche),
        ("ancla", "L5 de cada noche", ep_vig, l5_noche),
        ("umbral", "Umbral −20 % (371,7 mg)", con_peso(vigorous_episodes(e1, validos; threshold = 0.8UMBRAL_VIG)), ref_h),
        ("umbral", "Umbral +20 % (557,5 mg)", con_peso(vigorous_episodes(e1, validos; threshold = 1.2UMBRAL_VIG)), ref_h),
        ("umbral", "AFMV (142,6 mg)", con_peso(vigorous_episodes(e1, validos; threshold = UMBRAL_AFMV)), ref_h),
        ("dias", "Solo laborables", con_peso(filter(r -> dayofweek(r.medio) <= 5, ep_vig), lab), ref_h),
        ("dias", "Solo fin de semana", con_peso(filter(r -> dayofweek(r.medio) >= 6, ep_vig), fds), ref_h),
    ]
    filas = NamedTuple[]
    zref = (NaN, NaN); θref = NaN
    for (grupo, nombre, ep, ancla) in variantes
        tm = timing(ep, ancla; peso = :peso)
        bt = bootstrap_days(ep, ancla; B = B, peso = :peso)
        if grupo == "referencia"
            zref = (tm.R_cos, tm.R_sin); θref = tm.θ_h
        end
        ancla_h = ancla isa Real ? ancla : media_circular_h([hour(t) + minute(t) / 60 for t in ancla])
        push!(filas, (id = id, grupo = grupo, variante = nombre, ancla_ref = nombre_ref, spt_retirada = spt_star,
                      ancla_h = ancla_h,
                      dias = length(unique(Date.(ep.medio))), episodios = tm.n,
                      θ_h = tm.θ_h, R = tm.R, R_cos = tm.R_cos, R_sin = tm.R_sin,
                      Δθ_h = dif_h(tm.θ_h, θref), Δz = hypot(tm.R_cos - zref[1], tm.R_sin - zref[2]),
                      θ_ic_desde_h = bt.θ_ic_h[1], θ_ic_hasta_h = bt.θ_ic_h[2], θ_ic_amplitud_h = bt.amplitud_h,
                      R_ic_inf = bt.R_ic[1], R_ic_sup = bt.R_ic[2]))
    end
    f = DataFrame(filas)
    a = f[f.grupo .== "ancla", :]
    @printf("%-12s %d noches (%s); ancla SPT − %s = %+.1f h; Δθ por el ancla de %+.1f a %+.1f h, |Δz| ≤ %.2f\n",
            id, nrow(noches), spt_star ? "SPT*" : "postura", nombre_ref, dif_h(spt_media, ref_h),
            minimum(a.Δθ_h), maximum(a.Δθ_h), maximum(a.Δz))
    return f
end

function tex_sensibilidad(df)
    ids = unique(df.id)
    io = IOBuffer()
    println(io, "% Generado por scripts/11_sensibilidad_timing.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}X *{$(5length(ids))}{r}@{}}")
    println(io, "\\toprule")
    print(io, " ")
    for (j, id) in enumerate(ids)
        r = df[df.id .== id, :][1, :]
        print(io, " & \\multicolumn{5}{c}{\\texttt{", replace(id, "_" => "\\_"), "}",
              r.ancla_ref == "L5" && !r.spt_retirada ? "" : "\\textsuperscript{*}", "}")
    end
    println(io, " \\\\")
    for j in eachindex(ids)
        a = 2 + 5(j - 1); print(io, "\\cmidrule(l){$(a)-$(a + 4)}")
    end
    println(io)
    print(io, "Variante")
    for _ in ids
        print(io, " & Episodios & \$\\bar\\theta\$ (h) & \$\\Delta\\bar\\theta\$ (h) & \$R\$ & \$|\\Delta z|\$")
    end
    println(io, " \\\\")
    titulos = Dict("ancla" => "Ancla", "umbral" => "Umbral (ancla de referencia)", "dias" => "Días (ancla y umbral de referencia)")
    vars = unique(df[!, [:grupo, :variante]])
    grupo_prev = ""
    for v in eachrow(vars)
        if v.grupo != grupo_prev
            println(io, "\\midrule")
            haskey(titulos, v.grupo) && println(io, "\\multicolumn{$(1 + 5length(ids))}{@{}l}{\\emph{", titulos[v.grupo], "}} \\\\")
            grupo_prev = v.grupo
        end
        print(io, "\\quad ", tex_nombre(v.variante))
        for id in ids
            r = df[(df.id .== id) .& (df.variante .== v.variante), :]
            if nrow(r) == 0
                print(io, " & --- & --- & --- & --- & ---"); continue
            end
            r = r[1, :]
            ref = v.grupo == "referencia"
            print(io, " & ", r.episodios, " & ", fmt(r.θ_h), " & ", ref ? "---" : fmt_signo(r.Δθ_h), " & ",
                  fmt(r.R; d = 2), " & ", ref ? "---" : fmt(r.Δz; d = 2))
        end
        println(io, " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

isfile(RITMO) || error("falta $(RITMO) (ejecuta scripts/06_circadiano.jl)")
ritmo = CSV.read(RITMO, DataFrame)
ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "angulos_5s.arrow")), readdir(INTERIM))))
res = filter(!isnothing, map(id -> sensibilidad(id, ritmo), ids))
df = reduce(vcat, res)
for k in ("tables", "figures"); mkpath(joinpath(OUT_T, k)); end
CSV.write(joinpath(OUT_T, "tables", "sensibilidad_timing.csv"), df)
write(joinpath(OUT_T, "tables", "sensibilidad_timing.tex"), tex_sensibilidad(df))
guardar(joinpath(OUT_T, "figures", "sensibilidad_timing"), plot_timing_sensitivity(df);
        procedencia = (; tipo = :datos, registros = unique(df.id), funcion = "plot_timing_sensitivity", entradas = ["data/interim/<ID>/epochs_1s.arrow", "data/interim/<ID>/angulos_5s.arrow"]))
println("Sensibilidad en outputs/timing/tables/sensibilidad_timing.{csv,tex} y outputs/timing/figures/")
