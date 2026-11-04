# 15 · Figuras complementarias y cifras de sus pies.
#
# Qué hace
#   1. Cifras de los pies de las figuras que H25 añadió al informe (vista semanal e
#      histogramas de otros registros, ventana típica, ventanas de epochs): recalcula, con
#      la misma elección que los scripts 03 y 05, los agregados que se citan en los pies,
#      para que ninguna cifra del texto salga de una inspección a ojo.
#   2. Figuras de QC, ancla y timing (H26), de los registros con días válidos, a partir de
#      las tablas que ya escriben 04, 06 y 08; la cobertura del reloj y el barrido del
#      umbral se recalculan con las mismas funciones y parámetros que 06.
#   3. Figuras de luz × actividad × sueño (H27), de los registros con luz por epoch
#      (scripts/09): luz anclada al sueño, AFMV según la luz, densidad luz × ENMO, estados
#      por hora y perfil de 24 h de la luz válida.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/15_figuras_complementarias.jl [ID ...]   # sin ID, todos
#   julia --project=scripts scripts/15_figuras_complementarias.jl --solo-figuras   # sin las cifras de H25
#   julia --project=scripts scripts/15_figuras_complementarias.jl --solo-luz       # solo H27
#
# Depende de
#   01 (data/interim/<ID>/{acc,gaps,lux,battery}.arrow), 04 (qc_<ID>_{dias,tramos}.csv,
#   qc_resumen.csv), 05 (epochs_{1,5,60}s.arrow), 06 (timing_resumen.csv,
#   timing_por_dia_<ID>.csv), 08 (sueno_<ID>_noches.csv) y 09 (luz_60s.arrow,
#   luz_<ID>_dias.csv).
#
# Entrada
#   data/interim/<ID>/
#
# Salida (individual; no se versiona)
#   outputs/gt3x/tables/pies_registro.csv     por ID: muestras, % entre 23 y 08 h, luz máxima
#                                             y % de registros a 0, batería inicial y final,
#                                             tramos de idle sleep (número, mediana, máximo)
#   outputs/gt3x/tables/pies_ventana_tipica.csv  por ID: minuto típico de día (scripts/03)
#   outputs/metrics/tables/pies_epochs_ventana.csv  por ID: ventana de 5 min de scripts/05
#   outputs/qc/figures/dos_registros_dias_qc         composición de las horas de cada día
#   outputs/qc/figures/dos_registros_capsense        tramos inmóviles: duración × contacto
#   outputs/timing/figures/dos_registros_ancla_noches  punto medio del periodo − L5, por noche
#   outputs/timing/figures/dos_registros_timing_dias   vector de cada día y de la semana
#   outputs/timing/figures/dos_registros_cobertura     cobertura del reloj y dosis corregida
#   outputs/timing/figures/dos_registros_umbral        θ, R y episodios según el umbral
#   outputs/timing/tables/umbral_barrido.csv           el barrido (umbral, episodios, θ, R)
#   outputs/timing/tables/pies_h26.csv                 cifras de los pies de H26, por ID
#   outputs/light/figures/dos_registros_luz_anclada    luz alta y ventanas ancladas al sueño
#   outputs/light/figures/dos_registros_afmv_luz       minutos de AFMV según la luz, por día
#   outputs/light/figures/dos_registros_luz_enmo       densidad luz × ENMO por minuto
#   outputs/light/figures/dos_registros_estados_hora   estado de la luz según la hora
#   outputs/light/figures/dos_registros_perfil_luz     perfil de 24 h de la luz válida y de ENMO
#   outputs/light/tables/pies_h27.csv                  cifras de los pies de H27, por ID
#
# Consola: ID y las mismas cifras. Solo agregados; ningún dato de info.txt.

using Brincadiano.GT3X, Brincadiano.Metrics, Brincadiano.Circadian, Brincadiano.Light, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_GT3X = joinpath(ROOT, "outputs", "gt3x", "tables")
const OUT_METRICS = joinpath(ROOT, "outputs", "metrics", "tables")
const OUTPUTS = joinpath(ROOT, "outputs")
const UMBRAL_VIGOROSO = 464.6      # mg, como en Circadian y scripts/06 (Hildebrand, cadera)
const UMBRAL_AFMV = 142.6          # mg, como en scripts/06 y 11
const MIN_USO = 600                # min de uso de un día válido, como en scripts/06
# Barrido del umbral (mg): rejilla de 10 en 10 más los umbrales exactos del análisis (la
# referencia, ±20 % y el de AFMV), para que la referencia de Δθ sea la de 06.
const UMBRALES = sort(unique(vcat(collect(140.0:10.0:800.0),
                                  [UMBRAL_AFMV, 0.8UMBRAL_VIGOROSO, UMBRAL_VIGOROSO, 1.2UMBRAL_VIGOROSO])))

# ─── 1. Cifras de los pies (H25) ──────────────────────────────────────────────

"""Agregados de la vista semanal (lo que dibuja `plot_overview`)."""
function pies_registro(id, r)
    h = hour.(r.acc.time)
    idle = r.gaps[r.gaps.kind .== "idle_sleep", :]
    dur = Dates.value.(idle.stop .- idle.start) ./ 1000
    return (id = id, muestras = nrow(r.acc),
            pct_muestras_23_08 = 100count((h .>= 23) .| (h .< 8)) / nrow(r.acc),
            luz_max_bruta = isempty(r.lux.lux_raw) ? missing : maximum(r.lux.lux_raw),
            pct_registros_luz_0 = isempty(r.lux.lux_raw) ? missing : 100count(==(0), r.lux.lux_raw) / nrow(r.lux),
            bateria_inicial_v = isempty(r.battery.mv) ? missing : first(r.battery.mv) / 1000,
            bateria_final_v = isempty(r.battery.mv) ? missing : last(r.battery.mv) / 1000,
            idle_tramos = nrow(idle), idle_mediana_s = isempty(dur) ? missing : median(dur),
            idle_max_min = isempty(dur) ? missing : maximum(dur) / 60)
end

"""Minuto típico de día: la misma elección que `scripts/03_explorar_sujeto.jl`."""
function pies_ventana_tipica(id, r)
    fs = r.info.sample_rate
    nrm = vector_norm(r.acc)
    m = binned_summary(r.acc.time, abs.(nrm .- 1), Minute(1))
    dia = (m.n .== 60fs) .& (8 .<= hour.(m.start) .< 21)
    any(dia) || return nothing
    v = m.vmean[dia]; idx = findall(dia)
    i = idx[sortperm(v)[ceil(Int, 0.75 * length(v))]]
    sel = (r.acc.time .>= m.start[i]) .& (r.acc.time .< m.start[i] + Minute(1))
    return (id = id, hora = Dates.format(m.start[i], "HH:MM"), minutos_candidatos = length(v),
            media_desvio_mg = 1000m.vmean[i], norma_min_g = minimum(nrm[sel]), norma_max_g = maximum(nrm[sel]))
end

"""Ventana de 5 min de `scripts/05_metricas.jl` (alrededor del minuto con más ENMO)."""
function pies_epochs_ventana(id)
    dir = joinpath(INTERIM, id)
    e1 = load_epochs(dir, 1); e5 = load_epochs(dir, 5); e60 = load_epochs(dir, 60)
    full = (e60.cov .== 1) .& .!e60.nonwear
    any(full) || return nothing
    k = findall(full)[argmax(e60.enmo[full])]
    t0 = e60.start[k] - Minute(2); t1 = t0 + Minute(5)
    w1 = (e1.start .>= t0) .& (e1.start .< t1) .& (e1.n .> 0)
    w5 = (e5.start .>= t0) .& (e5.start .< t1) .& (e5.n .> 0)
    return (id = id, hora_minuto = Dates.format(e60.start[k], "HH:MM"), enmo1_max_mg = maximum(e1.enmo[w1]),
            enmo5_max_mg = maximum(e5.enmo[w5]), enmo60_minuto_mg = e60.enmo[k],
            segundos_vigorosos = count(e1.enmo[w1] .>= UMBRAL_VIGOROSO))
end

function cifras_pies(ids)
    reg = DataFrame(); tip = DataFrame(); epw = DataFrame()
    for id in ids
        r = read_interim(joinpath(INTERIM, id))
        a = pies_registro(id, r); push!(reg, a; promote = true)
        @printf("%s: %d muestras, %.1f %% entre 23 y 08 h; luz máx. %s (%s %% a 0); batería %s → %s V; idle sleep %d tramos (mediana %s s, máx. %s min)\n",
                id, a.muestras, a.pct_muestras_23_08, a.luz_max_bruta, a.pct_registros_luz_0 === missing ? "—" : @sprintf("%.0f", a.pct_registros_luz_0),
                a.bateria_inicial_v, a.bateria_final_v, a.idle_tramos, a.idle_mediana_s, a.idle_max_min === missing ? "—" : @sprintf("%.0f", a.idle_max_min))
        t = pies_ventana_tipica(id, r)
        t === nothing || (push!(tip, t; promote = true);
            @printf("  minuto típico %s: |‖a‖−1| medio %.0f mg; norma %.2f–%.2f g (%d candidatos)\n", t.hora, t.media_desvio_mg, t.norma_min_g, t.norma_max_g, t.minutos_candidatos))
        if isfile(joinpath(INTERIM, id, "epochs_60s.arrow"))
            e = pies_epochs_ventana(id)
            e === nothing || (push!(epw, e; promote = true);
                @printf("  ventana de epochs (%s): máx. ENMO 1 s %.0f mg, 5 s %.0f mg; minuto %.0f mg; %d s vigorosos\n",
                        e.hora_minuto, e.enmo1_max_mg, e.enmo5_max_mg, e.enmo60_minuto_mg, e.segundos_vigorosos))
        end
    end
    mkpath(OUT_GT3X); mkpath(OUT_METRICS)
    CSV.write(joinpath(OUT_GT3X, "pies_registro.csv"), reg)
    CSV.write(joinpath(OUT_GT3X, "pies_ventana_tipica.csv"), tip)
    CSV.write(joinpath(OUT_METRICS, "pies_epochs_ventana.csv"), epw)
end

# ─── 2. Figuras de QC, ancla y timing (H26) ───────────────────────────────────

leer(sub, archivo) = CSV.read(joinpath(OUTPUTS, sub, "tables", archivo), DataFrame)
procedencia(ids, funcion, entradas) = (; tipo = :datos, registros = ids, funcion, entradas)

"Ancla, días válidos y cobertura del reloj exactamente como en scripts/06_circadiano.jl."
function contexto_timing(id)
    dir = joinpath(INTERIM, id)
    e1 = load_epochs(dir, 1); e60 = load_epochs(dir, 60)
    dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    rm_ = rar_metrics(profile24h(e60, validos; nonwear_as = :missing))
    ancla = rm_.estimable ? rm_ : rar_metrics(profile24h(e60, validos; nonwear_as = :rest))
    cob = clock_coverage(e1, validos; nonwear_as = rm_.estimable ? :missing : :rest)
    return (; e1, validos, ancla_h = ancla.L5_medio_h, cob, estimable = rm_.estimable)
end

function figuras_h26(ids)
    res = leer("timing", "timing_resumen.csv")
    ids = [id for id in ids if id in res.id]          # solo registros con días válidos (06)
    isempty(ids) && return println("H26: ningún registro con días válidos")
    qcr = leer("qc", "qc_resumen.csv")
    pies = DataFrame(); barrido = DataFrame()
    set_theme!(tema_informe())

    # QC: composición de los días y tramos inmóviles.
    dias = [(id, leer("qc", "qc_$(id)_dias.csv")) for id in ids]
    guardar(joinpath(OUTPUTS, "qc", "figures", "dos_registros_dias_qc"), plot_qc_days(dias);
            procedencia = procedencia(ids, "plot_qc_days", ["outputs/qc/tables/qc_<ID>_dias.csv"]))
    tramos = [(id, leer("qc", "qc_$(id)_tramos.csv"), qcr.contacto_ref_pct[findfirst(==(id), qcr.id)] / 100) for id in ids]
    guardar(joinpath(OUTPUTS, "qc", "figures", "dos_registros_capsense"), plot_capsense_runs(tramos);
            procedencia = procedencia(ids, "plot_capsense_runs", ["outputs/qc/tables/qc_<ID>_tramos.csv", "outputs/qc/tables/qc_resumen.csv"]))

    # Ancla noche a noche.
    noches = [(id, leer("sleep", "sueno_$(id)_noches.csv")) for id in ids]
    guardar(joinpath(OUTPUTS, "timing", "figures", "dos_registros_ancla_noches"), plot_anchor_nights(noches);
            procedencia = procedencia(ids, "plot_anchor_nights", ["outputs/sleep/tables/sueno_<ID>_noches.csv"]))

    # Timing día a día frente al de la semana.
    tdias = map(ids) do id
        r = res[findfirst(==(id), res.id), :]
        (id, leer("timing", "timing_por_dia_$(id).csv"), r.θ_h, r.R, (r.θ_ic_desde_h, r.θ_ic_hasta_h))
    end
    guardar(joinpath(OUTPUTS, "timing", "figures", "dos_registros_timing_dias"), plot_timing_days(tdias);
            procedencia = procedencia(ids, "plot_timing_days", ["outputs/timing/tables/timing_por_dia_<ID>.csv", "outputs/timing/tables/timing_resumen.csv"]))

    # Cobertura del reloj y barrido del umbral (mismas funciones y parámetros que 06).
    cobs = []
    for id in ids
        c = contexto_timing(id)
        ep = vigorous_episodes(c.e1, c.validos; threshold = UMBRAL_VIGOROSO)
        ep.peso = coverage_weights(ep, c.cob)
        push!(cobs, (id, c.cob, ep, length(c.validos), c.ancla_h))
        for u in UMBRALES
            e = vigorous_episodes(c.e1, c.validos; threshold = u)
            if nrow(e) == 0
                push!(barrido, (id = id, umbral_mg = u, episodios = 0, θ_h = NaN, R = NaN)); continue
            end
            e.peso = coverage_weights(e, c.cob)
            tm = timing(e, c.ancla_h; peso = :peso)
            push!(barrido, (id = id, umbral_mg = u, episodios = tm.n, θ_h = tm.θ_h, R = tm.R))
        end
        # Control: con el umbral de referencia, el θ debe ser el de timing_resumen.csv.
        tm = timing(ep, c.ancla_h; peso = :peso)
        r = res[findfirst(==(id), res.id), :]
        abs(tm.θ_h - r.θ_h) < 1e-6 || @warn "$(id): θ ($(tm.θ_h)) distinto del de 06 ($(r.θ_h))"
        b = barrido[barrido.id .== id, :]
        δ(θ) = abs(mod(θ - tm.θ_h + 12, 24) - 12)
        en(lo, hi) = (lo .<= b.umbral_mg .<= hi) .& isfinite.(b.θ_h)
        n = noches[findfirst(x -> x[1] == id, noches)][2]
        d = tdias[findfirst(x -> x[1] == id, tdias)][2]
        okb = isfinite.(b.θ_h); kmax = findall(okb)[argmax(δ.(b.θ_h[okb]))]
        push!(pies, (id = id, θ_h = tm.θ_h, θ_coincide_06 = abs(tm.θ_h - r.θ_h) < 1e-6,
                     umbral_max_δθ_mg = b.umbral_mg[kmax], episodios_max_δθ = b.episodios[kmax],
                     max_δθ_hasta_372_h = maximum(δ.(b.θ_h[en(0, 372)])),
                     max_δθ_370_560_h = maximum(δ.(b.θ_h[en(370, 560)])),
                     max_δθ_140_800_h = maximum(δ.(b.θ_h[en(140, 800)])),
                     R_min_140_800 = minimum(b.R[en(140, 800)]), R_max_140_800 = maximum(b.R[en(140, 800)]),
                     episodios_140 = b.episodios[1], episodios_800 = b.episodios[end],
                     cobertura_min = minimum(c.cob), cobertura_max = maximum(c.cob), dias_validos = length(c.validos),
                     dias_con_episodios = nrow(d), dias_R_bajo_0_3 = count(<(0.3), d.R),
                     spt_L5_min_h = minimum(filter(isfinite, n.spt_menos_L5_h)),
                     spt_L5_max_h = maximum(filter(isfinite, n.spt_menos_L5_h)),
                     pct_spt_idle_mediana = (v = filter(isfinite, n.pct_spt_idle); isempty(v) ? NaN : median(v))))
    end
    guardar(joinpath(OUTPUTS, "timing", "figures", "dos_registros_cobertura"), plot_clock_coverage(cobs);
            procedencia = procedencia(ids, "plot_clock_coverage", ["data/interim/<ID>/epochs_1s.arrow", "data/interim/<ID>/epochs_60s.arrow"]))
    guardar(joinpath(OUTPUTS, "timing", "figures", "dos_registros_umbral"),
            plot_threshold_sweep([(id, barrido[barrido.id .== id, :]) for id in ids];
                                 referencia = UMBRAL_VIGOROSO, marcas = (UMBRAL_AFMV, 0.8UMBRAL_VIGOROSO, 1.2UMBRAL_VIGOROSO));
            procedencia = procedencia(ids, "plot_threshold_sweep", ["data/interim/<ID>/epochs_1s.arrow", "data/interim/<ID>/epochs_60s.arrow"]))
    CSV.write(joinpath(OUTPUTS, "timing", "tables", "umbral_barrido.csv"), barrido)
    CSV.write(joinpath(OUTPUTS, "timing", "tables", "pies_h26.csv"), pies)
    for r in eachrow(pies)
        @printf("%s: θ %.2f h (= 06: %s); |Δθ| máx. %.2f h con 370–560 mg y %.2f h con 140–800 mg; R %.2f–%.2f; episodios %d (140 mg) → %d (800 mg)\n",
                r.id, r.θ_h, r.θ_coincide_06, r.max_δθ_370_560_h, r.max_δθ_140_800_h, r.R_min_140_800, r.R_max_140_800,
                r.episodios_140, r.episodios_800)
        @printf("  cobertura del reloj %d–%d días (de %d válidos); %d días con episodios, %d con R < 0,3; SPT − L5 de %.2f a %.2f h; %% de SPT en idle (mediana) %.1f\n",
                r.cobertura_min, r.cobertura_max, r.dias_validos, r.dias_con_episodios, r.dias_R_bajo_0_3,
                r.spt_L5_min_h, r.spt_L5_max_h, r.pct_spt_idle_mediana)
    end
end

# ─── 3. Figuras de luz × actividad × sueño (H27) ──────────────────────────────

function figuras_h27(ids)
    ids = [id for id in ids if isfile(joinpath(INTERIM, id, "luz_60s.arrow")) &&
                               isfile(joinpath(OUTPUTS, "light", "tables", "luz_$(id)_dias.csv"))]
    isempty(ids) && return println("H27: ningún registro con luz por epoch")
    set_theme!(tema_informe())
    serie = [(id, load_light(joinpath(INTERIM, id))) for id in ids]
    dias = [(id, leer("light", "luz_$(id)_dias.csv")) for id in ids]
    fig(nombre, f, funcion, entradas) = guardar(joinpath(OUTPUTS, "light", "figures", nombre), f;
                                                  procedencia = procedencia(ids, funcion, entradas))
    e_dias = ["outputs/light/tables/luz_<ID>_dias.csv"]; e_serie = ["data/interim/<ID>/luz_60s.arrow"]
    fig("dos_registros_luz_anclada", plot_light_anchored(dias), "plot_light_anchored", e_dias)
    fig("dos_registros_afmv_luz", plot_mvpa_light(dias), "plot_mvpa_light", e_dias)
    fig("dos_registros_luz_enmo", plot_light_enmo_density(serie), "plot_light_enmo_density", e_serie)
    fig("dos_registros_estados_hora", plot_light_states_by_hour(serie), "plot_light_states_by_hour", e_serie)
    fig("dos_registros_perfil_luz", plot_light_profile(serie), "plot_light_profile", e_serie)

    pies = DataFrame()
    for ((id, d), (_, dd)) in zip(serie, dias)
        h = hour.(d.start)
        pv = [begin n = count(h .== k); n == 0 ? NaN : 100count((h .== k) .& (d.estado .== "válido")) / n end for k in 0:23]
        kmax = argmax(replace(pv, NaN => -1.0)) - 1
        dp = (d.estado .!= "sueño") .& (d.estado .!= "retirado")
        conocida = dp .& (d.lux_obs .>= 0.5) .& isfinite.(d.lux_media)
        push!(pies, (id = id, dias = nrow(dd), dias_con_luz_alta = count(isfinite, dd.primera_alta_h),
                     hora_max_valido = kmax, pct_valido_hora_max = pv[kmax + 1],
                     horas_con_10_validos = count(k -> count((h .== k) .& (d.estado .== "válido")) >= 10, 0:23),
                     min_despierto_luz_conocida = count(conocida),
                     pct_luz_0_con_luz_conocida = 100count(conocida .& (d.lux_media .< 1)) / max(count(conocida), 1),
                     afmv_dia_medio = sum(dd.afmv_min) / nrow(dd),
                     afmv_oscuro_dia_medio = sum(dd.afmv_oscuro) / nrow(dd)))
    end
    CSV.write(joinpath(OUTPUTS, "light", "tables", "pies_h27.csv"), pies)
    for r in eachrow(pies)
        @printf("%s: %d de %d días con luz alta; hora con más luz válida %02d h (%.0f %%); %d horas con ≥ 10 min válidos; luz conocida a 0 el %.0f %% de %d min despierto\n",
                r.id, r.dias_con_luz_alta, r.dias, r.hora_max_valido, r.pct_valido_hora_max, r.horas_con_10_validos,
                r.pct_luz_0_con_luz_conocida, r.min_despierto_luz_conocida)
    end
end

ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "meta.arrow")), readdir(INTERIM))))
if "--solo-luz" in ARGS
    figuras_h27(ids)
else
    "--solo-figuras" in ARGS || cifras_pies(ids)
    figuras_h26(ids)
    figuras_h27(ids)
end
