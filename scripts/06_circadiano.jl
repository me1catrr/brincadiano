# 06 · Ritmo de reposo-actividad y timing circular de la actividad vigorosa.
#
# Qué hace
#   Con los epochs ya calibrados y marcados del script 05 calcula dos cosas. El ritmo:
#   perfil de 24 h, L5, M10, amplitud relativa e IS/IV, en días con al menos 600 min de
#   uso. El L5 se estima con el no uso como hora ausente; si ese perfil cubre menos del
#   90 % del día, el ancla pasa a un L5* que trata el no uso como reposo (supuesto que
#   la tabla declara). El timing: episodios con ENMO de 1 s ≥ 464,6 mg (Hildebrand
#   2014, cadera, 7–11 años), pausas de hasta 2 s y al menos 3 s vigorosos, y su media
#   circular ponderada por dosis (θ respecto al ancla, R, R·cos θ, R·sin θ, Rayleigh y
#   bootstrap de días). Desde H20, cada episodio pesa su dosis dividida por la cobertura
#   del reloj en su instante medio (días en que se observó esa hora; Circadian.
#   coverage_weights): es la media circular del perfil medio de dosis. El resumen guarda
#   también θ y R sin corregir. Una pasada de sensibilidad repite los episodios a 142,6 mg
#   (AFMV) y entra en el resumen, no en el CSV de episodios.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/06_circadiano.jl [ID ...]    # sin ID, quien tenga epochs_60s.arrow
#
# Depende de
#   05: epochs_1s.arrow (episodios) y epochs_60s.arrow (perfil, IS/IV, días válidos).
#   El no uso ya viene marcado; no vuelve a leer qc_tramos.arrow.
#
# Lo usan
#   08, si existe outputs/circadian/tables/ritmo_resumen.csv, para restar el ancla del
#   punto medio del periodo de reposo nocturno (si no existe, el sueño se calcula igual).
#   10 lee ritmo_resumen, timing_resumen y episodios_<ID>.csv.
#   92 lee los tres para L5/M10, IS/IV y la media circular frente a numpy.
#   El 95 repite este cálculo sobre la señal imputada y no lee estas tablas.
#   Capítulos 10 y 11 del informe, vía `make sync-assets`.
#
# Entrada
#   data/interim/<ID>/epochs_1s.arrow
#   data/interim/<ID>/epochs_60s.arrow
#
# Salida (individual; no se versiona)
#   outputs/circadian/tables/ritmo_resumen.{csv,tex}     L5, L5*, M10, RA, IS, IV; una fila por registro
#   outputs/circadian/figures/<ID>_perfil24h.{pdf,png}   las dos variantes del perfil
#   outputs/timing/tables/timing_resumen.{csv,tex}       θ, R, R·cos θ, R·sin θ, Rayleigh, bootstrap
#   outputs/timing/tables/episodios_<ID>.csv             cada episodio vigoroso (umbral 464,6 mg),
#                                                        con `dosis` y `peso` (dosis / cobertura)
#   outputs/timing/tables/timing_por_dia_<ID>.csv        θ y R de cada día válido
#   outputs/timing/figures/<ID>_diana.{pdf,png}          diana; se omite si no hay episodios
#
# Un registro sin días válidos se omite y no entra en los resúmenes.
# Consola: solo el ID, el ancla elegido y los recuentos. Detalle en los capítulos 10 y 11.

using Brincadiano.Metrics, Brincadiano.Circadian, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_C = joinpath(ROOT, "outputs", "circadian")
const OUT_T = joinpath(ROOT, "outputs", "timing")
const MIN_USO = 600
const UMBRAL_VIG = 464.6
const UMBRAL_AFMV = 142.6

# Coma decimal (convención del informe).
fmt(x; d = 2) = (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
hhmm(h) = isnan(h) ? "---" : Dates.format(Time(0) + Minute(round(Int, 60h) % 1440), "HH:MM")

function circadiano(id)
    dir = joinpath(INTERIM, id)
    isfile(joinpath(dir, "epochs_60s.arrow")) || error("$id: faltan los epochs (ejecuta scripts/05_metricas.jl)")
    e1 = load_epochs(dir, 1); e60 = load_epochs(dir, 60)
    dd = daily_summary(e60)
    validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    if isempty(validos)
        println(rpad(id, 12), " sin días válidos: se omite")
        return nothing
    end
    pm = profile24h(e60, validos; nonwear_as = :missing)
    pr = profile24h(e60, validos; nonwear_as = :rest)
    rm_, rr = rar_metrics(pm), rar_metrics(pr)
    ancla, nombre_ancla = rm_.estimable ? (rm_, "L5") : (rr, "L5 (no uso = reposo)")
    ivm = is_iv(e60, validos; nonwear_as = :missing); ivr = is_iv(e60, validos; nonwear_as = :rest)

    # Cobertura del reloj con la misma regla del no uso que el ancla: cada episodio pesa su
    # dosis dividida por los días en que se observó su hora (media circular del perfil medio).
    cob = clock_coverage(e1, validos; nonwear_as = rm_.estimable ? :missing : :rest)
    ep = vigorous_episodes(e1, validos; threshold = UMBRAL_VIG)
    ep.peso = coverage_weights(ep, cob)
    tm = timing(ep, ancla.L5_medio_h; peso = :peso)
    tm0 = timing(ep, ancla.L5_medio_h)                                   # sin corregir (sensibilidad)
    bt = bootstrap_days(ep, ancla.L5_medio_h; B = 2000, peso = :peso)
    CSV.write(joinpath(OUT_T, "tables", "timing_por_dia_$(id).csv"), timing_by_day(ep, ancla.L5_medio_h))
    epa = vigorous_episodes(e1, validos; threshold = UMBRAL_AFMV)
    epa.peso = coverage_weights(epa, cob)
    tma = timing(epa, ancla.L5_medio_h; peso = :peso)

    for d in (OUT_C, OUT_T), k in ("tables", "figures"); mkpath(joinpath(d, k)); end
    ep.fase_h = phase_to_hours.(tm.fase)
    CSV.write(joinpath(OUT_T, "tables", "episodios_$(id).csv"), ep)

    fig = plot_profile24h([("no uso como ausente", pm), ("no uso como reposo", pr)], ancla)
    guardar(joinpath(OUT_C, "figures", "$(id)_perfil24h"), fig;
            procedencia = (; tipo = :datos, registros = [id], funcion = "plot_profile24h", entradas = ["data/interim/$id/epochs_60s.arrow"]))
    orden = Dict(d => i for (i, d) in enumerate(sort(collect(validos))))
    if nrow(ep) > 0
        fig = plot_timing_circular(tm.fase, ep.peso, [orden[Date(t)] for t in ep.medio], tm.θ, tm.R;
                                   ancla = nombre_ancla == "L5" ? "punto medio de L5" : "punto medio de L5*")
        guardar(joinpath(OUT_T, "figures", "$(id)_diana"), fig;
                procedencia = (; tipo = :datos, registros = [id], funcion = "plot_timing_circular", entradas = ["data/interim/$id/epochs_1s.arrow", "data/interim/$id/epochs_60s.arrow"]))
    end

    ritmo = (id = id, dias_validos = length(validos), cobertura_perfil = rm_.cobertura_perfil,
             L5_estimable = rm_.estimable, L5 = rm_.L5, L5_medio_h = rm_.L5_medio_h,
             L5_reposo = rr.L5, L5_reposo_medio_h = rr.L5_medio_h,
             M10 = ancla.M10, M10_medio_h = ancla.M10_medio_h, RA = ancla.RA,
             IS = ivm.IS, IV = ivm.IV, IS_reposo = ivr.IS, IV_reposo = ivr.IV, ancla = nombre_ancla)
    tim = (id = id, ancla = nombre_ancla, ancla_h = ancla.L5_medio_h, episodios = tm.n,
           seg_vig_dia = sum(ep.segundos_vigorosos; init = 0) / max(length(validos), 1),
           θ_h = tm.θ_h, R = tm.R, R_cos = tm.R_cos, R_sin = tm.R_sin,
           R_sin_ponderar = tm.R_sin_ponderar, rayleigh_p = tm.rayleigh_p, R2 = tm.R2, θ2_h = tm.θ2_h,
           θ_ic_desde_h = bt.θ_ic_h[1], θ_ic_hasta_h = bt.θ_ic_h[2], θ_ic_amplitud_h = bt.amplitud_h,
           R_ic_inf = bt.R_ic[1], R_ic_sup = bt.R_ic[2],
           θ_informable = tm.n >= 10 && bt.amplitud_h < 12,     # cap. 11, casos degenerados
           episodios_afmv = tma.n, θ_h_afmv = tma.θ_h, R_afmv = tma.R,
           cobertura_min = minimum(cob), cobertura_max = maximum(cob),
           θ_h_sin_corregir = tm0.θ_h, R_sin_corregir = tm0.R)
    @printf("%-12s %d días válidos; L5 %s (cobertura %.0f %%), ancla %s a las %s; RA %.2f IS %.2f IV %.2f\n",
            id, length(validos), rm_.estimable ? "estimable" : "NO estimable", 100rm_.cobertura_perfil,
            nombre_ancla, hhmm(ancla.L5_medio_h), ancla.RA,
            rm_.estimable ? ivm.IS : ivr.IS, rm_.estimable ? ivm.IV : ivr.IV)
    @printf("%-12s %d episodios vigorosos (%.0f s/día); θ = %.1f h tras el ancla, R = %.2f, R₂ = %.2f (eje %.1f/%.1f h); Rayleigh p = %.2g; AFMV: θ = %.1f h, R = %.2f\n",
            id, tm.n, tim.seg_vig_dia, tm.θ_h, tm.R, tm.R2, tm.θ2_h, mod(tm.θ2_h + 12, 24), tm.rayleigh_p, tma.θ_h, tma.R)
    @printf("%-12s IC 95 %% (bootstrap de días): θ de %.1f a %.1f h (%.1f h de amplitud); R de %.2f a %.2f\n",
            id, bt.θ_ic_h[1], bt.θ_ic_h[2], bt.amplitud_h, bt.R_ic[1], bt.R_ic[2])
    return ritmo, tim
end

function tex_ritmo(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/06_circadiano.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}l *{10}{>{\\raggedleft\\arraybackslash}X}@{}}")
    println(io, "\\toprule")
    println(io, "ID & Días válidos & Cobertura del perfil (\\%) & L5 (mg) & Punto medio de L5 & L5* (mg) & ",
            "Punto medio de L5* & M10 (mg) & RA & IS / IV \\\\")
    println(io, "\\midrule")
    for r in eachrow(df)
        is_, iv_ = r.L5_estimable ? (r.IS, r.IV) : (r.IS_reposo, r.IV_reposo)
        println(io, replace(r.id, "_" => "\\_"), " & ", r.dias_validos, " & ", fmt(100r.cobertura_perfil; d = 0), " & ",
                r.L5_estimable ? fmt(r.L5; d = 1) : "---", " & ", r.L5_estimable ? hhmm(r.L5_medio_h) : "no estimable", " & ",
                fmt(r.L5_reposo; d = 1), " & ", hhmm(r.L5_reposo_medio_h), " & ", fmt(r.M10; d = 1), " & ",
                fmt(r.RA), " & ", fmt(is_), " / ", fmt(iv_), r.L5_estimable ? "" : "\\textsuperscript{*}", " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

function tex_timing(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/06_circadiano.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}l *{10}{>{\\raggedleft\\arraybackslash}X}@{}}")
    println(io, "\\toprule")
    println(io, "ID & Ancla & Episodios & Segundos vigorosos por día & \$\\bar\\theta\$ (h tras el ancla) & \$R\$ & \$R_2\$ & ",
            "IC 95\\,\\% de \$\\bar\\theta\$ (h) & \$R\\cos\\bar\\theta\$ / \$R\\sin\\bar\\theta\$ & Rayleigh \$p\$ & AFMV: \$\\bar\\theta\$ (h) / \$R\$ \\\\")
    println(io, "\\midrule")
    for r in eachrow(df)
        println(io, replace(r.id, "_" => "\\_"), " & ", hhmm(r.ancla_h), r.ancla == "L5" ? "" : "\\textsuperscript{*}", " & ",
                r.episodios, " & ", fmt(r.seg_vig_dia; d = 0), " & ", fmt(r.θ_h; d = 1), " & ", fmt(r.R), " & ", fmt(r.R2), " & ",
                fmt(r.θ_ic_desde_h; d = 1), "--", fmt(r.θ_ic_hasta_h; d = 1), " & ",
                fmt(r.R_cos), " / ", fmt(r.R_sin), " & ", r.rayleigh_p < 0.001 ? "\$<0{,}001\$" : fmt(r.rayleigh_p; d = 3), " & ",
                fmt(r.θ_h_afmv; d = 1), " / ", fmt(r.R_afmv), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "epochs_60s.arrow")), readdir(INTERIM))))
res = filter(!isnothing, map(circadiano, ids))
ritmo = DataFrame(first.(res)); tim = DataFrame(last.(res))
CSV.write(joinpath(OUT_C, "tables", "ritmo_resumen.csv"), ritmo)
write(joinpath(OUT_C, "tables", "ritmo_resumen.tex"), tex_ritmo(ritmo))
CSV.write(joinpath(OUT_T, "tables", "timing_resumen.csv"), tim)
write(joinpath(OUT_T, "tables", "timing_resumen.tex"), tex_timing(tim))
println("Resúmenes en outputs/circadian/tables y outputs/timing/tables")
