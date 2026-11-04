# 09 · Luz: estados de calidad y métricas de exposición.
#
# Qué hace
#   Reconstruye la luz a epochs de 60 s (cada registro vale hasta el siguiente y es
#   desconocida dentro de un hueco) y le asigna un estado por minuto, en este orden
#   de prioridad: retirado, sueño, sin dato, oscuro e inactivo, oscuro y activo,
#   válido. El sueño sale de las noches del script 08 y el no uso, de los epochs de
#   60 s del 05. Resume por día y por registro los minutos por encima de 100 y de
#   1000, la hora media de luz alta y las ventanas ancladas a la SPT (2 h tras
#   despertar, 3 h antes de dormir). Las unidades son brutas: la conversión a lux
#   del wGT3X-BT no está resuelta, así que los umbrales no son lux.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/09_luz.jl [ID ...]    # sin ID, quien tenga angulos_5s.arrow
#
# Depende de
#   01 (lux.arrow y gaps.arrow), 05 (epochs_60s.arrow) y 08 (angulos_5s.arrow, de
#   donde se reconstruye el periodo de reposo nocturno). No lee las tablas de outputs/sleep/.
#
# Lo usan
#   10 (luz_60s.arrow, para el actograma) y 94 (la misma serie, exportada a LightLogR
#   para contrastar las métricas). `make sync-assets` copia tablas y figuras al
#   capítulo 13.
#
# Entrada
#   data/interim/<ID>/{lux,gaps,meta}.arrow
#   data/interim/<ID>/epochs_60s.arrow
#   data/interim/<ID>/angulos_5s.arrow
#
# Salida (individual; no se versiona; los días van numerados, sin fecha de calendario)
#   data/interim/<ID>/luz_60s.arrow              luz bruta, umbrales y estado de cada minuto
#   outputs/light/tables/luz_<ID>_dias.csv       una fila por día
#   outputs/light/tables/luz_resumen.{csv,tex}   una fila por registro; se reescribe entera
#   outputs/light/figures/<ID>_luz.{pdf,png}     luz por día, estados y periodo de reposo nocturno
#
# Un registro sin noches completas se omite. Consola: solo el ID y los porcentajes
# de cada estado. Umbrales fijos: 100 y 1000 en unidades brutas.

using Brincadiano.GT3X, Brincadiano.Metrics, Brincadiano.Sleep, Brincadiano.Light, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_T = joinpath(ROOT, "outputs", "light", "tables")
const OUT_F = joinpath(ROOT, "outputs", "light", "figures")
const THR = (100, 1000)
const GE = "\$\\ge\$"

fmt(x; d = 1) = (x isa Real && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
hhmm(h) = isnan(h) ? "---" : Dates.format(Time(0) + Minute(round(Int, 60h) % 1440), "HH:MM")

function luz(id)
    dir = joinpath(INTERIM, id)
    isfile(joinpath(dir, "angulos_5s.arrow")) || error("$id: faltan los ángulos (ejecuta scripts/08_sueno.jl)")
    r = read_interim(dir)
    e60 = load_epochs(dir, 60)
    epf = load_angles(dir)
    eje = longitudinal_axis(epf; mask = .!epf.nonwear)
    noches = sleep_nights(epf; axis = eje, nonwear = epf.nonwear, rec_start = first(r.acc.time), rec_stop = last(r.acc.time))
    if nrow(noches) == 0
        println(rpad(id, 12), " sin noches completas: se omite")
        return nothing
    end
    le = light_epochs(r.lux, r.gaps; origin = floor(first(r.acc.time), Day), stop = last(r.acc.time), thr = THR)
    st = light_states(le, e60, noches)
    st = st[(st.start .>= floor(first(r.acc.time), Minute)) .& (st.start .<= last(r.acc.time)), :]
    save_light(dir, st)
    dd = light_daily(st, noches; thr = THR)
    s = light_summary(dd)

    mkpath(OUT_T); mkpath(OUT_F)
    out = select(dd, Not(:dia)); out.dia = 1:nrow(out)            # sin fechas de calendario
    CSV.write(joinpath(OUT_T, "luz_$(id)_dias.csv"), select(out, :dia, Not(:dia)))
    pct(e) = 100 * count(==(e), st.estado) / nrow(st)
    fila = (id = id, (Symbol("pct_total_" * replace(e, " " => "_")) => pct(e) for e in ESTADOS_LUZ)..., s...)

    fig = plot_light_days(st; thr = THR)
    guardar(joinpath(OUT_F, "$(id)_luz"), fig;
            procedencia = (; tipo = :datos, registros = [id], funcion = "plot_light_days", entradas = ["data/interim/$id/"]))

    @printf("%-12s estados (%% del registro): válido %.0f, oscuro y activo %.0f, oscuro e inactivo %.0f, sin dato %.0f, sueño %.0f, retirado %.0f\n",
            id, (pct(e) for e in ESTADOS_LUZ)...)
    @printf("%-12s %d días; despierto: válido %.0f %%, oscuro y activo %.0f %%, sin dato %.0f %%; ≥%d: %.0f min/día, ≥%d: %.0f min/día (diario %.0f, finde %.0f)\n",
            id, s.dias, s.pct_valido, s.pct_oscuro_activo, s.pct_sin_dato, THR[1], s.min_sobre_bajo, THR[2], s.min_sobre_alto,
            s.min_sobre_alto_diario, s.min_sobre_alto_finde)
    @printf("%-12s segundo a segundo: %.0f min/día ≥%d, %.1f min/día ≥%d\n", id, s.seg_sobre_bajo, THR[1], s.seg_sobre_alto, THR[2])
    @printf("%-12s luz alta: hora media %s (DE %.1f h, %d días), %.1f entradas/día; AFMV %.0f min/día, %.0f con luz alta, %.0f a oscuras, %.0f sin luz\n",
            id, hhmm(s.hora_media_alta_h), s.hora_media_alta_de_h, s.dias_con_alta, s.entradas_alta, s.afmv_min, s.afmv_exterior,
            s.afmv_oscuro, s.afmv_sin_luz)
    @printf("%-12s 2 h tras despertar: %.0f %% válido, %.0f min ≥%d; 3 h antes de dormir: %.0f %% válido, %.0f min ≥%d\n",
            id, s.tras_despertar_pct_valido, s.tras_despertar_min_bajo, THR[1], s.antes_dormir_pct_valido, s.antes_dormir_min_bajo, THR[1])
    return fila
end

function tabla_tex(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/09_luz.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X *{$(nrow(df))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(replace.(df.id, "_" => "\\_"), " & "), " \\\\")
    println(io, "\\midrule")
    f(n, g) = println(io, n, " & ", join([g(r) for r in eachrow(df)], " & "), " \\\\")
    sec(t) = println(io, "\\multicolumn{$(nrow(df) + 1)}{@{}l}{\\textit{$t}} \\\\")
    sec("Calidad (\\% del tiempo despierto con el dispositivo puesto)")
    f("Días analizados", r -> string(r.dias))
    f("Válido / oscuro y activo / oscuro e inactivo / sin dato", r -> join(fmt.([r.pct_valido, r.pct_oscuro_activo, r.pct_oscuro_inactivo, r.pct_sin_dato]; d = 0), " / "))
    sec("Exposición (unidades brutas)")
    f("Minutos al día con media $GE $(THR[1]) / $(THR[2])", r -> fmt(r.min_sobre_bajo; d = 0) * " / " * fmt(r.min_sobre_alto; d = 0))
    f("Minutos acumulados segundo a segundo $GE $(THR[1]) / $(THR[2])", r -> fmt(r.seg_sobre_bajo; d = 0) * " / " * fmt(r.seg_sobre_alto; d = 1))
    f("Minutos $GE $(THR[2]): diario / fin de semana", r -> fmt(r.min_sobre_alto_diario; d = 0) * " / " * fmt(r.min_sobre_alto_finde; d = 0))
    f("Hora media con luz alta (DE entre días, h)", r -> hhmm(r.hora_media_alta_h) * " (" * fmt(r.hora_media_alta_de_h) * ")")
    f("Entradas en luz alta ($GE 5 min) al día", r -> fmt(r.entradas_alta))
    sec("Integración con la actividad (minutos de AFMV al día)")
    f("Total / con luz alta / a oscuras / sin luz", r -> join(fmt.([r.afmv_min, r.afmv_exterior, r.afmv_oscuro, r.afmv_sin_luz]; d = 0), " / "))
    sec("Ventanas ancladas al sueño")
    f("2\\,h tras despertar: \\% válido / minutos $GE $(THR[1])", r -> fmt(r.tras_despertar_pct_valido; d = 0) * " / " * fmt(r.tras_despertar_min_bajo; d = 0))
    f("3\\,h antes de dormir: \\% válido / minutos $GE $(THR[1])", r -> fmt(r.antes_dormir_pct_valido; d = 0) * " / " * fmt(r.antes_dormir_min_bajo; d = 0))
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

ids = filter(a -> !startswith(a, "--"), ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "angulos_5s.arrow")), readdir(INTERIM))))
res = filter(!isnothing, map(luz, ids))
df = DataFrame(res)
CSV.write(joinpath(OUT_T, "luz_resumen.csv"), df)
write(joinpath(OUT_T, "luz_resumen.tex"), tabla_tex(df))
println("Resumen en ", relpath(joinpath(OUT_T, "luz_resumen.{csv,tex}"), ROOT))
