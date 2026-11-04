# 94 · Métricas de luz de Brincadiano frente a LightLogR.
#
# Qué hace
#   Contrasta las métricas, no la lectura. LightLogR no abre ActiGraph: recibe la serie
#   por minuto que ya construyó el script 09. Por día compara minutos ≥ 100 y ≥ 1000
#   (unidades brutas) y la primera, la última y la hora media con luz ≥ 100. En el
#   registro completo compara IS e IV de la luz en vigilia (Circadian.is_iv frente a
#   interdaily_stability / intradaily_variability). Los minutos que no son «válido»
#   van como hueco en la serie de luz válida; sueño, retirada y sin dato van como
#   hueco en la serie de vigilia.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/94_validar_luz_lightlogr.jl --confirmo
#   Sin --confirmo termina tras un aviso. Hace falta autorización: lanza R.
#
# Depende de
#   09  data/interim/<ID>/luz_60s.arrow (también filtra qué registros entran)
#   08  angulos_5s.arrow, para reconstruir las noches que pide light_daily
#   01  el registro, por el inicio y el fin
#   scripts/R/lightlogr_metricas.R y el renv de scripts/R/ (LightLogR 0.10.6)
#
# Lo usan
#   El capítulo 14: `make sync-assets` copia outputs/validation/tables/ a
#   report/assets/13_validacion/brincadiano/tables/. El CSV por registro
#   (outputs/validation/luz/) no entra en esa copia. Ningún script de Julia, tampoco
#   el 10, lee la tabla de diferencias.
#
# Entrada
#   data/interim/<ID>/luz_60s.arrow y angulos_5s.arrow
#
# Salida (individual; no se versiona)
#   outputs/validation/luz/<ID>.csv
#       hora local ISO, luz_valida y luz_vigilia (NA fuera de su definición). Es la
#       entrada de LightLogR y se conserva.
#   outputs/validation/tables/luz_lightlogr.csv     una fila por registro: diferencias máximas
#   outputs/validation/tables/luz_lightlogr.tex
#   Los CSV por día y global que devuelve R son temporales y se borran al terminar.
#
# Consola: solo el ID y las diferencias. La serie exportada no lleva nombre ni fecha
# de nacimiento; la hora local sí, porque LightLogR la necesita para el día.

using Brincadiano.GT3X, Brincadiano.Sleep, Brincadiano.Light, Brincadiano.Circadian
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const VAL = joinpath(ROOT, "outputs", "validation")
const OUT = joinpath(VAL, "tables")
const RENV = joinpath(ROOT, "scripts", "R")
const VIGILIA = ("válido", "oscuro y activo", "oscuro e inactivo")

fmt(x; d = 2) = (x isa Real && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
mean_(v) = isempty(v) ? NaN : sum(v) / length(v)
maxabs(v) = isempty(v) ? NaN : maximum(abs.(v))
num(v) = [ismissing(x) || x == "NA" ? NaN : Float64(x) for x in v]

function comparar(id)
    dir = joinpath(INTERIM, id)
    r = read_interim(dir)
    st = load_light(dir)
    epf = load_angles(dir)
    eje = longitudinal_axis(epf; mask = .!epf.nonwear)
    noches = sleep_nights(epf; axis = eje, nonwear = epf.nonwear, rec_start = first(r.acc.time), rec_stop = last(r.acc.time))
    # 1. Exportación para LightLogR.
    mkpath(joinpath(VAL, "luz"))
    csv = joinpath(VAL, "luz", "$(id).csv")
    CSV.write(csv, DataFrame(Datetime = Dates.format.(st.start, "yyyy-mm-ddTHH:MM:SS"),
                             luz_valida = [e == "válido" ? x : missing for (e, x) in zip(st.estado, st.lux_media)],
                             luz_vigilia = [e in VIGILIA ? x : missing for (e, x) in zip(st.estado, st.lux_media)]))
    dias_ll = tempname() * ".csv"; glob_ll = tempname() * ".csv"
    try
        run(pipeline(Cmd(`Rscript lightlogr_metricas.R $csv $dias_ll $glob_ll`; dir = RENV); stdout = devnull, stderr = devnull))
        ll = CSV.read(dias_ll, DataFrame; missingstring = "NA"); lg = CSV.read(glob_ll, DataFrame; missingstring = "NA")
        # 2. Brincadiano, con las mismas definiciones.
        d1 = light_daily(st, noches; thr = (100, 1000))
        d2 = light_daily(st, noches; thr = (10, 100))          # timing de la luz ≥ 100
        o = Date(floor(first(st.start), Day))
        d1.dia_desde_inicio = [Dates.value(d - o) for d in d1.dia]
        nos = innerjoin(select(d1, :dia_desde_inicio, :min_sobre_bajo, :min_sobre_alto),
                        select(DataFrame(dia_desde_inicio = d1.dia_desde_inicio, p = d2.primera_alta_h,
                                         u = d2.ultima_alta_h, m = d2.hora_media_alta_h), :), on = :dia_desde_inicio)
        j = innerjoin(nos, ll, on = :dia_desde_inicio)
        pr, ul, me = num(j.primera_100_h), num(j.ultima_100_h), num(j.media_100_h)
        okt = isfinite.(j.p) .& isfinite.(pr)
        luz = DataFrame(start = st.start, nonwear = falses(nrow(st)),
                        lux = [e in VIGILIA ? x : NaN for (e, x) in zip(st.estado, st.lux_media)])
        iv = is_iv(luz, Set(Date.(st.start)); col = :lux)
        return (id = id, dias = nrow(j),
                dif_min100 = maxabs(j.min_sobre_bajo .- j.min_100), dif_min1000 = maxabs(j.min_sobre_alto .- j.min_1000),
                dias_timing = count(okt), dif_primera_min = 60maxabs(j.p[okt] .- pr[okt]),
                dif_ultima_min = 60maxabs(j.u[okt] .- ul[okt]), dif_media_min = 60maxabs(j.m[okt] .- me[okt]),
                IS_nuestro = iv.IS, IS_lightlogr = Float64(lg.IS[1]), IV_nuestro = iv.IV, IV_lightlogr = Float64(lg.IV[1]))
    finally
        rm(dias_ll; force = true); rm(glob_ll; force = true)
    end
end

function tabla_tex(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/94_validar_luz_lightlogr.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X *{$(nrow(df))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(replace.(df.id, "_" => "\\_"), " & "), " \\\\")
    println(io, "\\midrule")
    f(n, g) = println(io, n, " & ", join([g(r) for r in eachrow(df)], " & "), " \\\\")
    f("Días comparados", r -> string(r.dias))
    f("Minutos \$\\ge\$ 100 / \$\\ge\$ 1000: diferencia máxima", r -> fmt(r.dif_min100; d = 0) * " / " * fmt(r.dif_min1000; d = 0))
    f("Días con luz \$\\ge\$ 100 (timing)", r -> string(r.dias_timing))
    f("Primera / última hora \$\\ge\$ 100: diferencia máxima (min)", r -> fmt(r.dif_primera_min; d = 0) * " / " * fmt(r.dif_ultima_min; d = 0))
    f("Hora media \$\\ge\$ 100: diferencia máxima (min)", r -> fmt(r.dif_media_min; d = 1))
    f("IS de la luz: \\Brincadiano{} / LightLogR", r -> fmt(r.IS_nuestro) * " / " * fmt(r.IS_lightlogr))
    f("IV de la luz: \\Brincadiano{} / LightLogR", r -> fmt(r.IV_nuestro) * " / " * fmt(r.IV_lightlogr))
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

if !("--confirmo" in ARGS)
    println("Contraste R4 (LightLogR). Lanza R: añade --confirmo.")
    exit(0)
end
ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "luz_60s.arrow")), readdir(INTERIM)))
res = DataFrame(map(comparar, ids))
mkpath(OUT)
CSV.write(joinpath(OUT, "luz_lightlogr.csv"), res)
write(joinpath(OUT, "luz_lightlogr.tex"), tabla_tex(res))
for x in eachrow(res)
    @printf("%-12s %d días: |Δ| min≥100 %.0f, min≥1000 %.0f; timing ≥100 (%d días): primera %.0f min, última %.0f min, media %.1f min; IS %.3f/%.3f, IV %.3f/%.3f\n",
            x.id, x.dias, x.dif_min100, x.dif_min1000, x.dias_timing, x.dif_primera_min, x.dif_ultima_min, x.dif_media_min,
            x.IS_nuestro, x.IS_lightlogr, x.IV_nuestro, x.IV_lightlogr)
end
println("Tablas en ", relpath(OUT, ROOT))
