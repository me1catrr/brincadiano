# 04 · Control de calidad por archivo.
#
# Qué hace
#   Para cada registro calcula el informe de Brincadiano.QC: cobertura del muestreo,
#   huecos, saturación, luz, contacto (CAPSENSE) y tramos inmóviles. El no uso sigue
#   GGIR 3.3.9 (método «2023»: ventanas de 15 min dentro de 60 min, rango < 50 mg y
#   desviación típica < 13 mg en al menos dos ejes), con la adaptación de este
#   proyecto al idle sleep (un tramo sin muestras cuenta como quieto). Un día es
#   válido si el uso alcanza el umbral (10 h por defecto). El actograma de esta
#   figura lleva los tramos encima; el del script 03 (<ID>_dias24h, sin sufijo) no.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/04_qc.jl [ID ...]            # sin ID, todos los que tengan meta.arrow
#   julia --project=scripts scripts/04_qc.jl --horas-validas 8   # umbral de día válido, en horas
#
# Depende de
#   01 (data/interim/<ID>/*.arrow). No lee las tablas del 02 ni las figuras del 03.
#
# Lo usan
#   05 y 08, que cargan data/interim/<ID>/qc_tramos.arrow para marcar el no uso en
#   epochs y ángulos. El 06 hereda esa marca a través de los epochs del 05. El 92
#   vuelve a calcular las ventanas de 15 min para compararlas con GGIR; no lee estas
#   tablas. `make sync-assets` copia tablas y figuras al capítulo 7.
#
# Entrada
#   data/interim/<ID>/{acc,gaps,lux,capsense,meta}.arrow
#
# Salida (individual; no se versiona)
#   data/interim/<ID>/qc_tramos.arrow              intermedio: tramos inmóviles y su etiqueta
#   outputs/qc/tables/qc_resumen.csv               una fila por registro (se reescribe entera)
#   outputs/qc/tables/qc_resumen.tex
#   outputs/qc/tables/qc_<ID>_dias.csv             uso por día y si el día es válido
#   outputs/qc/tables/qc_<ID>_tramos.csv           cada tramo inmóvil
#   outputs/qc/figures/<ID>_dias24h_qc.{pdf,png}   actograma con los tramos
#
# Consola: solo el ID, horas, porcentaje muestreado, tramos y días válidos.

using Brincadiano.GT3X, Brincadiano.QC, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_T = joinpath(ROOT, "outputs", "qc", "tables")
const OUT_F = joinpath(ROOT, "outputs", "qc", "figures")

function argval(flag, default)
    i = findfirst(==(flag), ARGS)
    return i === nothing ? default : parse(Float64, ARGS[i + 1])
end
const HORAS_VALIDAS = argval("--horas-validas", 10.0)

# Coma decimal (convención del informe): «{,}» evita el espacio que LaTeX añade tras la coma en modo matemático.
fmt(x; d = 1) = ismissing(x) || (x isa AbstractFloat && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»

function qc_id(id)
    r = read_interim(joinpath(INTERIM, id))
    q = qc_report(r; valid_hours = HORAS_VALIDAS)
    mkpath(OUT_T); mkpath(OUT_F)
    CSV.write(joinpath(OUT_T, "qc_$(id)_dias.csv"), q.dias)
    CSV.write(joinpath(OUT_T, "qc_$(id)_tramos.csv"), q.tramos)
    save_runs(joinpath(INTERIM, id), q.tramos)
    fig = plot_days24h(r; runs = q.tramos)
    guardar(joinpath(OUT_F, "$(id)_dias24h_qc"), fig;
            procedencia = (; tipo = :datos, registros = [id], funcion = "plot_days24h", entradas = ["data/interim/$id/"]))

    g = q.huecos
    gh(k, col) = (i = findfirst(==(k), g.kind); i === nothing ? 0.0 : g[i, col])
    tr = q.tramos
    th(lab) = (sel = tr.etiqueta .== lab; (count(sel), sum(tr.horas[sel]; init = 0.0)))
    nu, ic = th("no uso probable"), th("inmóvil con contacto")
    fila = (id = id, horas = q.muestreo.duracion_h, pct_muestreado = q.muestreo.pct_tiempo_muestreado,
            fs_efectiva = q.muestreo.fs_efectiva, pct_seg_completos = q.muestreo.pct_segundos_completos,
            idle_n = Int(gh("idle_sleep", :n)), idle_h = gh("idle_sleep", :horas),
            idle_max_h = gh("idle_sleep", :max_h),
            otros_huecos_n = sum(g.n[g.kind .!= "idle_sleep"]; init = 0),
            pct_saturadas = q.saturacion.pct_saturadas,
            luz_registros = q.luz.registros, luz_pct_ceros = q.luz.pct_ceros,
            luz_max_bruto = q.luz.max_bruto, luz_pct_en_idle = q.luz.pct_en_idle,
            contacto_ref_pct = 100 * q.contacto_ref,
            no_uso_n = nu[1], no_uso_h = nu[2], inmovil_contacto_n = ic[1], inmovil_contacto_h = ic[2],
            dias = nrow(q.dias), dias_validos = count(q.dias.valido))
    @printf("%-12s %6.1f h  muestreado %5.1f %%  tramos: %d no uso (%.1f h), %d con contacto (%.1f h)  días válidos %d/%d\n",
            id, fila.horas, fila.pct_muestreado, nu[1], nu[2], ic[1], ic[2], fila.dias_validos, fila.dias)
    return fila
end

function tabla_tex(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/04_qc.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}l *{9}{>{\\raggedleft\\arraybackslash}X}@{}}")
    println(io, "\\toprule")
    println(io, "ID & Horas & Muestreado (\\%) & \\emph{Idle sleep} (n / h) & Saturadas (\\%) & ",
            "Luz: ceros (\\%) & Luz en \\emph{idle} (\\%) & Contacto en uso (\\%) & No uso probable (n / h) & Días válidos \\\\")
    println(io, "\\midrule")
    for r in eachrow(df)
        println(io, replace(r.id, "_" => "\\_"), " & ", fmt(r.horas), " & ", fmt(r.pct_muestreado), " & ",
                r.idle_n, " / ", fmt(r.idle_h), " & ", fmt(r.pct_saturadas; d = 3), " & ",
                fmt(r.luz_pct_ceros), " & ", fmt(r.luz_pct_en_idle), " & ", fmt(r.contacto_ref_pct), " & ",
                r.no_uso_n, " / ", fmt(r.no_uso_h), " & ", r.dias_validos, "/", r.dias, " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

ids = filter(a -> !startswith(a, "--") && tryparse(Float64, a) === nothing, ARGS)
isempty(ids) && (ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "meta.arrow")), readdir(INTERIM))))
res = DataFrame(map(qc_id, ids))
CSV.write(joinpath(OUT_T, "qc_resumen.csv"), res)
write(joinpath(OUT_T, "qc_resumen.tex"), tabla_tex(res))
println("Resumen en ", relpath(joinpath(OUT_T, "qc_resumen.{csv,tex}"), ROOT),
        " (umbral de día válido: ", HORAS_VALIDAS, " h)")
