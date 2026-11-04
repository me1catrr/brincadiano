# 03 · Exploración visual de un registro.
#
# Qué hace
#   Dibuja la señal en bruto de un registro ya leído: la semana completa, el
#   actograma de 24 h y de 48 h, los histogramas de los tres ejes y cuatro ventanas
#   cortas elegidas de forma determinista. Sirve para mirar la señal antes del
#   control de calidad; no clasifica uso, sueño ni intensidad. Las figuras de días
#   van a outputs/qc/figures/ junto a las del script 04, que son otras (el 04 superpone
#   los tramos inmóviles y se llama <ID>_dias24h_qc).
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/03_explorar_sujeto.jl [ID ...]   # sin ID, solo 6005_V00
#   julia --project=scripts scripts/03_explorar_sujeto.jl --fixture  # ilustración versionable
#   julia --project=scripts scripts/03_explorar_sujeto.jl --solo-vista 6021_V00
#       solo la semana completa y los histogramas (H25: las únicas de otros registros
#       que entran en el informe; así no se generan figuras que no se usan)
#
# Depende de
#   01 (data/interim/<ID>/meta.arrow y el resto del registro). No usa 02 ni 04.
#
# Lo usan
#   Ningún script posterior. `make sync-assets` copia las figuras individuales al
#   informe. Con --fixture la ilustración se escribe ya en el árbol versionado del
#   capítulo 2, a partir del archivo público de ActiGraph
#   (test/fixtures/legacy_sample.gt3x, MIT; no es de ningún participante).
#
# Entrada
#   data/interim/<ID>/*.arrow
#   o, con --fixture, test/fixtures/legacy_sample.gt3x (se lee directo, sin interim)
#
# Salida (figuras individuales: no se versionan; PDF y PNG)
#   outputs/gt3x/figures/<ID>_vista_semana
#   outputs/gt3x/figures/<ID>_histogramas
#   outputs/gt3x/figures/<ID>_ventana_maxima      minuto de mayor |‖a‖ − 1 g|, muestreo completo
#   outputs/gt3x/figures/<ID>_ventana_tipica      percentil 75 de esa media entre las 08:00 y las 21:00
#   outputs/gt3x/figures/<ID>_ventana_reposo      minuto de menor media
#   outputs/gt3x/figures/<ID>_ventana_idle        entrada al primer idle sleep de al menos 10 min
#   outputs/qc/figures/<ID>_dias24h
#   outputs/qc/figures/<ID>_dias48h
#   Ventanas de 20 s (30 s en la de idle sleep). Si un criterio no tiene candidato, esa figura se omite.
# Salida versionada (--fixture)
#   report/assets/02_acelerometria/report/figures/senal-raw-ejemplo.pdf
#
# Consola: solo el ID y recuentos de muestras, huecos y registros de luz.

using Brincadiano.GT3X, Brincadiano.Viz
using CairoMakie, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_GT3X = joinpath(ROOT, "outputs", "gt3x", "figures")
const OUT_QC = joinpath(ROOT, "outputs", "qc", "figures")

function guardar_fig(fig, dir, nombre, id, funcion)
    mkpath(dir)
    guardar(joinpath(dir, nombre), fig;
            procedencia = (; tipo = :datos, registros = [id], funcion, entradas = ["data/interim/$id/"]))
    println("  ", relpath(joinpath(dir, nombre), ROOT), ".{pdf,png}")
end

# Inicio de una ventana de `dur` centrada en el segundo más intenso del minuto `m0`.
function centrar(r, nrm, m0, dur)
    sel = (r.acc.time .>= m0) .& (r.acc.time .< m0 + Minute(1))
    t = r.acc.time[sel]; d = abs.(nrm[sel] .- 1)
    s = binned_summary(t, d, Second(1); origin = m0, stop = m0 + Second(59))
    k = argmax(replace(s.vmean, NaN => -Inf))
    return clamp(s.start[k] - dur / 2, m0, m0 + Minute(1) - dur)
end

function explorar(id; solo_vista::Bool = false)
    dir = joinpath(INTERIM, id)
    isfile(joinpath(dir, "meta.arrow")) || error("$id: no está en data/interim/ (ejecuta scripts/01_read_gt3x.jl)")
    r = read_interim(dir)
    fs = r.info.sample_rate
    println(id, ": ", nrow(r.acc), " muestras, ", nrow(r.gaps), " huecos, ", nrow(r.lux), " registros de luz")

    guardar_fig(plot_overview(r), OUT_GT3X, "$(id)_vista_semana", id, "plot_overview")
    guardar_fig(plot_histograms(r), OUT_GT3X, "$(id)_histogramas", id, "plot_histograms")
    solo_vista && return
    guardar_fig(plot_days24h(r), OUT_QC, "$(id)_dias24h", id, "plot_days24h")
    guardar_fig(plot_days24h(r; double = true), OUT_QC, "$(id)_dias48h", id, "plot_days24h")

    nrm = vector_norm(r.acc)
    m = binned_summary(r.acc.time, abs.(nrm .- 1), Minute(1))
    full = m.n .== 60fs
    dia = full .& (8 .<= hour.(m.start) .< 21)
    dur = Second(20)
    if any(full)
        i_max = findall(full)[argmax(m.vmean[full])]
        i_min = findall(full)[argmin(m.vmean[full])]
        guardar_fig(plot_raw_window(r, centrar(r, nrm, m.start[i_max], dur), dur), OUT_GT3X, "$(id)_ventana_maxima", id, "plot_raw_window")
        guardar_fig(plot_raw_window(r, m.start[i_min] + Second(20), dur), OUT_GT3X, "$(id)_ventana_reposo", id, "plot_raw_window")
    end
    if any(dia)
        v = m.vmean[dia]; idx = findall(dia)
        i_tip = idx[sortperm(v)[ceil(Int, 0.75 * length(v))]]
        guardar_fig(plot_raw_window(r, centrar(r, nrm, m.start[i_tip], dur), dur), OUT_GT3X, "$(id)_ventana_tipica", id, "plot_raw_window")
    end
    largos = findall((r.gaps.kind .== "idle_sleep") .& (r.gaps.stop .- r.gaps.start .>= Minute(10)))
    if !isempty(largos)
        g0 = r.gaps.start[first(largos)]
        guardar_fig(plot_raw_window(r, g0 - Second(25), Second(30)), OUT_GT3X, "$(id)_ventana_idle", id, "plot_raw_window")
    end
end

function fixture()
    r = read_gt3x(joinpath(ROOT, "test", "fixtures", "legacy_sample.gt3x"))
    nrm = vector_norm(r.acc)
    # 20 s con muestreo completo y más movimiento del registro de demostración (84 s).
    s = binned_summary(r.acc.time, abs.(nrm .- 1), Second(20); origin = first(r.acc.time))
    k = argmax(ifelse.(s.n .== 20r.info.sample_rate, s.vmean, -Inf))
    fig = plot_raw_window(r, s.start[k], Second(20))
    dir = joinpath(ROOT, "report", "assets", "02_acelerometria", "report", "figures")
    mkpath(dir)
    guardar(joinpath(dir, "senal-raw-ejemplo"), fig; png = false,
            procedencia = (; tipo = :fixture, funcion = "plot_raw_window", entradas = ["test/fixtures/legacy_sample.gt3x"]))
    println("  ", relpath(joinpath(dir, "senal-raw-ejemplo.pdf"), ROOT))
end

if "--fixture" in ARGS
    fixture()
else
    ids = filter(a -> !startswith(a, "--"), ARGS)
    sv = "--solo-vista" in ARGS
    foreach(id -> explorar(id; solo_vista = sv), isempty(ids) ? ["6005_V00"] : ids)
end
