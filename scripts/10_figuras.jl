# 10 · Figuras de síntesis a partir de lo ya calculado.
#
# Qué hace
#   No recalcula métricas. Lee los Arrow y las tablas de los scripts anteriores y
#   dibuja las figuras que comparan registros o que el informe no puede componer
#   con una sola salida. Entra todo registro que esté en ritmo_resumen.csv y tenga
#   luz_60s.arrow. Las tres figuras de validación se omiten, sin error, si aún no
#   se han ejecutado el 92 o el 93.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/10_figuras.jl
#   Orden práctico: 01 → 04 → 05 → 06 → 08 → 09, y el 07 si se quiere el plano
#   circular. El 92 y el 93 son opcionales y exigen autorización.
#
# Depende de
#   01  capsense (actograma)
#   05  epochs_5s.arrow y epochs_60s.arrow
#   06  ritmo_resumen.csv, timing_resumen.csv, episodios_<ID>.csv
#   07  ejemplo_modelo.csv (escenario «con efecto»); si falta, el script se detiene
#       al llegar al plano circular
#   08  angulos_5s.arrow (recalcula las noches) y sueno_<ID>_noches.csv (el L5 de
#       cada noche, en el mismo orden)
#   09  luz_60s.arrow
#   92  outputs/validation/ggir_<ID>_enmo5s_cruce.csv y ggir_<ID>_nouso15min_cruce.csv
#       (opcionales)
#   93  outputs/validation/tables/ggir_sueno_noches.csv, filas config = ejeY (opcional)
#
# Lo usan
#   Solo el informe, tras `make sync-assets`. Ningún script lee estas figuras.
#
# Entrada
#   Las rutas de «Depende de». No abre data/raw/.
#
# Salida (individual; no se versiona; PDF y PNG; sin títulos: los pone el pie de LaTeX)
#   outputs/light/figures/<ID>_actograma                  caps. 12 y 15; luz, sueño, vigorosa, CAPSENSE
#   outputs/metrics/figures/dos_registros_enmo_ecdf       cap. 7; ENMO de 5 s en uso, con idle sleep = 0
#   outputs/vigorous/figures/episodios_vigorosos          cap. 8; el episodio de más dosis, a 30 Hz
#   outputs/sleep/figures/dos_registros_ventanas          cap. 9
#   outputs/circadian/figures/dos_registros_perfiles24h   cap. 10; marca L5* si el ancla no era estimable
#   outputs/timing/figures/dos_registros_diana            cap. 11
#   outputs/metabolic/figures/plano_circular_dos_registros  cap. 14
#   outputs/validation/figures/bland_altman_enmo5s        cap. 13; solo si existe el cruce del 92
#   outputs/validation/figures/no_uso_ggir                cap. 13; solo si existe el cruce del 92
#   outputs/validation/figures/ventanas_sueno_ggir        cap. 13; solo si existe la tabla del 93
#   outputs/validation/figures/bland_altman_mad5s         cap. 13; si el cruce del 92 trae MAD (H28)
#   outputs/validation/figures/imputacion_dumbbell        cap. 13; si existe la tabla del 95
#   outputs/validation/figures/intensidad_ggir            cap. 13; si existe la tabla del 96
#
# La ventana de 30 Hz del episodio vigoroso recalibra la señal del script 01 con la
# misma regla del 05; no reescribe los epochs. Consola: solo las rutas escritas.

using Brincadiano.GT3X, Brincadiano.Metrics, Brincadiano.Circadian, Brincadiano.Sleep, Brincadiano.Light, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT = joinpath(ROOT, "outputs")
const MIN_USO = 600

tabla(sub, nombre) = CSV.read(joinpath(OUT, sub, "tables", nombre), DataFrame)
# Cada figura declara su procedencia (registros, función de Viz y entradas): Viz.guardar
# escribe la ficha <nombre>.procedencia.tex que el informe muestra bajo la figura.
function salida(sub, nombre, fig; registros, funcion, entradas)
    dir = joinpath(OUT, sub, "figures"); mkpath(dir)
    guardar(joinpath(dir, nombre), fig;
            procedencia = (; tipo = :datos, registros = collect(String, registros), funcion, entradas))
    println("  ", relpath(joinpath(dir, nombre), ROOT), ".{pdf,png}")
end

ritmo = tabla("circadian", "ritmo_resumen.csv")
timing_r = tabla("timing", "timing_resumen.csv")
ids = sort(filter(id -> id in ritmo.id && isfile(joinpath(INTERIM, id, "luz_60s.arrow")), readdir(INTERIM)))

# Datos por registro.
D = Dict{String,Any}()
for id in ids
    dir = joinpath(INTERIM, id)
    r = read_interim(dir)
    e60 = load_epochs(dir, 60); e5 = load_epochs(dir, 5)
    epf = load_angles(dir)
    eje = longitudinal_axis(epf; mask = .!epf.nonwear)
    noches = sleep_nights(epf; axis = eje, nonwear = epf.nonwear, rec_start = first(r.acc.time), rec_stop = last(r.acc.time))
    sn = tabla("sleep", "sueno_$(id)_noches.csv")             # mismo orden que `noches`
    noches.L5_noche_h = sn.L5_noche_h
    noches.dia_desde_inicio = [Dates.value(d - Date(first(r.acc.time))) for d in noches.noche]
    ep = CSV.read(joinpath(OUT, "timing", "tables", "episodios_$(id).csv"), DataFrame;
                  types = Dict(:inicio => DateTime, :fin => DateTime, :medio => DateTime))
    rr = ritmo[ritmo.id .== id, :][1, :]; tr = timing_r[timing_r.id .== id, :][1, :]
    dd = daily_summary(e60); validos = Set(dd.dia[dd.minutos_uso .>= MIN_USO])
    perfil = profile24h(e60, validos; nonwear_as = rr.L5_estimable ? :missing : :rest)
    D[id] = (r = r, e60 = e60, e5 = e5, luz = load_light(dir), noches = noches, ep = ep, rr = rr, tr = tr,
             perfil = perfil, dias = length(validos))
end

# 1. Actograma multimodal por registro.
for id in ids
    d = D[id]
    salida("light", "$(id)_actograma", plot_actogram(d.e60, d.luz, d.noches, d.ep, d.r.capsense);
           registros = [id], funcion = "plot_actogram",
           entradas = ["data/interim/$id/", "outputs/sleep/tables/sueno_$(id)_noches.csv", "outputs/timing/tables/episodios_$(id).csv"])
end

# 2. Comparativas entre registros.
anc(rr) = rr.L5_estimable ? rr.L5_medio_h : rr.L5_reposo_medio_h
salida("circadian", "dos_registros_perfiles24h",
       plot_profiles_compare([(id * (D[id].rr.L5_estimable ? "" : " (L5*: no uso = reposo)"), D[id].perfil, anc(D[id].rr)) for id in ids]);
       registros = ids, funcion = "plot_profiles_compare",
       entradas = ["data/interim/<ID>/epochs_60s.arrow", "outputs/circadian/tables/ritmo_resumen.csv"])
salida("metrics", "dos_registros_enmo_ecdf",
       plot_enmo_ecdf([(id, D[id].e5.enmo_idle0[.!D[id].e5.nonwear]) for id in ids]);
       registros = ids, funcion = "plot_enmo_ecdf", entradas = ["data/interim/<ID>/epochs_5s.arrow"])
salida("timing", "dos_registros_diana",
       plot_timing_compare([(id, 2π .* D[id].ep.fase_h ./ 24, D[id].ep.dosis, 2π * D[id].tr.θ_h / 24, D[id].tr.R) for id in ids]);
       registros = ids, funcion = "plot_timing_compare",
       entradas = ["outputs/timing/tables/episodios_<ID>.csv", "outputs/timing/tables/timing_resumen.csv"])
sueno(d) = select(d.noches, :dia_semana, :libre, :metodo, :inicio_h, :fin_h, :L5_noche_h, :dia_desde_inicio)
salida("sleep", "dos_registros_ventanas", plot_sleep_raster([(id, sueno(D[id])) for id in ids]);
       registros = ids, funcion = "plot_sleep_raster",
       entradas = ["data/interim/<ID>/angulos_5s.arrow", "outputs/sleep/tables/sueno_<ID>_noches.csv"])

# 3. Actividad vigorosa: el episodio de más dosis del registro con más episodios.
idv = ids[argmax([nrow(D[id].ep) for id in ids])]
d = D[idv]
k = argmax(d.ep.dosis)
e0, e1 = d.ep.inicio[k], d.ep.fin[k] + Second(1)
ac = autocalibrate(d.r.acc, d.r.info.sample_rate)
off, sc = ac.aplicar ? (ac.offset, ac.scale) : ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
salida("vigorous", "episodios_vigorosos",
       plot_vigorous_summary(d.r, (e0 - Second(8), e1 + Second(8), e0, e1), [(id, D[id].ep, D[id].dias) for id in ids];
                             offset = off, scale = sc);
       registros = ids, funcion = "plot_vigorous_summary",
       entradas = ["data/interim/$idv/acc.arrow", "outputs/timing/tables/episodios_<ID>.csv"])

# 4. Plano circular del modelo sintético con los dos registros.
ej = tabla("metabolic", "ejemplo_modelo.csv"); ej = ej[ej.escenario .== "con efecto", :][1, :]
salida("metabolic", "plano_circular_dos_registros",
       plot_circular_plane(ej.A, ej.θ0_h, [(id, D[id].tr.R, D[id].tr.θ_h) for id in ids]);
       registros = ids, funcion = "plot_circular_plane",
       entradas = ["outputs/timing/tables/timing_resumen.csv", "outputs/metabolic/tables/ejemplo_modelo.csv"])

# 5. Validación frente a GGIR (cruces guardados por scripts/92 y 93).
val(n) = CSV.read(joinpath(OUT, "validation", n), DataFrame; types = Dict(:start => DateTime))
cr = [(id, val("ggir_$(id)_enmo5s_cruce.csv")) for id in ids if isfile(joinpath(OUT, "validation", "ggir_$(id)_enmo5s_cruce.csv"))]
isempty(cr) || salida("validation", "bland_altman_enmo5s", plot_bland_altman([(id, c.enmo, c.enmo_ggir) for (id, c) in cr]);
                      registros = first.(cr), funcion = "plot_bland_altman",
                      entradas = ["outputs/validation/ggir_<ID>_enmo5s_cruce.csv"])
if !isempty(cr) && hasproperty(cr[1][2], :mad)
    salida("validation", "bland_altman_mad5s",
           plot_bland_altman([(id, c.mad, c.mad_ggir) for (id, c) in cr]);
           registros = first.(cr), funcion = "plot_bland_altman",
           entradas = ["outputs/validation/ggir_<ID>_enmo5s_cruce.csv"])
end
impf = joinpath(OUT, "validation", "tables", "imputacion.csv")
if isfile(impf)
    imp = CSV.read(impf, DataFrame)
    salida("validation", "imputacion_dumbbell", plot_imputation_dumbbell(imp);
           registros = unique(imp.id), funcion = "plot_imputation_dumbbell",
           entradas = ["outputs/validation/tables/imputacion.csv"])
end
intf = joinpath(OUT, "validation", "tables", "ggir_intensidad_dias.csv")
if isfile(intf)
    di = CSV.read(intf, DataFrame)
    bandas = ["sedentario", "ligera", "moderada", "vigorosa"]
    mediaf(v) = (x = filter(isfinite, Float64.(collect(skipmissing(v)))); isempty(x) ? NaN : sum(x) / length(x))
    filas = NamedTuple[]
    for id in unique(di.id), b in bandas
        s = di[di.id .== id, :]
        push!(filas, (id = id, banda = b, nuestro = mediaf(s[!, b * "_nuestro"]), ggir = mediaf(s[!, b * "_ggir"])))
    end
    salida("validation", "intensidad_ggir", plot_intensity_ggir(DataFrame(filas));
           registros = unique(di.id), funcion = "plot_intensity_ggir",
           entradas = ["outputs/validation/tables/ggir_intensidad_dias.csv"])
end
nw = [(id, val("ggir_$(id)_nouso15min_cruce.csv")) for id in ids if isfile(joinpath(OUT, "validation", "ggir_$(id)_nouso15min_cruce.csv"))]
isempty(nw) || salida("validation", "no_uso_ggir", plot_agreement_timeline([(id, c.start, c.nuestro, c.ggir) for (id, c) in nw]);
                      registros = first.(nw), funcion = "plot_agreement_timeline",
                      entradas = ["outputs/validation/ggir_<ID>_nouso15min_cruce.csv"])
gs = joinpath(OUT, "validation", "tables", "ggir_sueno_noches.csv")
if isfile(gs)
    g = CSV.read(gs, DataFrame); g = g[g.config .== "ejeY", :]
    datos = map(ids) do id
        s = sueno(D[id]); gi = g[g.id .== id, :]
        m = Dict(n - 1 => (a, b) for (n, a, b) in zip(gi.noche, gi.inicio_ggir, gi.fin_ggir))
        s.inicio_ggir = [get(m, x, (NaN, NaN))[1] for x in s.dia_desde_inicio]
        s.fin_ggir = [get(m, x, (NaN, NaN))[2] for x in s.dia_desde_inicio]
        (id, s)
    end
    salida("validation", "ventanas_sueno_ggir", plot_sleep_raster(datos; ggir = true);
           registros = ids, funcion = "plot_sleep_raster",
           entradas = ["outputs/validation/tables/ggir_sueno_noches.csv", "outputs/sleep/tables/sueno_<ID>_noches.csv"])
end
