# 24 · Auditoría de los config.csv de GGIR enviados por el grupo.
#
# Qué hace
#   Lee cada config_*.csv de data/raw/ y sus subcarpetas (formato de GGIR: argument, value, context), comprueba que
#   lo es y los pone uno al lado del otro: argumento, valor en cada archivo, contexto (params_*) y
#   si difieren. Marca los parámetros que deciden si un resultado de GGIR es comparable con el
#   pipeline (epoch, calibración, no uso, imputación, días válidos, umbrales, colocación, sueño,
#   luz) y los compara con lo que hace Brincadiano. Las rutas (datadir, outputdir, archivos de
#   limpieza…) se enmascaran: pueden llevar nombres de carpetas o de personas.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/24_config_ggir.jl
#
# Salida (no se versiona)
#   outputs/agd/tables/config_ggir.csv          todos los argumentos, valor por archivo
#   outputs/agd/tables/config_ggir_clave.csv    los parámetros que deciden la comparabilidad
#
# Consola: versión de GGIR, número de archivos procesados (f1) y diferencias.

using CSV, DataFrames, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = joinpath(ROOT, "data", "raw")
const OUT = joinpath(ROOT, "outputs", "agd", "tables")

# Parámetro de GGIR => (qué hace Brincadiano, referencia en el informe o el código)
const CLAVE = [
    "GGIRversion" => ("—", ""),
    "f1" => ("—", "número de archivos procesados en la tanda"),
    "dataFormat" => ("raw (.gt3x)", "01"),
    "windowsizes" => ("epochs de 1, 5 y 60 s; no uso en 15/60 min", "05, 04"),
    "do.cal" => ("autocalibración de van Hees sin temperatura", "05"),
    "spherecrit" => ("0,3", "Metrics"), "minloadcrit" => ("—", ""),
    "acc.metric" => ("ENMO (y MAD)", "05"),
    "nonwear_approach" => ("2023 (15/60 min, 13 mg / 50 mg, ≥ 2 ejes)", "04"),
    "imputeTimegaps" => ("no: los huecos (idle sleep) se marcan, no se rellenan", "01"),
    "do.imp" => ("no: no uso ausente; L5* solo como sensibilidad", "06"),
    "strategy" => ("todo el registro", "04"), "ndayswindow" => ("—", ""),
    "includedaycrit" => ("≥ 10 h de uso por día", "04"),
    "includedaycrit.part5" => ("—", ""), "includenightcrit" => ("—", ""),
    "threshold.lig" => ("63,3 mg (Hildebrand, cadera)", "ACT-01"),
    "threshold.mod" => ("142,6 mg", "ACT-01"), "threshold.vig" => ("464,6 mg", "ACT-01"),
    "boutcriter.mvpa" => ("episodios vigorosos: pausas ≤ 2 s, ≥ 3 s (no bouts de GGIR)", "06"),
    "sensor.location" => ("cadera (derecha, según PRO-04 para COHORTE_A)", "PRO-04"),
    "HASPT.algo" => ("HorAngle (cadera), eje longitudinal Y", "08"),
    "HorAngle_threshold" => ("60°", "08"), "longitudinal_axis" => ("2 (Y)", "08"),
    "HASIB.algo" => ("vanHees2015 (5°, 5 min)", "08"),
    "anglethreshold" => ("5", "08"), "timethreshold" => ("5", "08"),
    "desiredtz" => ("hora local del dispositivo", "01"),
    "winhr" => ("L5 y M10", "06"), "M5L5res" => ("—", ""),
    "IVIS_windowsize_minutes" => ("60 min", "06"), "IVIS_acc_threshold" => ("—", ""),
    "LUXthresholds" => ("100 y 1000 brutos (lux_raw, sin convertir)", "09"),
    "LUX_cal_constant" => ("sin conversión a lux (no resuelta)", "Light"),
    "LUX_cal_exponent" => ("—", ""), "LUX_day_segments" => ("ventanas ancladas a la SPT", "09"),
    "save_ms5rawlevels" => ("—", ""), "part5_agg2_60seconds" => ("—", ""),
]

es_ruta(v) = occursin(r"[/\\]|^[A-Za-z]:", v)
enmascarar(v) = ismissing(v) ? missing : (es_ruta(v) ? "<ruta enmascarada>" : v)

rutas = Dict{String,String}()
for (dir, _, fs) in walkdir(RAW), f in fs
    startswith(f, "config") && endswith(lowercase(f), ".csv") && (rutas[f] = joinpath(dir, f))
end
archivos = sort(collect(keys(rutas)))
isempty(archivos) && error("No hay config*.csv en data/raw/")
tablas = Dict{String,DataFrame}()
for f in archivos
    t = CSV.read(rutas[f], DataFrame; types = String)
    names(t) == ["argument", "value", "context"] || error("$f no tiene el formato de config.csv de GGIR")
    tablas[f] = t
end
nombre(f) = replace(splitext(f)[1], "config_" => "")
args = sort(unique(vcat([t.argument for t in values(tablas)]...)))
df = DataFrame(argumento = args)
for f in archivos
    d = Dict(zip(tablas[f].argument, coalesce.(tablas[f].value, "")))
    df[!, nombre(f)] = [haskey(d, a) ? enmascarar(d[a]) : missing for a in args]
end
ctx = Dict{String,String}()
for t in values(tablas), r in eachrow(t)
    ismissing(r.context) || (ctx[r.argument] = r.context)
end
df.contexto = [get(ctx, a, "") for a in args]
cols = nombre.(archivos)
df.difiere = [length(unique(skipmissing([r[c] for c in cols]))) > 1 || any(ismissing, [r[c] for c in cols]) for r in eachrow(df)]
mkpath(OUT)
CSV.write(joinpath(OUT, "config_ggir.csv"), df)

clave = DataFrame(argumento = first.(CLAVE), brincadiano = first.(last.(CLAVE)), donde = last.(last.(CLAVE)))
clave = leftjoin(clave, df[:, [:argumento; Symbol.(cols); :contexto]]; on = :argumento, order = :left)
CSV.write(joinpath(OUT, "config_ggir_clave.csv"), clave)

for f in archivos
    d = Dict(zip(tablas[f].argument, coalesce.(tablas[f].value, "")))
    @printf("%-22s GGIR %-6s  f1 = %-5s  %d argumentos\n", f, get(d, "GGIRversion", "?"), get(d, "f1", "?"), nrow(tablas[f]))
end
println("Argumentos que difieren o faltan en alguno: ", count(df.difiere))
println("Salida en ", relpath(OUT, ROOT))
