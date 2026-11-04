# 17 · Inventario automático: las filas del embudo que calcula el pipeline.
#
# Qué hace
#   Rellena las dos últimas filas de la tabla-embudo del informe (`tab:embudo`, ficha DAT-02):
#   «con señal RAW legible» y «con días válidos suficientes». Las cinco primeras las aporta el
#   grupo (cohorte, `.gt3x`, HOMA-IR, Tanner, misma visita) y no se calculan aquí. Además cuenta,
#   por modelo y firmware (campos de la lista blanca), cuántos registros hay de cada generación
#   (PRO-06) y cuántos tienen L5 estimable. No recalcula nada: lee `data/interim/<ID>/` (01) y
#   las tablas `qc_resumen.csv` (04) y `ritmo_resumen.csv` (06).
#
# Criterios (parámetros, declarados)
#   · Señal RAW legible: existe `diagnostics.arrow`, 0 registros con checksum incorrecto y 0 bytes
#     truncados al final. Un archivo que falla en el lote (16) cuenta como no legible.
#   · Días suficientes: `dias_validos` ≥ umbral (≥ 10 h de uso por día, regla de 04). El umbral
#     principal es 4 días (`--umbral`); la tabla da también 3, 5, 6 y 7 como sensibilidad.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/17_inventario.jl                         # todos los de data/interim
#   julia --project=scripts scripts/17_inventario.jl --cohorte COHORTE_A        # etiqueta única
#   julia --project=scripts scripts/17_inventario.jl --cohortes ids_cohorte.csv   # CSV id,cohorte
#   julia --project=scripts scripts/17_inventario.jl --umbral 5
#
# Depende de
#   01 (o el lote 16) y 04; el 06 es opcional (L5 estimable).
#
# Salida (no se versiona)
#   outputs/inventario/inventario_ids.csv       una fila por registro (solo el ID)
#   outputs/inventario/inventario_cohorte.csv   recuentos por cohorte (agregados)
#   outputs/inventario/inventario_modelos.csv   recuentos por cohorte, modelo y firmware
#   outputs/inventario/inventario_embudo.tex    filas del `tab:embudo` calculadas por el pipeline
#
# Consola: solo recuentos agregados y, para los que fallan, el ID.

using Brincadiano.GT3X
using CSV, DataFrames, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const SAL = joinpath(ROOT, "outputs", "inventario")
const UMBRALES = [3, 4, 5, 6, 7]

function opcion(flag, defecto)
    i = findfirst(==(flag), ARGS)
    return i === nothing ? defecto : ARGS[i + 1]
end
umbral = parse(Int, opcion("--umbral", "4"))
umbral in UMBRALES || (UMBRALES = sort(unique([UMBRALES; umbral])))
etiqueta = opcion("--cohorte", "sin_cohorte")
f_cohortes = opcion("--cohortes", "")

# Un registro del lote tiene meta.arrow (leído) o lote_estado.csv (falló en el 16); así no cuentan
# data/interim/agd/ (script 23) ni otras carpetas que no son registros.
ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "meta.arrow")) || isfile(joinpath(INTERIM, d, "lote_estado.csv")),
                  readdir(INTERIM)))
isempty(ids) && error("No hay registros en data/interim/")

cohorte = Dict(id => etiqueta for id in ids)
if !isempty(f_cohortes)
    mapa = CSV.read(f_cohortes, DataFrame; types = String)
    for r in eachrow(mapa)
        haskey(cohorte, r.id) && (cohorte[r.id] = r.cohorte)
    end
end

tabla_opcional(ruta) = isfile(ruta) ? CSV.read(ruta, DataFrame) : nothing
qc = tabla_opcional(joinpath(ROOT, "outputs", "qc", "tables", "qc_resumen.csv"))
ritmo = tabla_opcional(joinpath(ROOT, "outputs", "circadian", "tables", "ritmo_resumen.csv"))
fila(df, id) = df === nothing ? nothing : (i = findfirst(==(id), df.id); i === nothing ? nothing : df[i, :])

# Arrow viene con Brincadiano.GT3X: se lee solo meta y diagnostics, sin cargar la aceleración.
leer(dir, nombre) = DataFrame(GT3X.Arrow.Table(joinpath(dir, nombre * ".arrow")); copycols = true)

filas = NamedTuple[]
for id in ids
    dir = joinpath(INTERIM, id)
    legible = false; chk = missing; trunc = missing; modelo = missing; fw = missing
    if isfile(joinpath(dir, "diagnostics.arrow")) && isfile(joinpath(dir, "meta.arrow"))
        d = leer(dir, "diagnostics")
        diag = Dict(zip(d.key, d.value))
        chk = tryparse(Int, get(diag, "checksum_errors", "")); trunc = tryparse(Int, get(diag, "truncated_tail_bytes", ""))
        m = leer(dir, "meta")[1, :]
        modelo = m.device_type; fw = m.firmware
        legible = chk == 0 && trunc == 0
    end
    q = fila(qc, id); r = fila(ritmo, id)
    push!(filas, (; id, cohorte = cohorte[id], modelo, firmware = fw, legible,
        checksum_errores = chk, bytes_truncados = trunc,
        dias = q === nothing ? missing : q.dias,
        dias_validos = q === nothing ? missing : q.dias_validos,
        L5_estimable = r === nothing ? missing : r.L5_estimable))
end
df = DataFrame(filas)
for u in UMBRALES
    df[!, Symbol("dias_ge$u")] = [ismissing(d) ? false : d >= u for d in df.dias_validos]
end
df.analizable = df.legible .& df[!, Symbol("dias_ge$umbral")]

mkpath(SAL)
CSV.write(joinpath(SAL, "inventario_ids.csv"), df)

por_cohorte = combine(groupby(df, :cohorte),
    :id => length => :registros,
    :legible => sum => :legibles,
    [Symbol("dias_ge$u") for u in UMBRALES] .=> sum .=> [Symbol("dias_ge$u") for u in UMBRALES],
    :analizable => sum => :analizables,
    :L5_estimable => (x -> count(skipmissing(x))) => :L5_estimable)
CSV.write(joinpath(SAL, "inventario_cohorte.csv"), por_cohorte)

modelos = combine(groupby(df, [:cohorte, :modelo, :firmware]), :id => length => :registros)
CSV.write(joinpath(SAL, "inventario_modelos.csv"), modelos)

# Filas del embudo, en el formato de tab:embudo (columnas = cohortes, en orden alfabético).
cs = sort(por_cohorte.cohorte)
cel(c, col) = por_cohorte[findfirst(==(c), por_cohorte.cohorte), col]
open(joinpath(SAL, "inventario_embudo.tex"), "w") do io
    println(io, "% Generado por scripts/17_inventario.jl (cohortes: ", join(cs, ", "), "). No editar a mano.")
    println(io, "con señal RAW legible & ", join([string(cel(c, :legibles)) for c in cs], " & "), " \\\\")
    println(io, "con días válidos suficientes (≥ $umbral) & ", join([string(cel(c, :analizables)) for c in cs], " & "), " \\\\")
end

println("Registros: ", nrow(df), "   cohortes: ", join(cs, ", "), "   umbral: ≥ $umbral días válidos")
for r in eachrow(por_cohorte)
    @printf("  %-14s legibles %d/%d   analizables %d   L5 estimable %d\n", r.cohorte, r.legibles, r.registros, r.analizables, r.L5_estimable)
end
for r in eachrow(df[.!df.legible, :])
    println("  no legible: ", r.id)
end
println("Salida en ", relpath(SAL, ROOT))
