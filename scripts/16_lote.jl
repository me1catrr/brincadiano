# 16 · Lectura por lotes, reanudable.
#
# Qué hace
#   Lee cada data/raw/*.gt3x con el script 01, uno por uno y en un proceso aparte, de modo que
#   un archivo que falle no tira el lote y que una interrupción se puede reanudar. Un registro
#   cuenta como leído si existe `diagnostics.arrow`, que `GT3X.write_interim` escribe el último;
#   un directorio a medias se rehace con `--force` solo para ese ID. Lleva un registro de estado
#   por ID. Con `--analisis` lanza después, una sola vez y con todos los IDs leídos, los pasos
#   02, 04, 05, 06, 08, 09, 12 y 13, porque esos scripts reescriben enteras sus tablas
#   `*_resumen.csv` y por tanto no se pueden reanudar registro a registro.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/16_lote.jl                    # lee los que falten
#   julia --project=scripts scripts/16_lote.jl 6005_V00 6021_V00  # solo esos IDs
#   julia --project=scripts scripts/16_lote.jl --proyecto COHORTE_A  # solo un proyecto (Brincadiano.Rutas)
#   julia --project=scripts scripts/16_lote.jl --analisis         # además, pasos 02–13 sobre todos
#   julia --project=scripts scripts/16_lote.jl --rehacer          # relee aunque estén leídos
#   julia --project=scripts scripts/16_lote.jl --simular          # dice qué haría, sin ejecutar nada
#
# Depende de
#   01 (y, con --analisis, 02, 04, 05, 06, 08, 09, 12 y 13). Lee data/raw/ en solo lectura.
#
# Salida (no se versiona)
#   data/interim/<ID>/lote_estado.csv   paso, estado (ok / error), fecha, segundos, mensaje
#   outputs/lote/lote_pasos.csv         estado de los pasos de análisis de la última tanda
#
# Consola: solo el ID, el paso y el estado. Del error se guarda la primera línea (≤ 160
# caracteres) del tipo de excepción; nunca la salida completa del subproceso.

using Brincadiano.Rutas: Rutas, raw_files, id_registro
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const INTERIM = joinpath(ROOT, "data", "interim")
const SALIDA = joinpath(ROOT, "outputs", "lote")
const PASOS_ANALISIS = ["02_gt3x_summary", "04_qc", "05_metricas", "06_circadiano",
                        "08_sueno", "09_luz", "12_espectro", "13_espectro_movimiento"]

simular = "--simular" in ARGS
rehacer = "--rehacer" in ARGS
analisis = "--analisis" in ARGS
i_proy = findfirst(==("--proyecto"), ARGS)
proyecto = i_proy === nothing ? nothing : ARGS[i_proy + 1]
pedidos = filter(a -> !startswith(a, "--") && a != proyecto, ARGS)

leido(id) = isfile(joinpath(INTERIM, id, "diagnostics.arrow"))

"Ejecuta un script en un subproceso; devuelve (ok, mensaje corto) sin volcar su salida."
function ejecutar(script, args = String[])
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "scripts")) $(joinpath(ROOT, "scripts", script * ".jl")) $args`
    err = tempname()
    proc = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = err))
    ok = success(proc)
    msg = ""
    if !ok
        # Solo el tipo de la excepción (p. ej. «ERROR: BoundsError»); no se guarda el resto.
        linea = first(filter(l -> startswith(l, "ERROR"), [readlines(err); "ERROR: fallo sin detalle"]))
        msg = linea[1:min(end, 160)]
    end
    rm(err; force = true)
    return ok, msg
end

function anotar(id, paso, ok, seg, msg)
    dir = joinpath(INTERIM, id)
    isdir(dir) || return
    f = joinpath(dir, "lote_estado.csv")
    fila = DataFrame(paso = [paso], estado = [ok ? "ok" : "error"], fecha = [string(now())],
                     segundos = [round(seg; digits = 1)], mensaje = [msg])
    previo = isfile(f) ? CSV.read(f, DataFrame; types = String) : DataFrame()
    isempty(previo) || (previo = previo[previo.paso .!= paso, :])
    CSV.write(f, isempty(previo) ? fila : vcat(previo, string.(fila); cols = :union))
end

ids = sort(id_registro.(raw_files("gt3x"; proyecto)))
isempty(ids) && error("No hay archivos .gt3x en data/raw/")
isempty(pedidos) || (ids = filter(in(pedidos), ids))
isempty(ids) && error("Ningún ID pedido está en data/raw/")

@printf("%-14s %-8s %s\n", "ID", "paso", "estado")
fallos = String[]
for id in ids
    if leido(id) && !rehacer
        @printf("%-14s %-8s %s\n", id, "01", "ya leído")
        continue
    end
    # Un directorio a medias (sin diagnostics.arrow) o --rehacer: 01 necesita --force para rehacerlo.
    args = isdir(joinpath(INTERIM, id)) ? [id, "--force"] : [id]
    if simular
        @printf("%-14s %-8s %s\n", id, "01", "se leería (simulación)")
        continue
    end
    t = @elapsed ((ok, msg) = ejecutar("01_read_gt3x", args))
    ok = ok && leido(id)
    ok || isempty(msg) && (msg = "no se escribió diagnostics.arrow")
    anotar(id, "01", ok, t, msg)
    @printf("%-14s %-8s %s (%.0f s)\n", id, "01", ok ? "ok" : "ERROR", t)
    ok || push!(fallos, id)
end

if analisis && !simular
    isempty(fallos) || @warn "Hay registros con error en la lectura; el análisis los omite: $(join(fallos, ", "))"
    leidos = filter(leido, ids)
    isempty(leidos) && error("Ningún registro leído: no hay nada que analizar")
    mkpath(SALIDA)
    tanda = DataFrame(paso = String[], estado = String[], segundos = Float64[], mensaje = String[])
    # Todos los IDs leídos en data/interim (no solo los pedidos): las tablas se reescriben enteras.
    for paso in PASOS_ANALISIS
        t = @elapsed ((ok, msg) = ejecutar(paso))
        push!(tanda, (paso, ok ? "ok" : "error", round(t; digits = 1), msg))
        @printf("%-14s %-8s %s (%.0f s)\n", "(todos)", split(paso, '_')[1], ok ? "ok" : "ERROR", t)
        ok || break
        CSV.write(joinpath(SALIDA, "lote_pasos.csv"), tanda)
    end
    CSV.write(joinpath(SALIDA, "lote_pasos.csv"), tanda)
elseif analisis
    println("(simulación) se lanzarían: ", join(PASOS_ANALISIS, ", "))
end

isempty(fallos) || (println("Con error: ", join(fallos, ", ")); exit(1))
