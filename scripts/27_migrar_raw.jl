# 27 · Migración reversible de data/raw/ a la estructura por proyecto (docs/datos_raw.md).
#
# Qué hace
#   Lee el manifiesto del 25 (`outputs/inventario/raw_inventario.csv`: ruta actual, destino y
#   SHA-256 de cada original) y mueve cada archivo que no esté ya en su destino. Por defecto
#   solo dice lo que haría. Antes de mover un archivo comprueba que su SHA-256 es el del
#   manifiesto y que el destino no existe; nunca sobrescribe ni borra. `mv` dentro del mismo
#   volumen renombra: conserva contenido, fechas, permisos (a-w) y atributos. Si el destino está
#   en otro volumen se niega (sería copiar y borrar). Cada movimiento se anota en
#   `outputs/inventario/migracion_<fecha>.csv` (origen, destino, sha256) en cuanto se hace, de
#   modo que `--revertir` puede deshacer también una migración interrumpida.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/25_inventario_raw.jl            # manifiesto previo
#   cp outputs/inventario/raw_inventario.csv outputs/inventario/raw_inventario_pre.csv
#   julia --project=scripts scripts/27_migrar_raw.jl                # simulación
#   julia --project=scripts scripts/27_migrar_raw.jl --ejecutar
#   julia --project=scripts scripts/25_inventario_raw.jl --verificar outputs/inventario/raw_inventario_pre.csv
#   julia --project=scripts scripts/27_migrar_raw.jl --revertir outputs/inventario/migracion_<fecha>.csv
#
# Consola: nombre de archivo (ID pseudonimizado), origen y destino relativos.

using Brincadiano.Rutas: Rutas
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SAL = joinpath(Rutas.outputs_dir(), "inventario")

sha256(path) = first(split(read(`shasum -a 256 $path`, String)))
mismo_volumen(a, b) = stat(a).device == stat(b).device

"Mueve `origen` → `destino` sin sobrescribir; devuelve true si lo ha movido."
function mover(origen, destino, hash; ejecutar)
    isfile(origen) || (println("  NO ESTÁ: ", basename(origen)); return false)
    ispath(destino) && (println("  DESTINO OCUPADO (no se toca): ", relpath(destino, ROOT)); return false)
    sha256(origen) == hash || (println("  SHA-256 DISTINTO del manifiesto (no se mueve): ", basename(origen)); return false)
    @printf("  %-34s %s → %s\n", basename(origen), relpath(dirname(origen), ROOT), relpath(dirname(destino), ROOT))
    ejecutar || return false
    mkpath(dirname(destino))
    mismo_volumen(origen, dirname(destino)) || error("$(basename(origen)): el destino está en otro volumen; no se mueve")
    mv(origen, destino)
    sha256(destino) == hash || error("$(basename(destino)): SHA-256 distinto tras mover; revisar a mano")
    return true
end

function anotar(log, origen, destino, hash)
    nuevo = !isfile(log)
    open(log, "a") do io
        nuevo && println(io, "origen,destino,sha256")
        println(io, "\"", relpath(origen, ROOT), "\",\"", relpath(destino, ROOT), "\",", hash)
    end
end

i = findfirst(==("--revertir"), ARGS)
if i !== nothing
    log = CSV.read(ARGS[i + 1], DataFrame; types = String)
    ejecutar = "--ejecutar" in ARGS
    println(ejecutar ? "Revirtiendo:" : "Simulación de la reversión (añade --ejecutar):")
    n = sum((mover(joinpath(ROOT, r.destino), joinpath(ROOT, r.origen), r.sha256; ejecutar) for r in reverse(eachrow(log))); init = 0)
    println(ejecutar ? "Devueltos $n archivos." : "")
    exit(0)
end

ejecutar = "--ejecutar" in ARGS
man_f = joinpath(SAL, "raw_inventario.csv")
isfile(man_f) || error("Falta el manifiesto: ejecutar antes scripts/25_inventario_raw.jl")
man = CSV.read(man_f, DataFrame; types = String)
pend = man[man.en_destino .!= "true", :]
println(ejecutar ? "Migrando $(nrow(pend)) archivos:" : "Simulación ($(nrow(pend)) archivos; añade --ejecutar para mover):")
function migrar(pend; ejecutar)
    log = joinpath(SAL, "migracion_" * Dates.format(now(), "yyyy-mm-ddTHHMMSS") * ".csv")
    n = 0
    for r in eachrow(pend)
        o = joinpath(ROOT, r.ruta); d = joinpath(ROOT, r.destino)
        if mover(o, d, r.sha256; ejecutar)
            anotar(log, o, d, r.sha256); n += 1
        end
    end
    ejecutar && println("Movidos $n de $(nrow(pend)). Registro: ", relpath(log, ROOT))
end
migrar(pend; ejecutar)
