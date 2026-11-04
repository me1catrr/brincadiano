# 25 · Inventario trazable de los archivos originales de data/raw/ (y manifiesto de migración).
#
# Qué hace
#   Recorre data/raw/ entero (subcarpetas incluidas) y escribe una fila por archivo: proyecto,
#   sujeto, visita, formato, frecuencia, epoch, fechas, colocación, procedencia, tamaño, SHA-256
#   y la carpeta que le corresponde en la estructura por proyecto (`destino`, docs/datos_raw.md).
#   La atribución sale de Brincadiano.Rutas.clasificar (el nombre y lo declarado por el grupo, con
#   su estado: «declarado (fuente)» o «hipótesis»). Los metadatos salen de la lista blanca de cada
#   formato (info.txt de los .gt3x con GT3X.read_info; settings de los .agd con AGD.read_agd_info);
#   de los RAW.csv, los config y los Excel solo se miran el nombre, el tamaño y el hash. Ningún
#   archivo se mueve ni se modifica.
#   Con `--verificar <manifiesto.csv>` comprueba que cada archivo del manifiesto sigue existiendo
#   con el mismo nombre y el mismo SHA-256, esté donde esté: sirve para validar una migración.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/25_inventario_raw.jl
#   julia --project=scripts scripts/25_inventario_raw.jl --verificar outputs/inventario/raw_inventario_pre.csv
#
# Salida (no se versiona)
#   outputs/inventario/raw_inventario.csv
#
# Consola: recuentos por proyecto y formato; en --verificar, los nombres que no cuadran.

using Brincadiano.GT3X, Brincadiano.AGD
using Brincadiano.Rutas: Rutas, clasificar
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()
const SAL = joinpath(Rutas.outputs_dir(), "inventario")

sha256(path) = first(split(read(`shasum -a 256 $path`, String)))

function inventario()
    filas = NamedTuple[]
    for (dir, dirs, fs) in walkdir(RAW)
        filter!(d -> !startswith(d, "."), dirs)
        for f in sort(fs)
            startswith(f, ".") && continue
            p = joinpath(dir, f); a = clasificar(f)
            fs_hz = missing; epoch_s = missing; inicio = missing; fin = missing; descarga = missing
            serie = missing; modelo = missing; firmware = missing; meta = "no se abre"
            try
                if a.formato == "gt3x"
                    i = read_info(p)
                    fs_hz = i.sample_rate; inicio = i.start; fin = something(i.stop, missing)
                    descarga = something(i.download, missing); serie = i.serial; modelo = i.device_type
                    firmware = i.firmware; meta = "info.txt (lista blanca)"
                elseif a.formato == "agd"
                    i = read_agd_info(p)
                    fs_hz = something(i.sample_rate, missing); epoch_s = i.epoch_s; inicio = i.start
                    fin = something(i.stop, missing); descarga = something(i.download, missing)
                    serie = i.serial; modelo = i.device; firmware = i.firmware; meta = "settings (lista blanca)"
                end
            catch err
                meta = "error al leer: " * first(sprint(showerror, err), 80)
            end
            push!(filas, (; archivo = f, ruta = relpath(p, ROOT), a.formato, a.proyecto, estado_proyecto = a.estado,
                a.sujeto, a.visita, frecuencia_hz = fs_hz, epoch_s, inicio, fin, descarga, serie, modelo, firmware,
                a.colocacion, a.procedencia, metadatos = meta, bytes = filesize(p), sha256 = sha256(p),
                destino = relpath(joinpath(RAW, a.destino), ROOT), en_destino = relpath(p, RAW) == a.destino))
        end
    end
    return DataFrame(filas)
end

function verificar(manifiesto)
    man = CSV.read(manifiesto, DataFrame; types = String)
    ahora = inventario()
    malos = 0
    for r in eachrow(man)
        k = findall(==(r.sha256), ahora.sha256)
        if isempty(k)
            println("  FALTA o CAMBIÓ: ", r.archivo); malos += 1
        elseif !any(ahora.archivo[k] .== r.archivo)
            println("  RENOMBRADO: ", r.archivo); malos += 1
        end
    end
    @printf("Manifiesto: %d archivos; %d con problemas; %d archivos ahora; %d en su carpeta de destino\n",
            nrow(man), malos, nrow(ahora), count(ahora.en_destino))
    return malos
end

i = findfirst(==("--verificar"), ARGS)
i === nothing || exit(verificar(ARGS[i + 1]) == 0 ? 0 : 1)

df = inventario()
df.copia_identica_de = [(k = findfirst(==(h), df.sha256); k < j ? df.ruta[k] : missing) for (j, h) in enumerate(df.sha256)]
mkpath(SAL)
CSV.write(joinpath(SAL, "raw_inventario.csv"), df)
for g in groupby(sort(df, [:proyecto, :formato]), [:proyecto, :formato])
    @printf("  %-12s %-17s %3d archivos  (%d en su carpeta; %d copias idénticas)\n", g.proyecto[1], g.formato[1],
            nrow(g), count(g.en_destino), count(!ismissing, g.copia_identica_de))
end
println("Total: ", nrow(df), " archivos. Salida en ", relpath(SAL, ROOT))
