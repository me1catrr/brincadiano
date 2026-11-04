# 01 · Lectura de los .gt3x originales y exportación a Arrow.
#
# Qué hace
#   Primer paso del pipeline. Lee cada data/raw/*.gt3x con Brincadiano.GT3X (log.bin
#   actual o activity.bin + lux.bin del formato antiguo) y lo deja en Arrow. Los huecos
#   (idle sleep, muestras ausentes, segundo incompleto, relleno de ceros) se guardan en
#   la tabla de huecos: este script no los rellena. Comprueba que las muestras
#   registradas más las que faltan coinciden con las esperadas a la frecuencia del
#   dispositivo, y avisa si se descartó algún registro por checksum.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/01_read_gt3x.jl            # solo los directorios que falten
#   julia --project=scripts scripts/01_read_gt3x.jl --force    # borra la lógica de «ya existe» y rehace todos
#   julia --project=scripts scripts/01_read_gt3x.jl 6005_V00   # solo esos IDs (lo usa el 16 para el lote)
#   julia --project=scripts scripts/01_read_gt3x.jl --proyecto MiopiaUNAV   # solo un proyecto
#   Los .gt3x se buscan en cualquier carpeta de data/raw/ (Brincadiano.Rutas.raw_files).
#
# Depende de
#   Nada del resto del pipeline. Lee data/raw/ en solo lectura.
#
# Lo usan
#   02 (tabla de cobertura), 03 (figuras de la señal) y 04 (control de calidad).
#   A partir de estos Arrow siguen 05, 06, 08, 09 y 10. El 92, el 93, el 94 y el 96
#   vuelven a abrir el mismo interim. El 95 relee estos Arrow para las tres variantes
#   de imputación. El 90 no los usa: vuelve a abrir el .gt3x para compararlo con
#   read.gt3x y pygt3x.
#
# Entrada
#   data/raw/*.gt3x
#
# Salida (no se versiona)
#   data/interim/<ID>/acc.arrow           aceleración triaxial, solo muestras registradas
#   data/interim/<ID>/gaps.arrow          huecos, con clase y número de muestras que faltan
#   data/interim/<ID>/lux.arrow           luz bruta (solo cuando el valor cambia)
#   data/interim/<ID>/battery.arrow       tensión de batería
#   data/interim/<ID>/capsense.arrow      contacto con la piel, un registro por minuto
#   data/interim/<ID>/meta.arrow          metadatos de la lista blanca (sin datos del sujeto)
#   data/interim/<ID>/diagnostics.arrow   formato, escala, errores de checksum
#
# Consola: solo el ID pseudonimizado, el formato, las horas, el porcentaje registrado
# y el recuento de huecos. Ningún campo personal de info.txt.

using Brincadiano.GT3X
using Brincadiano.Rutas: Rutas, raw_files, id_registro
using Brincadiano.GT3X: expected_samples
using DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const INTERIM = joinpath(ROOT, "data", "interim")
force = "--force" in ARGS

i_proy = findfirst(==("--proyecto"), ARGS)
proyecto = i_proy === nothing ? nothing : ARGS[i_proy + 1]
files = raw_files("gt3x"; proyecto)          # en cualquier carpeta de data/raw/
ids_pedidos = filter(a -> !startswith(a, "--") && a != proyecto, ARGS)
isempty(ids_pedidos) || (files = filter(f -> id_registro(f) in ids_pedidos, files))
isempty(files) && error("No hay archivos .gt3x en data/raw/ que coincidan")
@printf("%-14s %8s %9s %10s %8s %7s %s\n", "ID", "formato", "horas", "registrado", "huecos", "ISM", "estado")
for f in files
    id = id_registro(f)
    out = joinpath(INTERIM, id)
    if isdir(out) && !force
        @printf("%-14s %s\n", id, "ya existe (usa --force para rehacer)")
        continue
    end
    r = read_gt3x(f)
    exp = expected_samples(r)
    ok = nrow(r.acc) + sum(r.gaps.n_missing) == exp
    write_interim(r, out)
    @printf("%-14s %8s %9.1f %9.1f%% %8d %7d %s\n", id, r.format, exp / r.info.sample_rate / 3600,
            100 * nrow(r.acc) / exp, nrow(r.gaps), count(==("idle_sleep"), r.gaps.kind),
            ok ? "ok" : "¡contabilidad no cuadra!")
    r.diagnostics["checksum_errors"] isa Int && r.diagnostics["checksum_errors"] > 0 &&
        @warn "$id: $(r.diagnostics["checksum_errors"]) registros con checksum incorrecto descartados"
end
