"""
    GT3X

Lectura nativa de archivos ActiGraph `.gt3x` (ZIP con `info.txt` y `log.bin`, o el
formato antiguo NHANES `activity.bin` + `lux.bin`), siguiendo las especificaciones
públicas de ActiGraph:

- <https://github.com/actigraph/GT3X-File-Format> (log.bin)
- <https://github.com/actigraph/NHANES-GT3X-File-Format> (activity.bin + lux.bin)

Privacidad: los campos de `info.txt` se leen con **lista blanca** ([`INFO_WHITELIST`](@ref));
nunca se conservan datos del sujeto (nombre, fecha de nacimiento, edad, sexo, talla,
peso…). El registro METADATA (0x06) de `log.bin`, que contiene datos del sujeto en
JSON, se salta sin decodificarlo.
"""
module GT3X

using Dates
using DataFrames
using Arrow
using ZipArchives: ZipReader, zip_names, zip_readentry

export GT3XInfo, GT3XRecording, read_info, parse_info, read_gt3x, write_interim, read_interim, impute_idle_sleep

include("GT3X/info.jl")
include("GT3X/logbin.jl")
include("GT3X/legacy.jl")
include("GT3X/recording.jl")
include("GT3X/interim.jl")
include("GT3X/impute.jl")

end # module
