# Lectura de info.txt con lista blanca de campos.
# Especificación: GT3X-File-Format/info.txt.md y NHANES-GT3X-File-Format/fileformats/info.txt.md.

"""
Campos de `info.txt` que se leen. Cualquier otra clave (p. ej. `Subject Name`,
`DateOfBirth`, `Age`, `Sex`, `Height`, `Mass`, `Limb`, `Side`, `Dominance`) se descarta
sin almacenarse.
"""
const INFO_WHITELIST = (
    "Serial Number", "Device Type", "Firmware", "Battery Voltage", "Sample Rate",
    "Start Date", "Stop Date", "Last Sample Time", "Download Date", "TimeZone",
    "Board Revision", "Unexpected Resets",
    "Acceleration Scale", "Acceleration Min", "Acceleration Max",
)

"""
    GT3XInfo

Metadatos del dispositivo leídos de `info.txt` (solo campos de [`INFO_WHITELIST`](@ref)).
Las fechas están en hora local del dispositivo, tal como las guarda ActiLife; `tz_offset`
es el desfase respecto a UTC declarado en `TimeZone` (`nothing` si no consta).
"""
Base.@kwdef struct GT3XInfo
    serial::String
    device_type::String = ""
    firmware::String = ""
    battery_voltage::Union{Float64,Nothing} = nothing
    sample_rate::Int
    start::DateTime
    stop::Union{DateTime,Nothing} = nothing
    last_sample::Union{DateTime,Nothing} = nothing
    download::Union{DateTime,Nothing} = nothing
    tz_offset::Union{Minute,Nothing} = nothing
    board_revision::Union{Int,Nothing} = nothing
    unexpected_resets::Union{Int,Nothing} = nothing
    accel_scale::Union{Float64,Nothing} = nothing
    accel_min::Union{Float64,Nothing} = nothing
    accel_max::Union{Float64,Nothing} = nothing
end

const DOTNET_EPOCH = DateTime(1, 1, 1)

"""
    ticks2datetime(ticks) -> Union{DateTime,Nothing}

Convierte *ticks* de .NET (intervalos de 100 ns desde 0001-01-01T00:00:00) a `DateTime`
con resolución de milisegundos (se trunca lo inferior al ms, p. ej. en `Download Date`).
`0` significa «sin fecha» y devuelve `nothing`.
"""
function ticks2datetime(ticks::Integer)
    ticks == 0 && return nothing
    return DOTNET_EPOCH + Millisecond(div(ticks, 10_000))
end

# Algunos firmwares escriben decimales con coma ("4,03").
_parsefloat(s) = parse(Float64, replace(s, ',' => '.'))

# "02:00:00", "-05:00:00" o "05:30" → Minute
function _parseoffset(s)
    neg = startswith(s, '-')
    parts = parse.(Int, split(lstrip(s, ['-', '+']), ':'))
    m = 60 * parts[1] + (length(parts) > 1 ? parts[2] : 0)
    return Minute(neg ? -m : m)
end

"""
    parse_info(text::AbstractString) -> GT3XInfo

Interpreta el contenido de `info.txt` (formato `Clave: valor` por línea). Solo se
conservan las claves de [`INFO_WHITELIST`](@ref).
"""
function parse_info(text::AbstractString)
    kv = Dict{String,String}()
    for line in eachline(IOBuffer(text))
        i = findfirst(':', line)
        i === nothing && continue
        key = strip(line[1:prevind(line, i)])
        key in INFO_WHITELIST || continue          # lista blanca: lo demás ni se guarda
        kv[key] = strip(line[nextind(line, i):end])
    end
    for k in ("Serial Number", "Sample Rate", "Start Date")
        haskey(kv, k) || error("info.txt sin el campo obligatorio \"$k\"")
    end
    get_or(k, f) = haskey(kv, k) && !isempty(kv[k]) ? f(kv[k]) : nothing
    dt(k) = get_or(k, s -> ticks2datetime(parse(Int64, s)))
    return GT3XInfo(;
        serial = kv["Serial Number"],
        device_type = get(kv, "Device Type", ""),
        firmware = get(kv, "Firmware", ""),
        battery_voltage = get_or("Battery Voltage", _parsefloat),
        sample_rate = round(Int, _parsefloat(kv["Sample Rate"])),
        start = dt("Start Date"),
        stop = dt("Stop Date"),
        last_sample = dt("Last Sample Time"),
        download = dt("Download Date"),
        tz_offset = get_or("TimeZone", _parseoffset),
        board_revision = get_or("Board Revision", s -> parse(Int, s)),
        unexpected_resets = get_or("Unexpected Resets", s -> parse(Int, s)),
        accel_scale = get_or("Acceleration Scale", _parsefloat),
        accel_min = get_or("Acceleration Min", _parsefloat),
        accel_max = get_or("Acceleration Max", _parsefloat),
    )
end

_zip(path::AbstractString) = ZipReader(read(path))

"""
    read_info(path) -> GT3XInfo

Lee `info.txt` de un archivo `.gt3x` (lista blanca de campos).
"""
read_info(path::AbstractString) = parse_info(zip_readentry(_zip(path), "info.txt", String))

"""
    default_accel_scale(serial) -> Float64

Factor de escala (LSB/g) por prefijo de número de serie, según la especificación:
NEO/CLE → 341 (±6 g); MOS → 256 (±8 g). Solo se usa si no hay registro PARAMETERS
ni `Acceleration Scale` en `info.txt`.
"""
function default_accel_scale(serial::AbstractString)
    startswith(serial, "MOS") && return 256.0
    (startswith(serial, "NEO") || startswith(serial, "CLE") || startswith(serial, "MRA")) && return 341.0
    error("No se conoce el factor de escala para el número de serie con prefijo \"$(first(serial, 3))\"")
end

function Base.show(io::IO, ::MIME"text/plain", i::GT3XInfo)
    println(io, "GT3XInfo: ", i.device_type, " ", i.serial, " · firmware ", i.firmware, " · ", i.sample_rate, " Hz")
    print(io, "  ", i.start, " → ", something(i.stop, "—"), " · TZ ", something(i.tz_offset, "—"))
end
