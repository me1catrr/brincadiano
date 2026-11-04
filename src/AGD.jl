"""
    AGD

Lectura de los `.agd` de ActiLife (base de datos SQLite con los *epochs* ya procesados:
counts de los tres ejes, pasos, lux e inclinómetro). Se usan como **referencia** frente a lo
que el pipeline calcula desde el `.gt3x` original, nunca como entrada del pipeline.

La base se abre con el `sqlite3` del sistema como proceso aparte, en solo lectura
(`mode=ro&immutable=1`), sin dependencias nuevas de Julia.

Privacidad: la tabla `settings` mezcla campos del dispositivo con datos del sujeto
(`subjectname`, `sex`, `height`, `mass`, `age`, `race`, `limb`, `side`, `dominance`,
`dateOfBirth`). La consulta filtra por [`AGD_SETTINGS_WHITELIST`](@ref) **dentro de SQLite**:
los demás valores no llegan a salir de la base. De `logDiaryTimes` (diario de retiradas) y
`proximity` (otros dispositivos) solo se cuentan las filas.
"""
module AGD

using Dates
using DataFrames
using ..GT3X: ticks2datetime

export AGDInfo, AGDRecording, read_agd, read_agd_info, AGD_SETTINGS_WHITELIST

"""
Campos de `settings` que se leen. Cualquier otro (datos del sujeto, `machinename`,
`osversion`, `notes`…) se queda en la base.
"""
const AGD_SETTINGS_WHITELIST = (
    "softwarename", "softwareversion", "finished", "devicename", "deviceserial",
    "deviceversion", "filter", "modenumber", "epochlength", "epochcount", "startdatetime",
    "stopdatetime", "downloaddatetime", "batteryvoltage", "original sample rate",
    "unexpectedResets", "datetimeformat", "decimal", "culturename", "proximityIntervalInSeconds",
)

"Columnas de la tabla `data` que se leen, si existen (los AGD antiguos no traen todas)."
const AGD_DATA_COLUMNS = ("axis1", "axis2", "axis3", "steps", "lux",
                          "inclineOff", "inclineStanding", "inclineSitting", "inclineLying")

"Tablas de las que solo se cuenta el número de filas."
const AGD_COUNTED_TABLES = ("filters", "capsense", "sleep", "awakenings", "logDiaryTimes",
                            "proximity", "crouterEpoch", "crouterMinute")

"""
    AGDInfo

Metadatos de un `.agd` (solo [`AGD_SETTINGS_WHITELIST`](@ref)). Las fechas están en hora
local, como en `info.txt`. `filter` es el filtro de los counts («Normal» o «Low Frequency
Extension»); `mode` es el `modenumber` de ActiLife, que no se decodifica: lo que contiene el
archivo se ve en las columnas de `data`.
"""
Base.@kwdef struct AGDInfo
    serial::String
    device::String = ""
    firmware::String = ""
    software::String = ""
    software_version::String = ""
    filter::String = ""
    mode::Union{Int,Nothing} = nothing
    epoch_s::Int
    epoch_count::Union{Int,Nothing} = nothing
    sample_rate::Union{Int,Nothing} = nothing
    start::DateTime
    stop::Union{DateTime,Nothing} = nothing
    download::Union{DateTime,Nothing} = nothing
    battery_voltage::Union{Float64,Nothing} = nothing
    unexpected_resets::Union{Int,Nothing} = nothing
    finished::Union{Bool,Nothing} = nothing
end

"""
    AGDRecording

Un `.agd` leído: `info`, `epochs` (una fila por *epoch*: `time` y las columnas de
[`AGD_DATA_COLUMNS`](@ref) presentes), `filters` (tramos de tiempo que ActiLife marcó en el
archivo, `start`/`stop`) y `rows` (filas de cada tabla de [`AGD_COUNTED_TABLES`](@ref);
`missing` si la tabla no existe).
"""
struct AGDRecording
    info::AGDInfo
    epochs::DataFrame
    filters::DataFrame
    rows::Dict{String,Union{Int,Missing}}
end

const _SEP = "\x1f"

function _uri(path)
    p = abspath(path)
    p = replace(p, "%" => "%25", " " => "%20", "?" => "%3F", "#" => "%23")
    return "file:" * p * "?mode=ro&immutable=1"
end

"Ejecuta una consulta y devuelve sus filas como vectores de `String`."
function _query(path, sql; sqlite3 = "sqlite3")
    out = read(`$sqlite3 -readonly -batch -noheader -separator $_SEP $(_uri(path)) $sql`, String)
    return [split(l, _SEP) for l in eachline(IOBuffer(out)) if !isempty(l)]
end

_tables(path; kw...) = Set(first.(_query(path, "SELECT name FROM sqlite_master WHERE type = 'table';"; kw...)))
_columns(path, t; kw...) = [r[2] for r in _query(path, "PRAGMA table_info($t);"; kw...)]

_parsefloat(s) = parse(Float64, replace(s, ',' => '.'))
_ticks(s) = ticks2datetime(parse(Int64, s))

"""
    parse_agd_settings(kv::AbstractDict) -> AGDInfo

Interpreta los pares `settingName => settingValue`. Solo se miran las claves de
[`AGD_SETTINGS_WHITELIST`](@ref); el resto, si llegara, se ignora.
"""
function parse_agd_settings(kv::AbstractDict)
    kv = Dict(String(k) => String(v) for (k, v) in kv if k in AGD_SETTINGS_WHITELIST)
    for k in ("deviceserial", "epochlength", "startdatetime")
        haskey(kv, k) || error("AGD sin el campo obligatorio \"$k\"")
    end
    get_or(k, f) = haskey(kv, k) && !isempty(kv[k]) ? f(kv[k]) : nothing
    return AGDInfo(;
        serial = kv["deviceserial"],
        device = get(kv, "devicename", ""),
        firmware = get(kv, "deviceversion", ""),
        software = get(kv, "softwarename", ""),
        software_version = get(kv, "softwareversion", ""),
        filter = get(kv, "filter", ""),
        mode = get_or("modenumber", s -> parse(Int, s)),
        epoch_s = parse(Int, kv["epochlength"]),
        epoch_count = get_or("epochcount", s -> parse(Int, s)),
        sample_rate = get_or("original sample rate", s -> round(Int, _parsefloat(s))),
        start = _ticks(kv["startdatetime"]),
        stop = get_or("stopdatetime", _ticks),
        download = get_or("downloaddatetime", _ticks),
        battery_voltage = get_or("batteryvoltage", _parsefloat),
        unexpected_resets = get_or("unexpectedResets", s -> parse(Int, s)),
        finished = get_or("finished", s -> lowercase(s) == "true"),
    )
end

function _settings(path; kw...)
    lista = join(("'" * k * "'" for k in AGD_SETTINGS_WHITELIST), ", ")
    filas = _query(path, "SELECT settingName, settingValue FROM settings WHERE settingName IN ($lista);"; kw...)
    return Dict(r[1] => (length(r) > 1 ? r[2] : "") for r in filas)
end

"""
    read_agd_info(path; sqlite3 = "sqlite3") -> AGDInfo

Lee solo los metadatos (lista blanca) de un `.agd`.
"""
read_agd_info(path::AbstractString; kw...) = parse_agd_settings(_settings(path; kw...))

"""
    read_agd(path; sqlite3 = "sqlite3") -> AGDRecording

Lee un `.agd` completo (metadatos, *epochs*, filtros y recuentos de tablas). Los valores
ausentes en SQLite (`NULL`) quedan como `missing`.
"""
function read_agd(path::AbstractString; kw...)
    isfile(path) || error("No existe el archivo")
    info = read_agd_info(path; kw...)
    tablas = _tables(path; kw...)
    cols = "data" in tablas ? filter(in(_columns(path, "data"; kw...)), collect(AGD_DATA_COLUMNS)) : String[]
    epochs = DataFrame(time = DateTime[])
    if "data" in tablas
        filas = _query(path, "SELECT " * join(["dataTimestamp"; cols], ", ") * " FROM data ORDER BY dataTimestamp;"; kw...)
        epochs = DataFrame(time = [_ticks(r[1]) for r in filas])
        for (j, c) in enumerate(cols)
            epochs[!, c] = [isempty(r[j + 1]) ? missing : _parsefloat(r[j + 1]) for r in filas]
        end
    end
    filters = DataFrame(start = DateTime[], stop = DateTime[])
    if "filters" in tablas
        for r in _query(path, "SELECT filterStartTimestamp, filterStopTimestamp FROM filters ORDER BY filterStartTimestamp;"; kw...)
            push!(filters, (_ticks(r[1]), _ticks(r[2])))
        end
    end
    rows = Dict{String,Union{Int,Missing}}(t => (t in tablas ?
        parse(Int, only(only(_query(path, "SELECT count(*) FROM $t;"; kw...)))) : missing) for t in AGD_COUNTED_TABLES)
    return AGDRecording(info, epochs, filters, rows)
end

function Base.show(io::IO, ::MIME"text/plain", r::AGDRecording)
    i = r.info
    println(io, "AGDRecording: ", i.device, " ", i.serial, " · firmware ", i.firmware, " · ",
            i.software, " ", i.software_version, " · epoch ", i.epoch_s, " s · filtro ", i.filter)
    print(io, "  ", i.start, " → ", something(i.stop, "—"), " · ", nrow(r.epochs), " epochs · columnas ",
          join(names(r.epochs)[2:end], ", "))
end

end # module
