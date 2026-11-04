# Exportación/lectura de un GT3XRecording a data/interim/<ID>/ en Arrow (compresión zstd).
#
#   acc.arrow        time, x, y, z            muestras registradas (g)
#   gaps.arrow       start, stop, n_missing, kind
#   lux.arrow        time, lux_raw [, lux]
#   battery.arrow    time, mv
#   capsense.arrow   time, signal, reference, worn, bursts
#   meta.arrow       una fila: metadatos del dispositivo (lista blanca) y del lector
#   diagnostics.arrow  key, value (texto)

const INTERIM_TABLES = ("acc", "gaps", "lux", "battery", "capsense")
const READER_VERSION = v"0.1.0"

function _meta_frame(r::GT3XRecording)
    i = r.info
    opt(x) = something(x, missing)
    return DataFrame(
        id = [r.id], format = [String(r.format)], scale = [r.scale], scale_source = [r.scale_source],
        serial = [i.serial], device_type = [i.device_type], firmware = [i.firmware],
        battery_voltage = [opt(i.battery_voltage)], sample_rate = [i.sample_rate],
        start = [i.start], stop = [opt(i.stop)], last_sample = [opt(i.last_sample)],
        download = [opt(i.download)],
        tz_offset_min = [i.tz_offset === nothing ? missing : Dates.value(i.tz_offset)],
        board_revision = [opt(i.board_revision)], unexpected_resets = [opt(i.unexpected_resets)],
        accel_scale_info = [opt(i.accel_scale)], accel_min = [opt(i.accel_min)], accel_max = [opt(i.accel_max)],
        reader_version = [string(READER_VERSION)],
    )
end

_flatten(d::AbstractDict, prefix = "") =
    reduce(vcat, (v isa AbstractDict ? _flatten(v, string(prefix, k, ".")) : [(string(prefix, k), string(v))]
                  for (k, v) in d); init = Tuple{String,String}[])

"""
    write_interim(rec, dir) -> dir

Escribe el registro en `dir` (se crea si no existe) como tablas Arrow comprimidas con
zstd. Pensado para `data/interim/<ID>/`.
"""
function write_interim(r::GT3XRecording, dir::AbstractString)
    mkpath(dir)
    for name in INTERIM_TABLES
        Arrow.write(joinpath(dir, name * ".arrow"), getfield(r, Symbol(name)); compress = :zstd)
    end
    Arrow.write(joinpath(dir, "meta.arrow"), _meta_frame(r))
    kv = sort(_flatten(r.diagnostics))
    Arrow.write(joinpath(dir, "diagnostics.arrow"), DataFrame(key = first.(kv), value = last.(kv)))
    return dir
end

"""
    read_interim(dir) -> GT3XRecording

Reconstruye un `GT3XRecording` desde `data/interim/<ID>/`. Las tablas se materializan
en memoria (`copycols`), independientes del archivo.
"""
function read_interim(dir::AbstractString)
    tbl(name) = DataFrame(Arrow.Table(joinpath(dir, name * ".arrow")); copycols = true)
    m = tbl("meta")[1, :]
    nm(x) = ismissing(x) ? nothing : x
    info = GT3XInfo(; serial = m.serial, device_type = m.device_type, firmware = m.firmware,
        battery_voltage = nm(m.battery_voltage), sample_rate = m.sample_rate, start = m.start,
        stop = nm(m.stop), last_sample = nm(m.last_sample), download = nm(m.download),
        tz_offset = ismissing(m.tz_offset_min) ? nothing : Minute(m.tz_offset_min),
        board_revision = nm(m.board_revision), unexpected_resets = nm(m.unexpected_resets),
        accel_scale = nm(m.accel_scale_info), accel_min = nm(m.accel_min), accel_max = nm(m.accel_max))
    d = tbl("diagnostics")
    diag = Dict{String,Any}(k => something(tryparse(Int, v), tryparse(Bool, v), v) for (k, v) in zip(d.key, d.value))
    t = Dict(name => tbl(name) for name in INTERIM_TABLES)
    return GT3XRecording(m.id, info, Symbol(m.format), m.scale, m.scale_source,
                         t["acc"], t["gaps"], t["lux"], t["battery"], t["capsense"], diag)
end
