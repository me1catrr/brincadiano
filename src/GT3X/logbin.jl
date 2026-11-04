# Parser de log.bin (formato actual de .gt3x).
# Especificación: GT3X-File-Format/README.md y LogRecords/*.md.
#
# Registro: 0x1E | tipo (1 B) | timestamp u32 LE (s) | tamaño u16 LE | payload | checksum.
# El checksum es el complemento a uno del XOR de cabecera y payload (valor inicial 0).
# Los timestamps están en hora local del dispositivo (coinciden con Start Date de
# info.txt, comprobado en los archivos del estudio).

const RECORD_SEPARATOR = 0x1E
const REC_ACTIVITY   = 0x00
const REC_BATTERY    = 0x02
const REC_EVENT      = 0x03
const REC_LUX        = 0x05
const REC_METADATA   = 0x06   # datos del sujeto en JSON: nunca se decodifica
const REC_CAPSENSE   = 0x0D
const REC_PARAMETERS = 0x15
const REC_ACTIVITY2  = 0x1A

const EVENT_IDLE_SLEEP_ENTER = 0x08
const EVENT_IDLE_SLEEP_EXIT  = 0x09

# Claves de PARAMETERS (espacio de direcciones, identificador) que se usan.
const PARAM_ACCEL_SCALE = (0x0000, 0x0037)   # float SSP, LSB/g
const PARAM_SAMPLE_RATE = (0x0001, 0x000A)   # entero, Hz
const PARAM_UTC_OFFSET  = (0x0001, 0x0028)   # entero con signo, s

@inline _u16(b, i) = UInt16(b[i]) | UInt16(b[i+1]) << 8
@inline _u32(b, i) = UInt32(b[i]) | UInt32(b[i+1]) << 8 | UInt32(b[i+2]) << 16 | UInt32(b[i+3]) << 24
@inline _i16(b, i) = reinterpret(Int16, _u16(b, i))

"""
    ssp_decode(v::UInt32) -> Float64

Decodifica un valor en coma flotante del registro PARAMETERS: fracción de 24 bits en
complemento a dos (escalada por 2²³) y exponente de 8 bits en complemento a dos
(LogRecords/Parameters.md, codificación «Simple Serial Protocol»).
"""
function ssp_decode(v::UInt32)
    v == 0x007FFFFF && return floatmax(Float64)
    v == 0x00800000 && return -floatmax(Float64)
    e = reinterpret(Int8, UInt8(v >> 24))
    m = Int32(v & 0x00FFFFFF)
    m & 0x00800000 != 0 && (m |= reinterpret(Int32, 0xFF000000))
    return (m / 8_388_608.0) * 2.0^e
end

"""
    unpack12!(dest, b, offset, nsamples)

Desempaqueta `nsamples` muestras triaxiales de 12 bits (orden **YXZ**, bits de más a
menos significativo, complemento a dos) a partir del byte `offset` de `b`, y añade a
`dest` las tuplas `(x, y, z)` como `Int16`. Común a ACTIVITY (log.bin) y activity.bin.
"""
function unpack12!(dest::Vector{NTuple{3,Int16}}, b::AbstractVector{UInt8}, offset::Int, nsamples::Int)
    k = 0                                    # índice del valor de 12 bits (0-based)
    for _ in 1:nsamples
        y = _get12(b, offset, k); x = _get12(b, offset, k + 1); z = _get12(b, offset, k + 2)
        push!(dest, (x, y, z))
        k += 3
    end
    return dest
end

@inline function _get12(b, offset, k)
    o = offset + (12k) >> 3
    v = iseven(k) ? (UInt16(b[o]) << 4) | (UInt16(b[o+1]) >> 4) :
                    (UInt16(b[o] & 0x0F) << 8) | UInt16(b[o+1])
    return v > 0x07FF ? reinterpret(Int16, v | 0xF000) : Int16(v)
end

"""
    LogBin

Contenido de `log.bin` tras el parseo, en unidades del dispositivo (sin escalar):
segundos de cada registro de actividad, muestras `Int16`, registro de cambios de lux,
batería, sensor capacitivo, eventos de *idle sleep*, parámetros y diagnóstico.
"""
Base.@kwdef struct LogBin
    act_seconds::Vector{UInt32} = UInt32[]            # segundo de cada registro ACTIVITY/ACTIVITY2
    act_nsamples::Vector{Int} = Int[]                 # muestras en ese registro
    samples::Vector{NTuple{3,Int16}} = NTuple{3,Int16}[]
    lux_seconds::Vector{UInt32} = UInt32[]
    lux_raw::Vector{UInt16} = UInt16[]
    battery_seconds::Vector{UInt32} = UInt32[]
    battery_mv::Vector{UInt16} = UInt16[]
    capsense_seconds::Vector{UInt32} = UInt32[]
    capsense::Vector{NTuple{4,Int}} = NTuple{4,Int}[]  # (signal, reference, state, bursts)
    ism_enter::Vector{UInt32} = UInt32[]
    ism_exit::Vector{UInt32} = UInt32[]
    params::Dict{Tuple{UInt16,UInt16},UInt32} = Dict{Tuple{UInt16,UInt16},UInt32}()
    record_counts::Dict{UInt8,Int} = Dict{UInt8,Int}()
    diagnostics::Dict{String,Int} = Dict("checksum_errors" => 0, "usb_events" => 0,
        "padding_bytes" => 0, "truncated_tail_bytes" => 0, "non_monotonic_seconds" => 0,
        "oversized_activity" => 0)
end

"""
    parse_logbin(b::AbstractVector{UInt8}; sample_rate) -> LogBin

Recorre todos los registros de `log.bin`. Los registros con *checksum* incorrecto se
descartan y se cuentan; los bytes de relleno entre registros (ceros, permitidos por la
especificación) se saltan; los registros de actividad de 1 byte (eventos USB) se ignoran.
"""
function parse_logbin(b::AbstractVector{UInt8}; sample_rate::Int)
    lb = LogBin()
    d = lb.diagnostics
    n = length(b)
    sizehint!(lb.samples, n ÷ 5)                     # ≈ 4,5 B por muestra en ACTIVITY
    pos = 1
    last_act = UInt32(0)
    while pos <= n
        if b[pos] != RECORD_SEPARATOR
            d["padding_bytes"] += 1
            pos += 1
            continue
        end
        if pos + 8 > n                                # cabecera incompleta al final
            d["truncated_tail_bytes"] += n - pos + 1
            break
        end
        type = b[pos+1]
        ts = _u32(b, pos + 2)
        sz = Int(_u16(b, pos + 6))
        pl = pos + 8                                  # primer byte del payload
        if pl + sz > n                                # payload o checksum fuera del archivo
            d["truncated_tail_bytes"] += n - pos + 1
            break
        end
        ck = 0x00
        @inbounds for i in pos:(pl + sz - 1)
            ck ⊻= b[i]
        end
        if ~ck != b[pl+sz]
            d["checksum_errors"] += 1
            pos = pl + sz + 1
            continue
        end
        lb.record_counts[type] = get(lb.record_counts, type, 0) + 1

        if type == REC_ACTIVITY || type == REC_ACTIVITY2
            if sz == 1                                # evento de conexión USB, no es actividad
                d["usb_events"] += 1
            else
                ns = type == REC_ACTIVITY ? (8sz) ÷ 36 : sz ÷ 6
                if ns > sample_rate
                    d["oversized_activity"] += 1
                    ns = sample_rate
                end
                ts <= last_act && last_act != 0 && (d["non_monotonic_seconds"] += 1)
                last_act = ts
                push!(lb.act_seconds, ts)
                push!(lb.act_nsamples, ns)
                if type == REC_ACTIVITY
                    unpack12!(lb.samples, b, pl, ns)
                else
                    @inbounds for k in 0:(ns - 1)
                        o = pl + 6k
                        push!(lb.samples, (_i16(b, o), _i16(b, o + 2), _i16(b, o + 4)))
                    end
                end
            end
        elseif type == REC_LUX && sz >= 2
            push!(lb.lux_seconds, ts); push!(lb.lux_raw, _u16(b, pl))
        elseif type == REC_EVENT && sz >= 1
            e = b[pl]
            e == EVENT_IDLE_SLEEP_ENTER && push!(lb.ism_enter, ts)
            e == EVENT_IDLE_SLEEP_EXIT && push!(lb.ism_exit, ts)
        elseif type == REC_BATTERY && sz >= 2
            push!(lb.battery_seconds, ts); push!(lb.battery_mv, _u16(b, pl))
        elseif type == REC_CAPSENSE && sz >= 6
            push!(lb.capsense_seconds, ts)
            push!(lb.capsense, (Int(_u16(b, pl)), Int(_u16(b, pl + 2)), Int(b[pl+4]), Int(b[pl+5])))
        elseif type == REC_PARAMETERS
            for k in 0:(sz ÷ 8 - 1)
                o = pl + 8k
                lb.params[(_u16(b, o), _u16(b, o + 2))] = _u32(b, o + 4)
            end
        end
        # REC_METADATA y el resto de tipos: se saltan sin leer su payload.
        pos = pl + sz + 1
    end
    return lb
end

"""
    activity_gaps(lb, sample_rate; first_second, end_second) -> DataFrame

Huecos en la secuencia de segundos con actividad, **sin rellenarlos**. Un hueco va de
`start` (incluido) a `stop` (excluido), en segundos Unix locales. Se clasifica como
`"idle_sleep"` si contiene un evento de entrada en *idle sleep* (0x08) entre el último
segundo registrado y el siguiente; si no, como `"missing"`. Los segundos con menos
muestras que `sample_rate` se marcan como `"partial"`.
`first_second`/`end_second` (excluido) delimitan el periodo esperado del registro.
"""
function activity_gaps(lb::LogBin, sample_rate::Int; first_second::Integer, end_second::Integer)
    starts = Int[]; stops = Int[]; nmiss = Int[]; kinds = String[]
    enter = sort(Int.(lb.ism_enter))
    has_ism(a, z) = (i = searchsortedfirst(enter, a); i <= length(enter) && enter[i] <= z)
    function addgap(a, z, prev)
        z > a || return
        push!(starts, a); push!(stops, z); push!(nmiss, (z - a) * sample_rate)
        push!(kinds, has_ism(prev, z) ? "idle_sleep" : "missing")
    end
    secs = Int.(lb.act_seconds)
    if isempty(secs)
        addgap(Int(first_second), Int(end_second), Int(first_second))
    else
        addgap(Int(first_second), secs[1], Int(first_second))
        for i in eachindex(secs)
            if lb.act_nsamples[i] < sample_rate
                push!(starts, secs[i]); push!(stops, secs[i] + 1)
                push!(nmiss, sample_rate - lb.act_nsamples[i]); push!(kinds, "partial")
            end
            i < length(secs) && addgap(secs[i] + 1, secs[i+1], secs[i])
        end
        addgap(secs[end] + 1, Int(end_second), secs[end])
    end
    return DataFrame(start = unix2datetime.(starts), stop = unix2datetime.(stops),
                     n_missing = nmiss, kind = kinds)
end
