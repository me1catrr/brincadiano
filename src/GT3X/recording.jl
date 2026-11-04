# Ensamblado de un registro completo a partir de info.txt + log.bin (o activity.bin + lux.bin).

"""
    GT3XRecording

Registro leído de un `.gt3x`. Todos los tiempos son `DateTime` en hora local del
dispositivo (el desfase UTC está en `info.tz_offset`).

- `id`: identificador pseudonimizado (nombre del archivo sin extensión).
- `format`: `:logbin` (formato actual) o `:nhanes` (activity.bin + lux.bin).
- `scale`: factor de escala usado (LSB/g) y `scale_source` su procedencia.
- `acc`: `time`, `x`, `y`, `z` (`Float32`, g). Solo muestras realmente registradas.
- `gaps`: huecos **marcados, no rellenados**: `start`, `stop` (excluido), `n_missing`,
  `kind` ∈ {`"idle_sleep"`, `"missing"`, `"partial"`, `"zero_padding"`} (este último solo
  en `:nhanes`: cola de ceros de relleno de `activity.bin`).
- `lux`: en `:logbin`, registro de cambios (`time`, `lux_raw`): el dispositivo solo
  escribe la luz cuando cambia. Conversión a lux sin documentar para wGT3X-BT
  (no resuelto). En `:nhanes`, una fila por segundo con `lux_raw` y `lux` convertido.
- `battery` (`time`, `mv`), `capsense` (`time`, `signal`, `reference`, `worn`, `bursts`).
- `diagnostics`: recuentos de registros, errores de *checksum*, eventos USB, etc.
"""
struct GT3XRecording
    id::String
    info::GT3XInfo
    format::Symbol
    scale::Float64
    scale_source::String
    acc::DataFrame
    gaps::DataFrame
    lux::DataFrame
    battery::DataFrame
    capsense::DataFrame
    diagnostics::Dict{String,Any}
end

"""
    read_gt3x(path; id = nombre del archivo) -> GT3XRecording

Lee un archivo `.gt3x` completo. Detecta el formato por el contenido del ZIP:
`log.bin` (actual) o `activity.bin` (NHANES, GT3X+/ActiSleep+ con firmware ≤ 2.5.0).
"""
function read_gt3x(path::AbstractString; id::AbstractString = splitext(basename(path))[1])
    zr = _zip(path)
    names = zip_names(zr)
    info = parse_info(zip_readentry(zr, "info.txt", String))
    if "log.bin" in names
        return _read_logbin(zr, info, String(id))
    elseif "activity.bin" in names
        return _read_nhanes(zr, info, String(id))
    else
        error("$(id): el .gt3x no contiene log.bin ni activity.bin")
    end
end

# Instantes de cada muestra: segundo + k/fs, con resolución de ms. Se usa la misma
# convención que ActiLife (techo: .034, .067, .134… a 30 Hz; comprobado frente a sus
# exportaciones RAW.csv) para que ambas series casen por tiempo exacto. El instante
# verdadero de la muestra k es k/fs; la discrepancia es < 1 ms.
function _sample_times(seconds::AbstractVector{<:Integer}, nsamples::AbstractVector{Int}, fs::Int)
    offs = [Millisecond(cld(1000k, fs)) for k in 0:(fs - 1)]
    t = Vector{DateTime}(undef, sum(nsamples))
    j = 1
    @inbounds for (s, ns) in zip(seconds, nsamples)
        t0 = unix2datetime(s)
        for k in 1:ns
            t[j] = t0 + offs[k]; j += 1
        end
    end
    return t
end

function _acc_frame(times, samples::Vector{NTuple{3,Int16}}, scale::Float64)
    # k/256 es exacto en Float32; para 341 LSB/g el error relativo es < 6e-8.
    return DataFrame(time = times,
                     x = Float32[s[1] / scale for s in samples],
                     y = Float32[s[2] / scale for s in samples],
                     z = Float32[s[3] / scale for s in samples])
end

function _read_logbin(zr, info::GT3XInfo, id::String)
    fs = info.sample_rate
    lb = parse_logbin(zip_readentry(zr, "log.bin"); sample_rate = fs)
    diag = Dict{String,Any}(lb.diagnostics)

    # Escala: PARAMETERS > info.txt > prefijo del número de serie.
    scale, source = if haskey(lb.params, PARAM_ACCEL_SCALE)
        ssp_decode(lb.params[PARAM_ACCEL_SCALE]), "PARAMETERS"
    elseif info.accel_scale !== nothing
        info.accel_scale, "info.txt"
    else
        default_accel_scale(info.serial), "número de serie"
    end
    if haskey(lb.params, PARAM_SAMPLE_RATE) && Int(lb.params[PARAM_SAMPLE_RATE]) != fs
        @warn "$(id): Sample Rate de info.txt ($fs) ≠ PARAMETERS ($(Int(lb.params[PARAM_SAMPLE_RATE])))"
        diag["sample_rate_mismatch"] = Int(lb.params[PARAM_SAMPLE_RATE])
    end

    samples = lb.samples
    # wGT3X-BT con firmware 1.6.0 rota los ejes 90° sobre Z (LogRecords/Activity.md).
    if startswith(info.serial, "MOS") && info.firmware == "1.6.0"
        samples = [(s[2], -s[1], s[3]) for s in samples]
        diag["firmware_1_6_0_axis_fix"] = true
        @warn "$(id): firmware 1.6.0 de wGT3X-BT; se corrige la rotación de ejes X/Y"
    end

    times = _sample_times(lb.act_seconds, lb.act_nsamples, fs)
    acc = _acc_frame(times, samples, scale)

    first_s = round(Int, datetime2unix(info.start))
    last_t = something(info.last_sample, info.stop,
                       isempty(lb.act_seconds) ? info.start : unix2datetime(lb.act_seconds[end] + 1))
    gaps = activity_gaps(lb, fs; first_second = first_s, end_second = round(Int, datetime2unix(last_t)))

    lux = DataFrame(time = unix2datetime.(lb.lux_seconds), lux_raw = lb.lux_raw)
    battery = DataFrame(time = unix2datetime.(lb.battery_seconds), mv = lb.battery_mv)
    capsense = DataFrame(time = unix2datetime.(lb.capsense_seconds),
                         signal = getindex.(lb.capsense, 1), reference = getindex.(lb.capsense, 2),
                         worn = getindex.(lb.capsense, 3) .== 1, bursts = getindex.(lb.capsense, 4))

    diag["record_counts"] = Dict(string("0x", string(k; base = 16, pad = 2)) => v for (k, v) in lb.record_counts)
    diag["idle_sleep_enter_events"] = length(lb.ism_enter)
    diag["idle_sleep_exit_events"] = length(lb.ism_exit)
    haskey(lb.params, PARAM_UTC_OFFSET) &&
        (diag["parameters_utc_offset_s"] = Int(reinterpret(Int32, lb.params[PARAM_UTC_OFFSET])))

    return GT3XRecording(id, info, :logbin, scale, source, acc, gaps, lux, battery, capsense, diag)
end

"""
    expected_samples(rec) -> Int

Muestras esperadas entre `info.start` y el final del registro a la frecuencia nominal.
"""
function expected_samples(rec::GT3XRecording)
    fs = rec.info.sample_rate
    rec.format == :nhanes && return cld(rec.diagnostics["samples"], fs) * fs   # flujo continuo sin fecha de fin
    stop = something(rec.info.last_sample, rec.info.stop)
    return round(Int, (stop - rec.info.start) / Millisecond(1000)) * fs
end

function Base.show(io::IO, ::MIME"text/plain", r::GT3XRecording)
    fs = r.info.sample_rate
    println(io, "GT3XRecording ", r.id, " (", r.format, ", ", fs, " Hz, ", r.scale, " LSB/g de ", r.scale_source, ")")
    println(io, "  ", nrow(r.acc), " muestras registradas · ", nrow(r.gaps), " huecos (",
            count(==("idle_sleep"), r.gaps.kind), " idle sleep) · ", nrow(r.lux), " registros de lux")
    isempty(r.acc) || print(io, "  ", first(r.acc.time), " → ", last(r.acc.time))
end
