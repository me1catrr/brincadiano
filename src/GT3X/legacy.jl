# Formato NHANES: activity.bin + lux.bin (GT3X+ «NEO» y ActiSleep+ «MRA» con firmware ≤ 2.5.0).
# Especificación: NHANES-GT3X-File-Format/fileformats/{activity.bin,lux.bin}.md y la
# librería de referencia de ActiGraph (src/GT3X.Parsing.Library/Gt3xFile.cs).
#
# - activity.bin: flujo continuo de muestras de 12 bits (YXZ) desde Start Date, sin
#   timestamps ni huecos explícitos. Si el número de muestras es impar sobra un nibble.
# - lux.bin: un UInt16 LE por segundo desde Start Date.

"""
    nhanes_lux(raw, serial) -> Float64

Conversión de `lux.bin` a lux según la especificación NHANES: lecturas < 20 o 0xFFFF → 0;
el resto × 1,25 (GT3X+) o × 3,25 (ActiSleep+, serie «MRA»), con máximo 2500 o 6000 lux,
redondeado al entero con empates lejos de cero (como la librería de referencia).
"""
function nhanes_lux(raw::UInt16, serial::AbstractString)
    actisleep = startswith(serial, "MRA")
    factor, maxlux = actisleep ? (3.25, 6000.0) : (1.25, 2500.0)
    (raw == 0xFFFF || raw < 20) && return 0.0
    return round(min(raw * factor, maxlux), RoundNearestTiesAway)
end

function _read_nhanes(zr, info::GT3XInfo, id::String)
    fs = info.sample_rate
    b = zip_readentry(zr, "activity.bin")
    ns = (8 * length(b)) ÷ 36
    samples = unpack12!(sizehint!(NTuple{3,Int16}[], ns), b, 1, ns)

    scale, source = info.accel_scale !== nothing ? (info.accel_scale, "info.txt") :
                                                   (default_accel_scale(info.serial), "número de serie")
    s0 = round(Int, datetime2unix(info.start))
    nsec = cld(ns, fs)

    # Cola de relleno: activity.bin puede terminar en ceros hasta completar la página de
    # memoria (p. ej. el ejemplo público de ActiGraph). Una muestra (0, 0, 0) es imposible
    # físicamente (la gravedad siempre está presente), así que una cola de ≥ 1 s de ceros
    # se recorta y se marca como hueco "zero_padding", sin imputarla.
    nkeep = something(findlast(!=((Int16(0), Int16(0), Int16(0))), samples), 0)
    ntrail = ns - nkeep
    if ntrail < fs
        nkeep = ns
    else
        resize!(samples, nkeep)
    end
    nsec_kept = cld(nkeep, fs)
    seconds = s0 .+ (0:(nsec_kept - 1))
    nsamples = [min(fs, nkeep - fs * i) for i in 0:(nsec_kept - 1)]
    acc = _acc_frame(_sample_times(seconds, nsamples, fs), samples, scale)

    # Sin timestamps no hay huecos detectables, salvo el relleno final o un último
    # segundo parcial.
    gaps = DataFrame(start = DateTime[], stop = DateTime[], n_missing = Int[], kind = String[])
    if nkeep < ns
        k = nkeep % fs                         # posición dentro del último segundo conservado
        tstart = unix2datetime(s0 + nkeep ÷ fs) + Millisecond(cld(1000k, fs))
        push!(gaps, (tstart, unix2datetime(s0 + nsec), nsec * fs - nkeep, "zero_padding"))
    elseif nsamples[end] < fs
        push!(gaps, (unix2datetime(seconds[end]), unix2datetime(seconds[end] + 1), fs - nsamples[end], "partial"))
    end

    lux = if "lux.bin" in zip_names(zr)
        lb = zip_readentry(zr, "lux.bin")
        raw = UInt16[_u16(lb, i) for i in 1:2:(length(lb) - 1)]   # todos los pares completos
        DataFrame(time = [unix2datetime(s0 + i) for i in 0:(length(raw) - 1)],
                  lux_raw = raw, lux = nhanes_lux.(raw, Ref(info.serial)))
    else
        DataFrame(time = DateTime[], lux_raw = UInt16[], lux = Float64[])
    end

    diag = Dict{String,Any}("activity_bin_bytes" => length(b), "samples" => ns,
                            "trailing_zero_samples" => ns - nkeep,
                            "lux_seconds" => nrow(lux))
    empty_battery = DataFrame(time = DateTime[], mv = UInt16[])
    empty_capsense = DataFrame(time = DateTime[], signal = Int[], reference = Int[], worn = Bool[], bursts = Int[])
    return GT3XRecording(id, info, :nhanes, scale, source, acc, gaps, lux, empty_battery, empty_capsense, diag)
end
