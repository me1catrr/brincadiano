using Brincadiano.GT3X: parse_logbin, unpack12!, ssp_decode, activity_gaps, read_gt3x,
                        expected_samples, REC_ACTIVITY, REC_ACTIVITY2, REC_EVENT, REC_LUX,
                        REC_METADATA, REC_PARAMETERS
using Dates, DataFrames

hexbytes(s) = parse.(UInt8, split(s); base = 16)

# Construye un registro de log.bin con checksum correcto (o incorrecto si `bad`).
function rec(type, ts, payload::AbstractVector{UInt8}; bad = false)
    payload = collect(payload)
    h = UInt8[0x1E, type, reinterpret(UInt8, [htol(UInt32(ts))])..., reinterpret(UInt8, [htol(UInt16(length(payload)))])...]
    ck = ~reduce(⊻, vcat(h, payload); init = 0x00)
    return vcat(h, payload, bad ? ck ⊻ 0xFF : ck)
end
# Empaqueta muestras (x, y, z) en 12 bits, orden YXZ, como ACTIVITY.
function pack12(samples)
    vals = UInt16[]
    for (x, y, z) in samples, v in (y, x, z)
        push!(vals, reinterpret(UInt16, Int16(v)) & 0x0FFF)
    end
    bits = BitVector(vcat([reverse(digits(Bool, v; base = 2, pad = 12)) for v in vals]...))
    append!(bits, falses(mod(-length(bits), 8)))
    return UInt8[foldl((a, b) -> (a << 1) | b, bits[i:i+7]; init = 0x00) for i in 1:8:length(bits)]
end

@testset "GT3X · log.bin" begin
    @testset "12 bits YXZ (ejemplo de la especificación)" begin
        b = hexbytes("00 60 08 EB D0 07 00 9E BF 00 70 08 EB F0")
        s = unpack12!(NTuple{3,Int16}[], b, 1, 3)
        @test s == [(8, 6, -323), (9, 7, -321), (8, 7, -321)]       # (x, y, z)
        @test round.(Float64.(collect(s[1])) ./ 341; digits = 3) == [0.023, 0.018, -0.947]
        # ida y vuelta con el empaquetador de test, incluidos extremos
        xs = [(0, 0, 0), (2047, -2048, 1), (-1, 1000, -1000), (123, -456, 789)]
        @test unpack12!(NTuple{3,Int16}[], pack12(xs), 1, length(xs)) == xs
    end

    @testset "PARAMETERS: coma flotante SSP" begin
        @test ssp_decode(0x09400000) == 256.0       # ACCEL_SCALE del ejemplo y de los wGT3X-BT
        @test ssp_decode(0x0C400000) == 2048.0      # IMU_ACCEL_SCALE del ejemplo
        @test ssp_decode(0x00000000) == 0.0
        @test ssp_decode(0x03419581) == 0x419581 * 2.0^3 / 2^23   # BATTERY_VOLTAGE, «4.1 V» en la especificación
    end

    @testset "log.bin sintético" begin
        t0 = 1_700_000_000
        one = pack12(fill((1, 2, 3), 30))
        b = vcat(
            rec(REC_METADATA, t0, Vector{UInt8}(codeunits("{\"SubjectName\":\"X\"}"))),
            rec(REC_PARAMETERS, t0, vcat(reinterpret(UInt8, htol.(UInt16[0x0000, 0x0037]))..., reinterpret(UInt8, [htol(0x09400000)])...)),
            rec(REC_ACTIVITY, t0, one),
            rec(REC_LUX, t0, UInt8[0x10, 0x00]),
            rec(REC_ACTIVITY, t0 + 1, one),
            rec(REC_EVENT, t0 + 2, UInt8[0x08]),            # entra en idle sleep
            rec(REC_EVENT, t0 + 5, UInt8[0x09]),            # sale
            zeros(UInt8, 3),                                # relleno permitido
            rec(REC_ACTIVITY, t0 + 5, one),
            rec(REC_ACTIVITY, t0 + 6, UInt8[0x01]),         # evento USB
            rec(REC_ACTIVITY, t0 + 6, one; bad = true),     # checksum incorrecto → se descarta
            rec(REC_ACTIVITY2, t0 + 8, reinterpret(UInt8, htol.(Int16[48, -52, 332, 0, 0, 0]))),
        )
        lb = parse_logbin(b; sample_rate = 30)
        @test lb.act_seconds == UInt32[t0, t0 + 1, t0 + 5, t0 + 8]
        @test lb.act_nsamples == [30, 30, 30, 2]
        @test lb.samples[1] == (1, 2, 3) && lb.samples[end-1] == (48, -52, 332)
        @test lb.diagnostics["checksum_errors"] == 1
        @test lb.diagnostics["usb_events"] == 1
        @test lb.diagnostics["padding_bytes"] == 3
        @test lb.lux_raw == UInt16[16]
        @test ssp_decode(lb.params[(0x0000, 0x0037)]) == 256.0
        g = activity_gaps(lb, 30; first_second = t0, end_second = t0 + 10)
        @test g.kind == ["idle_sleep", "missing", "partial", "missing"]
        @test datetime2unix.(g.start) == [t0 + 2, t0 + 6, t0 + 8, t0 + 9]
        @test datetime2unix.(g.stop) == [t0 + 5, t0 + 8, t0 + 9, t0 + 10]
        @test g.n_missing == [90, 60, 28, 30]
        # Contabilidad: registradas + huecos = esperadas
        @test sum(lb.act_nsamples) + sum(g.n_missing) == 10 * 30
    end

    f = raw("VTL-02-008.gt3x"); fcsv = raw("VTL-02-008RAW.csv")
    if isfile(f) && isfile(fcsv)
        @testset "archivo real VTL-02-008 frente a ActiLife" begin
            r = read_gt3x(f)
            @test r.id == "VTL-02-008" && r.format == :logbin
            @test r.scale == 256.0 && r.scale_source == "PARAMETERS"
            @test r.diagnostics["checksum_errors"] == 0
            @test issorted(r.acc.time) && allunique(r.acc.time)
            # Contabilidad exacta: muestras registradas + huecos = duración × fs
            @test nrow(r.acc) + sum(r.gaps.n_missing) == expected_samples(r) == 9 * 3600 * 30
            @test count(==("idle_sleep"), r.gaps.kind) == r.diagnostics["idle_sleep_enter_events"]

            a = read_actilife_raw(fcsv)
            @test nrow(a) == expected_samples(r)          # ActiLife rellena los huecos
            c = compare_actilife(r, a)
            @info "VTL-02-008 frente a ActiLife (redondeo lejos de cero)" c.maxdiff_recorded c.maxdiff_filled
            @test c.times_ok && c.idx_ok && c.n_recorded == nrow(r.acc)
            @test c.maxdiff_recorded == 0                 # coincidencia exacta muestra a muestra
            @test c.maxdiff_filled == 0                   # ActiLife repite la última muestra en los huecos
        end
    else
        @info "VTL-02-008 (.gt3x o RAW.csv) no disponible; se omiten tests con archivo real"
    end
end
