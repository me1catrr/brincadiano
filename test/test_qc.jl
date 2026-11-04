using Brincadiano.QC
using DataFrames, Dates

# Ruido determinista sin dependencias (generador congruencial lineal; suma de 4 uniformes
# centradas y escaladas a DE 1). Evita añadir Random a las dependencias de test.
mutable struct LCG; s::UInt64; end
function randn!(g::LCG, n)
    out = Vector{Float64}(undef, n)
    for i in 1:n
        acc = 0.0
        for _ in 1:4
            g.s = g.s * 0x5851f42d4c957f2d + 0x14057b7ef767814f
            acc += (g.s >> 11) / 2.0^53 - 0.5
        end
        out[i] = acc * sqrt(3)                       # var(suma de 4 U(-½,½)) = 1/3
    end
    return out
end

# Registro sintético a 30 Hz: 1 h de movimiento, 1 h de idle sleep (sin muestras),
# 1 h quieto pero muestreando (ruido de 1 mg) y 1 h de movimiento.
function sintetico(; fs = 30, rng = LCG(1))
    t0 = DateTime(2024, 1, 6, 10)                       # sábado 10:00
    tt(h0, h1) = [t0 + Hour(h0) + Millisecond(cld(1000k, fs)) for k in 0:(fs * 3600 * (h1 - h0) - 1)]
    mov(n) = (0.2 .* randn!(rng, n), 0.2 .* randn!(rng, n), 1 .+ 0.2 .* randn!(rng, n))
    qui(n) = (0.001 .* randn!(rng, n), 0.001 .* randn!(rng, n), 1 .+ 0.001 .* randn!(rng, n))
    a = tt(0, 1); c = tt(2, 3); d = tt(3, 4)
    xa, ya, za = mov(length(a)); xc, yc, zc = qui(length(c)); xd, yd, zd = mov(length(d))
    acc = DataFrame(time = vcat(a, c, d), x = Float32.(vcat(xa, xc, xd)),
                    y = Float32.(vcat(ya, yc, yd)), z = Float32.(vcat(za, zc, zd)))
    gaps = DataFrame(start = [t0 + Hour(1)], stop = [t0 + Hour(2)], n_missing = [fs * 3600],
                     kind = ["idle_sleep"])
    lux = DataFrame(time = [t0 + Minute(5), t0 + Minute(90), t0 + Minute(200)], lux_raw = UInt16[0, 40, 2595])
    # CAPSENSE cada minuto: contacto en las horas de movimiento, sin contacto el resto.
    ct = [t0 + Minute(m) for m in 0:239]
    cap = DataFrame(time = ct, worn = [(m < 60) || (m >= 180) for m in 0:239])
    return acc, gaps, lux, cap, fs
end

@testset "QC" begin
    acc, gaps, lux, cap, fs = sintetico()

    @testset "muestreo, huecos, saturación y luz" begin
        s = sampling_summary(acc, fs)
        @test s.duracion_h ≈ 4 atol = 1e-6
        @test s.pct_segundos_completos == 100 && s.fs_efectiva == fs
        @test s.pct_tiempo_muestreado ≈ 75 atol = 1e-4
        g = gap_summary(gaps)
        @test g.kind == ["idle_sleep"] && g.n == [1] && g.horas ≈ [1.0]
        c = clipping_summary(acc; threshold = 0.5)
        @test c.muestras_saturadas > 0                    # el ruido de 0,2 g supera 0,5 g a veces
        @test clipping_summary(acc; threshold = 7.5).muestras_saturadas == 0
        l = lux_summary(lux, gaps)
        @test l.registros == 3 && l.max_bruto == 2595
        @test l.pct_ceros ≈ 100 / 3 && l.pct_en_idle ≈ 100 / 3   # el de los 90 min cae en idle sleep
    end

    @testset "no uso por ventanas" begin
        w = nonwear_windows(acc, gaps, fs)
        mid(k) = w.start[k] + Minute(7) + Second(30)
        en(h0, h1) = [k for k in 1:nrow(w) if w.start[k] + Minute(7) >= w.start[1] + Hour(h0) + Minute(30) &&
                                                   w.start[k] + Minute(7) < w.start[1] + Hour(h1) - Minute(30)]
        @test all(w.nonwear[en(1, 2)])                   # idle sleep: quieto por definición
        @test all(w.nonwear[en(2, 3)])                   # muestreando pero DE < 13 mg
        @test !any(w.nonwear[en(0, 1)]) && !any(w.nonwear[en(3, 4)])
        @test all(w.score[en(1, 2)] .== 3)
    end

    @testset "tramos, etiquetas y días" begin
        r = qc_report(acc, gaps, lux, cap, fs; valid_hours = 1)
        @test nrow(r.tramos) == 1
        @test r.tramos.horas[1] ≈ 2.0                    # 11:00–13:00 completo (método «2023»)
        @test r.contacto_ref ≈ 1.0
        @test r.tramos.etiqueta == ["no uso probable"]   # sin contacto en el tramo inmóvil
        @test nrow(r.dias) == 1 && r.dias.valido == [true]
        @test r.dias.uso_h[1] ≈ r.dias.cubiertas_h[1] - r.dias.no_uso_h[1]
        # Con contacto durante el tramo, pasa a «inmóvil con contacto».
        cap2 = copy(cap); cap2.worn .= true
        @test qc_report(acc, gaps, lux, cap2, fs).tramos.etiqueta == ["inmóvil con contacto"]
    end

    @testset "resumen de huecos sin huecos" begin
        vacio = DataFrame(start = DateTime[], stop = DateTime[], n_missing = Int[], kind = String[])
        @test nrow(gap_summary(vacio)) == 0
    end
end
