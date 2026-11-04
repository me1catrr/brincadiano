using Brincadiano.Metrics
using DataFrames, Dates

@testset "Metrics" begin
    fs = 30; t0 = DateTime(2024, 1, 6, 10)
    tt(s0, s1) = [t0 + Second(s0) + Millisecond(cld(1000k, fs)) for k in 0:(fs * (s1 - s0) - 1)]

    @testset "ENMO muestra a muestra" begin
        @test enmo([0.0, 0.0], [0.0, 0.0], [1.2, 0.8]) ≈ Float32[0.2, 0.0]
    end

    @testset "epochs: valores conocidos y huecos" begin
        # 0–10 s: norma alterna 1,2 / 0,8 g (ENMO medio 100 mg, MAD 200 mg).
        # 10–20 s: idle sleep (sin muestras). 20–30 s: reposo exacto a 1 g.
        a = tt(0, 10); c = tt(20, 30)
        z = vcat([isodd(k) ? 1.2 : 0.8 for k in 1:length(a)], ones(length(c)))
        acc = DataFrame(time = vcat(a, c), x = zeros(Float32, length(z)), y = zeros(Float32, length(z)),
                        z = Float32.(z))
        gaps = DataFrame(start = [t0 + Second(10)], stop = [t0 + Second(20)], n_missing = [300],
                         kind = ["idle_sleep"])
        e5 = epoch_metrics(acc, gaps, fs, Second(5); origin = t0)
        @test nrow(e5) == 6
        @test e5.enmo[1] ≈ 100 atol = 1e-3
        @test e5.mad[1] ≈ 200 atol = 1e-3
        @test e5.cov[1:2] ≈ [1.0, 1.0]
        @test isnan(e5.enmo[3]) && e5.idle_frac[3] ≈ 1.0 && e5.enmo_idle0[3] == 0.0
        @test e5.enmo[5] ≈ 0 atol = 1e-4
        @test e5.mad[5] ≈ 0 atol = 1e-4
        # 10 s: mitad movimiento, mitad idle sleep → enmo_idle0 es la mitad de enmo.
        e10 = epoch_metrics(acc, gaps, fs, Second(10); origin = t0 + Second(5))
        @test e10.idle_frac[1] ≈ 0.5
        @test e10.enmo_idle0[1] ≈ e10.enmo[1] / 2 atol = 1e-3
        # Un hueco sin idle sleep anula enmo_idle0.
        g2 = copy(gaps); g2.kind .= "missing"
        @test isnan(epoch_metrics(acc, g2, fs, Second(5); origin = t0).enmo_idle0[3])
    end

    @testset "no uso por epoch y resumen diario" begin
        a = tt(0, 600)
        acc = DataFrame(time = a, x = zeros(Float32, length(a)), y = zeros(Float32, length(a)),
                        z = Float32[isodd(k) ? 1.1 : 0.9 for k in 1:length(a)])
        gaps = DataFrame(start = DateTime[], stop = DateTime[], n_missing = Int[], kind = String[])
        e = epoch_metrics(acc, gaps, fs, Minute(1); origin = t0)
        runs = DataFrame(start = [t0 + Minute(5)], stop = [t0 + Minute(10)], etiqueta = ["no uso probable"])
        flag_nonwear!(e, runs)
        @test e.nonwear == [falses(5); trues(5)]
        d = daily_summary(e)
        @test d.minutos_uso == [5]
        @test d.enmo_medio_mg[1] ≈ 50 atol = 1e-3
    end

    @testset "calibración en reposo" begin
        a = tt(0, 60)
        acc = DataFrame(time = a, x = zeros(Float32, length(a)), y = zeros(Float32, length(a)),
                        z = fill(1.02f0, length(a)))
        c = calibration_check(acc, fs)
        @test c.ventanas == 6
        @test c.sesgo_mg ≈ 20 atol = 1e-3
    end
end

@testset "Autocalibración (van Hees 2014)" begin
    # Puntos de reposo en 26 orientaciones de la esfera, distorsionados con un desplazamiento
    # y una escala conocidos: raw = v / scale − offset. Cada orientación, 10 s constantes.
    fs = 30; t0 = DateTime(2024, 1, 6)
    dirs = [(x, y, z) for x in -1:1, y in -1:1, z in -1:1 if (x, y, z) != (0, 0, 0)]
    off = (0.02, -0.015, 0.01); sc = (1.01, 0.99, 1.005)
    rows = DataFrame(time = DateTime[], x = Float32[], y = Float32[], z = Float32[])
    for (i, d) in enumerate(dirs)
        v = d ./ sqrt(sum(abs2, d))
        raw = ntuple(k -> v[k] / sc[k] - off[k], 3)
        for k in 0:(10fs - 1)
            push!(rows, (t0 + Second(10 * (i - 1)) + Millisecond(cld(1000k, fs)), Float32.(raw)...))
        end
    end
    c = autocalibrate(rows, fs)
    @test c.puntos == length(dirs) && c.esfera_ok && c.aplicar
    @test all(isapprox.(c.offset, off; atol = 1e-4)) && all(isapprox.(c.scale, sc; atol = 1e-4))
    @test c.error_despues_mg < 0.1 < c.error_antes_mg
    # Con una sola orientación, la esfera no está cubierta: no se aplica.
    uno = rows[1:10fs, :]; uno.time = [t0 + Millisecond(cld(1000k, fs)) for k in 0:(10fs - 1)]
    @test !autocalibrate(vcat(uno, uno), fs).aplicar

    @testset "tiempo por intensidad" begin
        t = [DateTime(2024, 1, 1) + Second(5k) for k in 0:(2 * 17280 - 1)]
        v = [hour(x) == 10 ? 200.0 : hour(x) == 11 ? 500.0 : hour(x) == 12 ? 100.0 : 10.0 for x in t]
        ep = DataFrame(start = t, enmo_idle0 = v, nonwear = [hour(x) == 13 for x in t])
        desp = [7 <= hour(x) < 22 for x in t]
        it = intensity_time(ep, desp)
        @test nrow(it) == 2
        @test it.min_despierto[1] ≈ 14 * 60                         # 07–22 h menos 1 h de no uso
        @test it.moderada[1] ≈ 60 && it.vigorosa[1] ≈ 60 && it.ligera[1] ≈ 60
        @test it.sedentario[1] ≈ 11 * 60
    end
end
