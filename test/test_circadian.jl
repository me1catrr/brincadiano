using Brincadiano.Circadian
using DataFrames, Dates

@testset "Circadian" begin
    @testset "estadística circular" begin
        c = circular_mean([0.0, π / 2])
        @test c.θ ≈ π / 4
        @test c.R ≈ sqrt(2) / 2
        @test circular_mean([0.0, π / 2], [3.0, 1.0]).θ ≈ atan(1, 3)
        @test hours_to_phase(15, 3) ≈ π
        @test phase_to_hours(π) ≈ 12
        @test hours_to_phase(1, 23) ≈ 2π * 2 / 24            # cruza la medianoche
        @test rayleigh(10, 0.0).p ≈ 1 atol = 1e-6
        r = rayleigh(100, 0.2)                               # Z = 4 → p ≈ e⁻⁴ (1 + …) ≈ 0,018
        @test r.Z ≈ 4
        @test 0.015 < r.p < 0.022
    end

    # Dos días de epochs de 60 s con actividad 50 + 40·cos(2π(h − 15)/24): máximo a las 15:00
    # y mínimo a las 3:00.
    d0 = DateTime(2024, 1, 6)
    t = [d0 + Minute(m) for m in 0:(2 * 1440 - 1)]
    act(x) = 50 + 40 * cos(2π * (hour(x) + minute(x) / 60 - 15) / 24)
    ep = DataFrame(start = t, enmo_idle0 = act.(t), nonwear = falses(length(t)))
    dias = Set([Date(d0), Date(d0) + Day(1)])

    @testset "perfil, L5/M10 y RA" begin
        p = profile24h(ep, dias)
        @test nrow(p) == 1440
        @test all(p.dias .== 2)
        r = rar_metrics(p)
        @test r.estimable
        @test r.L5_medio_h ≈ 3.0 atol = 1 / 60
        @test r.M10_medio_h ≈ 15.0 atol = 1 / 60
        @test r.RA ≈ (r.M10 - r.L5) / (r.M10 + r.L5)
        # No uso de 22:00 a 08:00: como ausente, el perfil no cubre el reposo → no estimable.
        ep2 = copy(ep); ep2.nonwear = [hour(x) >= 22 || hour(x) < 8 for x in ep2.start]
        @test !rar_metrics(profile24h(ep2, dias; nonwear_as = :missing)).estimable
        rr = rar_metrics(profile24h(ep2, dias; nonwear_as = :rest))
        @test rr.estimable
        @test rr.L5 == 0.0
        # Meseta de ceros de 22:00 a 08:00: el centro de las ventanas de 5 h empatadas
        # (inicios de 22:00 a 03:00) da un punto medio a las 3:00, el centro de la meseta.
        @test rr.L5_medio_h ≈ 3.0 atol = 1 / 60
    end

    @testset "IS e IV" begin
        x = is_iv(ep, dias)
        @test x.IS ≈ 1.0 atol = 1e-12                       # días idénticos
        @test x.IV ≈ (2π / 24)^2 rtol = 0.02                # coseno horario: IV = (2π/24)²
        ep3 = copy(ep); ep3.enmo_idle0 = [isodd(hour(s)) ? 100.0 : 0.0 for s in ep3.start]
        @test is_iv(ep3, dias).IV > 1.9                     # alterna cada hora
        @test is_iv(ep, Set([Date(d0)])).horas == 24
    end

    @testset "episodios vigorosos y timing" begin
        s0 = DateTime(2024, 1, 6, 9)
        v = [0, 500, 500, 0, 0, 500, 0, 0, 0, 500, 500, 500, 0]
        e1 = DataFrame(start = [s0 + Second(k) for k in 0:length(v)-1], enmo = Float64.(v),
                       nonwear = falses(length(v)))
        eps = vigorous_episodes(e1, Set([Date(s0)]))
        @test nrow(eps) == 2
        @test eps.inicio == [s0 + Second(1), s0 + Second(9)]
        @test eps.segundos_vigorosos == [3, 3]
        @test eps.duracion_s == [5, 3]
        @test eps.dosis == [1500.0, 1500.0]
        # Todos los episodios 6 h después del ancla → θ = 6 h y R = 1.
        tm = timing(eps, 3.0)
        @test tm.n == 2
        @test tm.θ_h ≈ 6.0 atol = 0.01
        @test tm.R ≈ 1.0 atol = 1e-6
        @test tm.R2 ≈ 1.0 atol = 1e-6                    # todos a la misma hora: también R₂ = 1
        # Dos grupos opuestos (+6 h y +18 h): R ≈ 0 pero R₂ ≈ 1.
        b = copy(eps); b.medio = [DateTime(2024, 1, 6, 9), DateTime(2024, 1, 6, 21)]
        tb = timing(b, 3.0)
        @test tb.R < 1e-6
        @test tb.R2 ≈ 1.0 atol = 1e-6
        # Bootstrap por días: tres días con episodios a +5, +6 y +7 h.
        bd = DataFrame(medio = [DateTime(2024, 1, 6 + k, 8 + k) for k in 0:2], dosis = [1.0, 1.0, 1.0])
        bt = bootstrap_days(bd, 3.0; B = 500)
        @test bt.B == 500
        @test 4.9 <= bt.θ_ic_h[1] <= 6.0 <= bt.θ_ic_h[2] <= 7.1
        @test bootstrap_days(bd, 3.0; B = 500) == bt          # reproducible con la misma semilla
        @test nrow(timing_by_day(bd, 3.0)) == 3
        # Ancla por noche: el ancla y los episodios se retrasan 1 h cada día. Con ancla
        # fija, θ se dispersa (R < 1); con el ancla de cada noche, θ = 5 h y R = 1.
        an = [DateTime(2024, 1, 6 + k, 3 + k) for k in 0:2]
        @test timing(bd, an).θ_h ≈ 5.0 atol = 1e-9
        @test timing(bd, an).R ≈ 1.0 atol = 1e-9
        @test timing(bd, 3.0).R < 0.99
        # Episodio antes de la primera ancla: se refiere a la siguiente; sin ancla a < 36 h, NaN.
        pre = DataFrame(medio = [DateTime(2024, 1, 5, 20), DateTime(2024, 1, 12, 9)], dosis = [1.0, 1.0])
        φ = episode_phases(pre, an)
        @test phase_to_hours(φ[1]) ≈ 17.0 atol = 1e-9
        @test isnan(φ[2])
        @test timing(pre, an).n == 1
        @test bootstrap_days(bd, an; B = 200).amplitud_h < 1e-9
        # Cobertura: el primer día solo se observa de 12 a 24 h. Un episodio a las 9 h (visto
        # 1 día) y otro a las 15 h (visto 2 días), con la misma dosis: corregidos, el de las 9
        # pesa el doble que el de las 15.
        c0 = DateTime(2024, 1, 1)
        ec = DataFrame(start = [c0 + Second(s) for s in 0:2*86400-1])
        ec.nonwear = [Date(t) == Date(c0) && hour(t) < 12 for t in ec.start]
        cob = clock_coverage(ec, Set(Date.(ec.start)))
        @test cob[9 * 3600 + 1] == 1 && cob[15 * 3600 + 1] == 2
        epc = DataFrame(medio = [c0 + Day(1) + Hour(9), c0 + Hour(15)], dosis = [10.0, 10.0])
        epc.peso = coverage_weights(epc, cob)
        @test epc.peso == [10.0, 5.0]
        @test timing(epc, 0.0; peso = :peso).θ_h < timing(epc, 0.0).θ_h   # el peso se desplaza hacia las 9 h
        @test tm.R_cos ≈ tm.R * cos(tm.θ)
        @test tm.R_sin ≈ tm.R * sin(tm.θ)
    end

    @testset "fiabilidad del timing (ICC)" begin
        # Tres niños × 4 días, todo sobre el eje C (θ = 0 h ⇒ C = R, S = 0). Medias 0,2 / 0,4 /
        # 0,6 y desvíos ±0,1: MSW = 12·0,01/9, MSB = 4·0,08/2 = 0,16, n₀ = 4 ⇒
        # σ²_entre = (0,16 − 0,01333)/4 = 0,03667 e ICC = 0,03667/(0,03667 + 0,01333) = 0,7333.
        d = DataFrame(id = repeat(["a", "b", "c"]; inner = 4), θ_h = zeros(12),
                      R = vcat((m .+ [0.1, -0.1, 0.1, -0.1] for m in (0.2, 0.4, 0.6))...))
        r = timing_reliability(d)
        @test (r.n_ninos, r.n_dias, r.k) == (3, 12, 4.0)
        @test r.var_dentro_C ≈ 0.12 / 9
        @test r.var_entre_C ≈ (0.16 - 0.12 / 9) / 4
        @test r.icc_C ≈ 0.7333333333 atol = 1e-9
        @test r.icc_vector ≈ r.icc_C                      # S ≡ 0 no aporta varianza
        @test isnan(r.icc_S)
        @test r.icc_k_C ≈ 4 * r.icc_C / (1 + 3 * r.icc_C)  # Spearman–Brown con k = 4 días
        @test timing_reliability(d; k = 1).icc_k_C ≈ r.icc_C
        # Mismo vector en el eje S (θ = 6 h): el resultado no depende de la orientación.
        d6 = transform(d, :θ_h => ByRow(_ -> 6.0) => :θ_h)
        @test timing_reliability(d6).icc_S ≈ r.icc_C atol = 1e-9
        # Sin variación dentro del niño, ICC = 1; sin diferencias entre niños, 0.
        dd = DataFrame(id = repeat(["a", "b"]; inner = 3), θ_h = zeros(6), R = [0.2, 0.2, 0.2, 0.6, 0.6, 0.6])
        @test timing_reliability(dd).icc_C ≈ 1
        de = DataFrame(id = repeat(["a", "b"]; inner = 2), θ_h = zeros(4), R = [0.3, 0.5, 0.3, 0.5])
        @test timing_reliability(de).icc_C == 0
        # Grupos desiguales: n₀ < media de días.
        du = DataFrame(id = ["a", "a", "a", "b", "b"], θ_h = zeros(5), R = [0.2, 0.3, 0.25, 0.6, 0.7])
        @test 0 < timing_reliability(du).icc_C < 1
        # Un solo niño: solo la varianza dentro del niño.
        r1 = timing_reliability(d[d.id .== "a", :])
        @test r1.n_ninos == 1 && isnan(r1.icc_C) && r1.var_dentro_C ≈ 0.04 / 3
    end
end
