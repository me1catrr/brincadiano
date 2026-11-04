using Brincadiano.Model
using DataFrames

@testset "Model" begin
    @testset "funciones especiales" begin
        @test Brincadiano.Model._lgamma(5.0) ≈ log(24) atol = 1e-12
        @test Brincadiano.Model._lgamma(0.5) ≈ log(sqrt(π)) atol = 1e-12
        # Con d1 = 2, P(F > f) = (1 + 2f/d2)^(−d2/2) exactamente.
        @test fdist_sf(3.0, 2, 10) ≈ (1 + 0.6)^-5 rtol = 1e-10
        @test fdist_sf(3.0873, 2, 100) ≈ 0.05 atol = 1e-4     # valor crítico tabulado
        @test fdist_sf(0.0, 2, 10) == 1.0
        @test normcdf(1.959964) ≈ 0.975 atol = 1e-6
        @test normcdf(0.0) ≈ 0.5 atol = 1e-7
        g = Rng(1)
        u = [rand(g) for _ in 1:20_000]; z = [randn(g) for _ in 1:20_000]
        @test 0.49 < sum(u) / length(u) < 0.51
        @test abs(sum(z) / length(z)) < 0.03 && 0.97 < sqrt(sum(abs2, z) / length(z)) < 1.03
    end

    @testset "mínimos cuadrados (QR)" begin
        x = collect(1.0:10.0)
        f = ols(hcat(ones(10), x), 1 .+ 2 .* x)
        @test f.β ≈ [1.0, 2.0] atol = 1e-10
        @test f.rss < 1e-18
        # Con ruido: coincide con las ecuaciones normales resueltas a mano (2 × 2).
        y = 3 .+ 0.5 .* x .+ [0.3, -0.2, 0.1, 0.4, -0.5, 0.2, -0.1, 0.0, 0.3, -0.4]
        f = ols(hcat(ones(10), x), y)
        sx, sy, sxx, sxy = sum(x), sum(y), sum(x .^ 2), sum(x .* y)
        b = (10sxy - sx * sy) / (10sxx - sx^2); a = (sy - b * sx) / 10
        @test f.β ≈ [a, b] atol = 1e-10
    end

    @testset "modelo de timing sobre cohortes sintéticas" begin
        # Insesgadez con dispersión entre niños de 4 h (DE de θ₀ ≈ 0,44 h con n = 5000).
        df = simulate_cohort(5000; A = 0.3, θ0_h = 18.0, θ_de_h = 4.0, seed = 3)
        m = fit_timing_model(df)
        @test m.p < 1e-6
        @test m.A ≈ 0.3 atol = 0.06
        @test m.θ0_h ≈ 18.0 atol = 1.5
        @test m.fase_favorable_h ≈ mod(m.θ0_h + 12, 24)
        b = bootstrap_theta0(simulate_cohort(1000; A = 0.3, θ_de_h = 4.0, seed = 4); B = 100)
        lo, hi = b.θ0_ic_h
        @test (lo <= 18 <= hi) || (hi < lo && (18 >= lo || 18 <= hi))   # el arco puede cruzar 24 h
        @test b.amplitud_h < 8
    end

    @testset "tobit por EM" begin
        g = Rng(9); n = 5000
        x = [randn(g) for _ in 1:n]
        ystar = [1 + 0.5x[i] + randn(g) for i in 1:n]
        c = 0.5; cens = ystar .<= c; y = max.(ystar, c)
        X = hcat(ones(n), x)
        @test 0.15 < count(cens) / n < 0.45
        @test abs(ols(X, y).β[2] - 0.5) > 0.05            # MCO sesgado por la censura
        t = tobit_em(X, y, cens)
        @test t.β ≈ [1.0, 0.5] atol = 0.05
        @test t.σ ≈ 1.0 atol = 0.05
    end

    @testset "potencia" begin
        pg = power_grid([300], [0.0, 0.4]; nsim = 60)
        @test pg.potencia[1] <= 0.15                       # error de tipo I ≈ 0,05
        @test pg.potencia[2] > pg.potencia[1]
    end

    @testset "error de medida en el plano (error_cs)" begin
        a = simulate_cohort(200; seed = 7)
        @test simulate_cohort(200; seed = 7, error_cs = 0.0) == a      # con 0, misma cohorte
        b = simulate_cohort(200; seed = 7, error_cs = 0.15)
        @test b.θ_verdadero_h == a.θ_verdadero_h && b.homa == a.homa  # solo cambia lo observado
        @test b.R != a.R
        # El error se añade en C y S: la diferencia de vectores tiene DE ≈ 0,15 por componente.
        φa = 2π .* a.θ_h ./ 24; φb = 2π .* b.θ_h ./ 24
        dc = b.R .* cos.(φb) .- a.R .* cos.(φa)
        @test 0.12 < sqrt(sum(abs2, dc) / length(dc)) < 0.18
    end
end
