using Brincadiano.Spectral, Brincadiano.Circadian, Brincadiano.Model
using DataFrames, Dates

@testset "Spectral" begin
    # Ruido determinista (sin Random): pseudoaleatorio en [−0,5, 0,5].
    ruido(i) = mod(sin(12.9898i) * 43758.5453, 1) - 0.5

    @testset "Lomb–Scargle generalizado" begin
        # 7 días cada 10 min: ritmo de 24 h (amplitud 2, máximo a las 15 h) y otro de 8 h.
        tc = collect(0:1/6:168 - 1/6)
        f(t) = 5 + 2cos(2π * (t - 15) / 24) + 0.5cos(2π * t / 8)
        yc = f.(tc)
        g = gls_periodogram(tc, yc, freq_grid(168; fmin = 0.5, fmax = 4))
        @test g.frecuencia[argmax(g.potencia)] ≈ 1.0 atol = 0.02
        g1 = gls_periodogram(tc, yc, [1.0, 3.0])
        @test g1.amplitud ≈ [2.0, 0.5] atol = 1e-9            # días completos: armónicos ortogonales
        @test g1.acrofase_h[1] ≈ 15.0 atol = 1e-9
        # Sin las noches (0–6 h, como un dispositivo retirado): una sinusoide pura se sigue
        # recuperando exacta, pero la componente de 8 h se fuga a 24 h (ventana espectral).
        t = [h for h in tc if mod(h, 24) >= 6]
        y24 = 5 .+ 2cos.(2π .* (t .- 15) ./ 24)
        @test gls_periodogram(t, y24, [1.0]).amplitud[1] ≈ 2.0 atol = 1e-9
        @test abs(gls_periodogram(t, f.(t), [1.0]).amplitud[1] - 2.0) > 0.05
        # Es un ajuste por mínimos cuadrados: misma RSS que `Model.ols` con media + cos + sin.
        y = f.(t)
        yr = y .+ [ruido(i) for i in eachindex(y)]
        ω = 2π / 24
        X = hcat(ones(length(t)), cos.(ω .* t), sin.(ω .* t))
        o = ols(X, yr)
        tss = sum((yr .- sum(yr) / length(yr)) .^ 2)
        @test gls_periodogram(t, yr, [1.0]).potencia[1] ≈ 1 - o.rss / tss rtol = 1e-8
        # Los NaN se ignoran.
        yn = copy(y24); yn[1:10] .= NaN
        @test gls_periodogram(t, yn, [1.0]).amplitud[1] ≈ 2.0 atol = 1e-9
    end

    @testset "ventana espectral" begin
        t = [h for h in 0:1/6:168 if mod(h, 24) >= 6]
        w = spectral_window(t, [0.0, 0.5, 1.0, 1.5])
        @test w[1] ≈ 1.0
        @test w[3] > 5w[2]                          # las noches sin dato crean un pico en 1 cpd
        treg = collect(0:1/6:168 - 1/6)
        @test spectral_window(treg, [1.0])[1] < 1e-20  # muestreo completo: sin pico
    end

    @testset "χ² y su relación con IS" begin
        @test chi2_quantile(10, 0.95) ≈ 18.307 rtol = 0.01
        @test chi2_quantile(100, 0.99) ≈ 135.807 rtol = 0.01
        # 7 días por minuto: actividad alta de 9 a 21 h, con ruido; dos horas ausentes.
        d0 = DateTime(2024, 1, 1)
        st = [d0 + Minute(m) for m in 0:7*1440-1]
        v = [(9 <= hour(s) < 21 ? 50.0 : 5.0) + 10ruido(i) for (i, s) in enumerate(st)]
        ep = DataFrame(start = st, enmo_idle0 = v, nonwear = falses(length(st)))
        ep.nonwear[600:720] .= true
        dias = Set(Date.(st))
        x = Circadian._hourly_series(ep, sort(collect(dias)))
        q = chi2_periodogram(x, 1.0; periods_h = 2:36)
        @test q.periodo_h[argmax(q.Qp .- q.Qp_crit)] == 24
        @test q.Qp[q.periodo_h .== 24][1] ≈ count(isfinite, x) * is_iv(ep, dias).IS rtol = 1e-10
        # IS por armónicos: con datos completos, la suma es IS; con la onda cuadrada, el
        # primer armónico domina y los pares (12, 6 h) casi no aportan.
        epc = copy(ep); epc.nonwear .= false
        h = harmonic_is(epc, dias)
        @test h.IS ≈ is_iv(epc, dias).IS
        @test h.IS_armonicos ≈ h.IS rtol = 1e-10
        @test h.armonicos.IS_k[1] > 0.7h.IS
        @test h.armonicos.IS_k[2] < 0.01h.IS
        @test h.armonicos.acrofase_h[1] ≈ 15.0 atol = 0.3          # centro del bloque 9–21 h
        # Con los 12 armónicos se reconstruye exactamente el perfil horario.
        rec = [h.media + sum(r.amplitud * cos(2π * r.k * (hh + 0.5 - r.acrofase_h) / 24) for r in eachrow(h.armonicos))
               for hh in 0:23]
        @test rec ≈ h.perfil rtol = 1e-10
        @test isfinite(harmonic_is(ep, dias).IS_armonicos)          # con horas ausentes, aproximada
    end

    @testset "momentos armónicos del timing" begin
        # Dosis concentrada a las 14:00 con ancla a las 3 h → θ₁ = 11 h, R₁ = 1.
        d0 = DateTime(2024, 1, 1)
        st = [d0 + Second(s) for s in 0:2*86400-1]
        v = [hour(s) == 14 && minute(s) == 0 && second(s) < 10 ? 600.0 : 20.0 for s in st]
        ep = DataFrame(start = st, enmo = v, nonwear = falses(length(st)))
        dias = Set(Date.(st))
        ht = harmonic_timing(ep, dias, 3.0; threshold = 464.6, K = 2)
        @test ht.segundos == 20
        @test ht.armonicos.R[1] ≈ 1.0 atol = 1e-3
        @test ht.armonicos.θ_h[1] ≈ 11.0 + 4.5 / 3600 atol = 1e-6
        @test ht.armonicos.R[2] ≈ 1.0 atol = 1e-3
        # Igual que la media circular de Circadian sobre esos segundos.
        φ = [hours_to_phase(hour(s) + minute(s) / 60 + second(s) / 3600, 3.0) for s in st[v .>= 464.6]]
        @test circular_mean(φ, v[v .>= 464.6]).R ≈ ht.armonicos.R[1]
        # Sin umbral, toda la actividad: R₁ pequeño (fondo plano de 20 m𝑔).
        @test harmonic_timing(ep, dias, 3.0).armonicos.R[1] < 0.05
        # Cobertura desigual: actividad constante, pero el primer día solo se observa de 12 a
        # 24 h. Sin corregir, la media apunta a la tarde (R₁ > 0); corregida, R₁ ≈ 0.
        epp = DataFrame(start = st, enmo = fill(20.0, length(st)),
                        nonwear = [Date(s) == Date(d0) && hour(s) < 12 for s in st])
        @test harmonic_timing(epp, dias, 3.0).armonicos.R[1] > 0.2
        @test harmonic_timing(epp, dias, 3.0; cobertura = true).armonicos.R[1] < 1e-9
        # Con cobertura completa, corregir no cambia nada.
        @test harmonic_timing(ep, dias, 3.0; threshold = 464.6, cobertura = true).armonicos.θ_h[1] ≈ ht.armonicos.θ_h[1]
    end
    @testset "espectro del movimiento" begin
        # 2 min a 30 Hz: el primer minuto, una oscilación de 2,5 Hz y 0,3 𝑔 en z; el segundo,
        # 12 Hz y 0,05 𝑔. Se quitan 10 muestras del segundo 70 (una ventana incompleta).
        fs = 30; t0 = DateTime(2024, 1, 1, 10)
        k = 0:(120fs - 1)
        zz = [i < 60fs ? 1 + 0.3sin(2π * 2.5 * i / fs) : 1 + 0.05sin(2π * 12 * i / fs) for i in k]
        acc = DataFrame(time = [t0 + Millisecond(round(Int, 1000i / fs)) for i in k],
                        x = zeros(length(k)), y = zeros(length(k)), z = zz)
        acc = acc[setdiff(1:nrow(acc), 70fs+1:70fs+10), :]
        ws = window_spectra(acc, fs, Second(5); psd = true)
        r = ws.rasgos
        @test nrow(r) == 23                                  # 24 ventanas, una incompleta
        @test !any(r.start .== t0 + Second(70))
        a = r[r.start .< t0 + Second(60), :]
        @test all(abs.(a.f_dom .- 2.5) .<= 30 / 256)          # cadencia: 2,5 Hz ± 1 bin
        @test all(a.frac_1 .> 0.99)
        @test a.potencia ≈ fill(0.3^2 / 2 * 1e6, nrow(a)) rtol = 0.02   # Parseval: varianza de la ventana
        @test all(a.entropia .< 0.3)
        b = r[r.start .>= t0 + Second(60), :]
        @test all(b.frac_3 .> 0.99)                          # 12 Hz: banda de 10–15 Hz
        @test size(ws.psd) == (23, 129)
        # Ruido blanco: entropía alta.
        wn = DataFrame(time = acc.time, x = zeros(nrow(acc)), y = zeros(nrow(acc)),
                       z = [1 + 0.01ruido(i) for i in 1:nrow(acc)])
        e = window_spectra(wn, fs).rasgos.entropia
        @test sum(e) / length(e) > 0.85
    end
end
