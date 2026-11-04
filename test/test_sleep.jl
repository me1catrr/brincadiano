using Brincadiano.Sleep
using DataFrames, Dates

@testset "Sleep" begin
    @testset "ángulos por epoch y mediana móvil" begin
        # 20 s a 30 Hz con el eje Y vertical (de pie) y 20 s con el eje Z vertical (tumbado),
        # más picos aislados que la mediana de 5 s debe eliminar.
        fs = 30
        t0 = DateTime(2024, 1, 1, 12)
        N = 40fs
        t = [t0 + Millisecond(ceil(Int, 1000k / fs)) for k in 0:N-1]
        y = [k < 20fs ? 1.0 : 0.0 for k in 0:N-1]
        z = [k < 20fs ? 0.0 : 1.0 for k in 0:N-1]
        x = zeros(N); x[100] = 5.0; x[1000] = -5.0
        acc = DataFrame(time = t, x = x, y = y, z = z)
        gaps = DataFrame(start = DateTime[], stop = DateTime[], n_missing = Int[], kind = String[])
        ep = epoch_angles(acc, gaps, fs, Second(5); origin = t0)
        @test nrow(ep) == 8
        @test all(ep.n .== 50)                                # 30 Hz diezmado a 10 Hz
        @test ep.angley[1] ≈ 90 atol = 1e-9
        @test ep.anglez[8] ≈ 90 atol = 1e-9
        @test all(abs.(ep.anglex) .< 1e-9)                    # los picos no sobreviven
        @test longitudinal_axis(ep[1:4, :]) == 2
    end

    @testset "herencia del ángulo en el idle sleep" begin
        ep = DataFrame(start = [DateTime(2024, 1, 1) + Second(5k) for k in 0:5],
                       n = [50, 0, 0, 0, 0, 50], idle_frac = [0, 1, 1, 0, 1, 0.0],
                       missing_frac = [0, 0, 0, 1, 0, 0.0],
                       anglex = [10, NaN, NaN, NaN, NaN, 20], angley = [80, NaN, NaN, NaN, NaN, 70],
                       anglez = [0, NaN, NaN, NaN, NaN, 5.0])
        f = fill_idle(ep)
        @test f.anglex[2:3] == [10, 10]
        @test f.imputado == [false, true, true, false, false, false]
        @test isnan(f.anglex[4])                              # hueco sin datos: no se rellena
        @test isnan(f.anglex[5])                              # y corta la herencia
    end

    @testset "SIB de van Hees 2015" begin
        # Cambios de 10° cada epoch salvo un tramo quieto de 10 min (120 epochs de 5 s).
        a = Float64[]
        for k in 1:50; push!(a, isodd(k) ? 0 : 10); end
        append!(a, fill(10.0, 120))
        for k in 1:50; push!(a, isodd(k) ? 0 : 10); end
        s = sib_vanhees2015(a)
        @test count(s) == 122                                 # de un cambio al siguiente
        @test !any(s[1:48]) && !any(s[end-48:end])
        @test all(sib_vanhees2015(fill(3.0, 200)))            # sin cambios: todo SIB (regla de GGIR)
        @test !any(sib_vanhees2015([isodd(k) ? 0.0 : 10.0 for k in 1:200]))
    end

    @testset "reglas de la ventana (HASPT)" begin
        es = 5; m = 60 ÷ es                                   # epochs por minuto
        v = falses(24 * 60 * m)
        v[1:20m] .= true                                      # 20 min: se descarta
        v[100m:400m] .= true                                  # bloque A
        v[440m:700m] .= true                                  # hueco de 40 min: se une a A
        v[800m:900m] .= true                                  # hueco de 100 min: no se une
        b = spt_block(v)
        @test b == (100m, 700m)
        @test spt_block(falses(1000)) === nothing
        @test spt_block(trues(10)) === nothing                # 50 s: demasiado corto
    end

    @testset "noches: postura, retirada y L5 nocturno" begin
        # Dos días de epochs de 5 s: de pie (Y a 90°) salvo tumbado de 22:00 a 07:00, con
        # vueltas en la cama (z cambia 30°) cada 30 min; la 2.ª noche, sin dispositivo.
        d0 = DateTime(2024, 1, 5)                             # viernes
        t = [d0 + Second(5k) for k in 0:(3 * 17280 - 1)]
        h = [Dates.value(x - d0) / 3_600_000 for x in t]
        lying = [(22 <= mod(x, 24) || mod(x, 24) < 7) for x in h]
        night2 = [(24 + 22 <= x < 48 + 7) for x in h]
        ang_y = [l ? 5.0 : 85.0 for l in lying]
        ang_z = [l ? (isodd(floor(Int, x * 2)) ? 30.0 : 0.0) : 0.0 for (l, x) in zip(lying, h)]
        n = length(t)
        ep = DataFrame(start = t, n = fill(50, n), idle_frac = zeros(n), missing_frac = zeros(n),
                       anglex = zeros(n), angley = ang_y, anglez = ang_z, imputado = falses(n))
        nw = collect(night2)
        nn = sleep_nights(ep; axis = 2, nonwear = nw)
        @test nrow(nn) == 2
        p1 = nn[1, :]
        @test p1.metodo == "postura"
        @test p1.inicio_h ≈ 22 atol = 1e-6
        @test p1.fin_h ≈ 31 atol = 1e-6
        @test p1.punto_medio_h ≈ 26.5 atol = 1e-6
        @test p1.libre                                        # noche de viernes
        @test 8.5 < p1.sueno_h <= 9                           # SIB desde la primera vuelta (22:30)
        @test nn[2, :metodo] == "retirada"
        @test nn[2, :duracion_h] ≈ 9 atol = 1e-6
        @test isnan(nn[2, :sueno_h])
        s = sleep_summary(nn)
        @test s.postura == 1 && s.retirada == 1
        @test s.punto_medio_h ≈ 26.5 atol = 1e-6

        # L5 nocturno con epochs de 60 s: actividad mínima de 01:00 a 06:00.
        t60 = [d0 + Minute(k) for k in 0:(2 * 1440 - 1)]
        act = [(1 <= hour(x) < 6) ? 1.0 : 50.0 for x in t60]
        e60 = DataFrame(start = t60, enmo_idle0 = act, nonwear = falses(length(t60)))
        @test night_l5(e60, Date(d0)) ≈ 24 + 3.5 atol = 1e-9

        # SRI: tres días idénticos → 100; el tercero con el sueño desplazado 12 h → negativo.
        tt = [d0 + Second(5k) for k in 0:(3 * 17280 - 1)]
        dorm = [(hour(x) < 7 || hour(x) >= 22) for x in tt]
        sr = sleep_regularity_index(tt, dorm)
        @test nrow(sr) == 2 && all(sr.SRI .≈ 100) && all(sr.frac_valida .≈ 1)
        inv = [Date(x) == Date(d0) + Day(2) ? !d : d for (x, d) in zip(tt, dorm)]
        sr2 = sleep_regularity_index(tt, inv)
        @test sr2.SRI[1] ≈ 100 && sr2.SRI[2] ≈ -100
        mis = Vector{Union{Missing,Bool}}(dorm); mis[1:8640] .= missing   # 12 h sin dato el día 1
        sr3 = sleep_regularity_index(tt, mis)
        @test sr3.frac_valida[1] ≈ 0.5 && sr3.SRI[1] ≈ 100
        # Retirada de 22:00 a 08:00 como reposo: meseta de ceros → centro (03:00).
        e60.nonwear .= [(22 <= hour(x) || hour(x) < 8) for x in t60]
        e60.enmo_idle0 .= 50.0
        @test night_l5(e60, Date(d0); nonwear_as = :rest) ≈ 27 atol = 1 / 60
    end

    @testset "submuestreo de noches" begin
        # Cuatro noches con inicio 22, 23, 22, 23 h (media 22,5) y una sin ventana que se ignora.
        n = DataFrame(metodo = ["postura", "postura", "postura", "postura", "sin ventana"],
                      inicio_h = [22.0, 23.0, 22.0, 23.0, NaN], fin_h = [30.0, 31.0, 30.0, 31.0, NaN],
                      punto_medio_h = [26.0, 27.0, 26.0, 27.0, NaN], duracion_h = [8.0, 8.0, 8.0, 8.0, NaN])
        s = night_subsample(n)
        @test sort(unique(s.k)) == [1, 2, 3] && all(s.subconjuntos[s.k .== 1] .== 4)
        i = s[(s.variable .== "inicio_h"), :]
        @test all(i.referencia .≈ 22.5)
        # k = 1: cada noche difiere 0,5 h de la media; k = 2: 6 pares con errores 0 y 0,5; k = 3: 0,167.
        @test i.max_abs[i.k .== 1] ≈ [0.5]
        @test i.mediana_abs[i.k .== 1] ≈ [0.5]
        @test i.max_abs[i.k .== 2] ≈ [0.5]
        @test i.mediana_abs[i.k .== 2][1] < 1e-12
        @test i.max_abs[i.k .== 3] ≈ [0.5 / 3]
        # Una variable constante no tiene error.
        d = s[s.variable .== "duracion_h", :]
        @test all(d.max_abs .≈ 0)
        @test_throws ArgumentError night_subsample(n[1:1, :])
    end

    @testset "siestas candidatas" begin
        # Dos días de 5 s: el ángulo oscila (movimiento) salvo de 13:00 a 14:00 del día 1 (quieto) y
        # de 03:00 a 04:00 (noche, fuera del reposo nocturno estimado → se ignora).
        t0 = DateTime(2024, 1, 1)
        t = [t0 + Second(5k) for k in 0:(2 * 17280 - 1)]
        ang = [isodd(k) ? 10.0 : 0.0 for k in 0:(2 * 17280 - 1)]
        for k in 1:length(t)
            h = hour(t[k]); d = Date(t[k])
            (d == Date(t0) && h == 13) && (ang[k] = 20.0)
            (d == Date(t0) + Day(1) && h == 3) && (ang[k] = 20.0)
        end
        ep = DataFrame(start = t, anglez = ang, nonwear = falses(length(t)), imputado = falses(length(t)))
        noches = DataFrame(inicio = [DateTime(2024, 1, 2, 2)], fin = [DateTime(2024, 1, 2, 5)], metodo = ["postura"])
        s = nap_candidates(ep, noches)
        @test nrow(s) == 1
        @test s.dia[1] == 1 && 59 <= s.duracion_min[1] <= 63 && 12.9 <= s.inicio_h[1] < 13.2
        @test s.pct_idle[1] == 0 && !s.junto_a_no_uso[1]
        # Sin la noche declarada, la quietud de las 03:00 queda fuera de la ventana horaria (9–18 h).
        @test nrow(nap_candidates(ep, noches[1:0, :])) == 1
        # Una duración mínima mayor que la siesta la descarta.
        @test nrow(nap_candidates(ep, noches; dur_min = (90, 240))) == 0
        # El no uso pegado al bloque lo marca.
        nw = copy(ep); nw.nonwear[findall(x -> x == DateTime(2024, 1, 1, 14, 5), nw.start)[1]:findall(x -> x == DateTime(2024, 1, 1, 14, 6), nw.start)[1]] .= true
        @test nap_candidates(nw, noches).junto_a_no_uso[1]
    end
end
