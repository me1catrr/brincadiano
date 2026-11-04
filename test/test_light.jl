using Brincadiano.Light
using DataFrames, Dates

@testset "Light" begin
    t0 = DateTime(2024, 1, 1, 12)
    @testset "serie por epoch: escalonada y huecos" begin
        lux = DataFrame(time = [t0, t0 + Second(90)], lux_raw = UInt16[10, 1000])
        gaps = DataFrame(start = [t0 + Second(150)], stop = [t0 + Second(170)], n_missing = [600], kind = ["idle_sleep"])
        le = light_epochs(lux, gaps; origin = t0, stop = t0 + Second(239))
        @test nrow(le) == 4
        @test le.lux_media[1] == 10 && le.lux_obs[1] == 1
        @test le.lux_media[2] ≈ 505                           # 30 s a 10 y 30 s a 1000
        @test le.lux_max[2] == 1000
        @test le.s_alto[2] == 30 && le.s_bajo[1] == 0 && le.s_alto[3] == 40
        @test le.lux_obs[3] ≈ 40 / 60                         # 20 s de idle sleep: desconocidos
        @test le.lux_media[3] == 1000
        le0 = light_epochs(lux, gaps; origin = t0 - Minute(1), stop = t0 + Second(59))
        @test le0.lux_obs[1] == 0 && isnan(le0.lux_media[1])  # antes del primer registro
    end

    @testset "estados y métricas diarias" begin
        # Dos días de epochs de 60 s desde el lunes a mediodía. Periodo de reposo nocturno de 22:00 a
        # 07:00; luz alta (1500) de 10:00 a 11:00 con AFMV; oscuro y activo de 16:00 a 16:30;
        # sin dato de 17:00 a 17:20; retirado de 19:00 a 19:30.
        d0 = DateTime(2024, 1, 1)
        t = [d0 + Minute(k) for k in 0:(3 * 1440 - 1)]
        h = [hour(x) + minute(x) / 60 for x in t]
        lm = [10 <= x < 11 ? 1500.0 : 16 <= x < 16.5 ? 0.0 : 200.0 for x in h]
        obs = [17 <= x < 17 + 1 / 3 ? 0.0 : 1.0 for x in h]
        le = DataFrame(start = t, lux_obs = obs, lux_media = lm, lux_max = lm)
        enmo = [(10 <= x < 11 || 16 <= x < 16.5) ? 200.0 : 20.0 for x in h]
        ep = DataFrame(start = t, enmo_idle0 = enmo, nonwear = [19 <= x < 19.5 for x in h])
        noches = DataFrame(noche = [Date(d0), Date(d0) + Day(1)], metodo = ["postura", "postura"],
                           inicio = [d0 + Hour(22), d0 + Day(1) + Hour(22)],
                           fin = [d0 + Day(1) + Hour(7), d0 + Day(2) + Hour(7)])
        st = light_states(le, ep, noches)
        e(hh) = st.estado[findfirst(==(d0 + Day(1) + Minute(round(Int, 60hh))), st.start)]
        @test e(3) == "sueño"
        @test e(10.5) == "válido"
        @test e(16.2) == "oscuro y activo"
        @test e(17.1) == "sin dato"
        @test e(19.2) == "retirado"
        dd = light_daily(st, noches)
        r = dd[dd.dia .== Date(d0) + Day(1), :][1, :]
        @test r.min_despierto == 15 * 60 - 30                 # 07:00–22:00 menos la retirada
        @test r.min_sobre_alto == 60
        @test r.primera_alta_h ≈ 10 && r.ultima_alta_h ≈ 10 + 59 / 60
        @test r.entradas_alta == 1
        @test r.afmv_min == 90 && r.afmv_exterior == 60 && r.afmv_oscuro == 30
        @test r.despertar_h ≈ 7 && r.dormir_h ≈ 22
        @test r.tras_despertar_pct_valido ≈ 100
        @test r.antes_dormir_min_bajo == 180 - 30             # 19:00–22:00 menos la retirada
        s = light_summary(dd; min_despierto = 600)
        @test s.dias >= 1
        @test s.hora_media_alta_h ≈ 10.5 atol = 0.02
    end
end
