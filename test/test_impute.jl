using Brincadiano.GT3X
using DataFrames, Dates

@testset "Imputación del idle sleep (solo sensibilidad)" begin
    fs = 30
    t0 = DateTime(2024, 1, 1, 12)
    ms(k) = Millisecond(cld(1000k, fs))
    # 1 s registrado, 2 s de idle sleep, 1 s registrado.
    t = vcat([t0 + ms(k) for k in 0:fs-1], [t0 + Second(3) + ms(k) for k in 0:fs-1])
    acc = DataFrame(time = t, x = fill(0.1f0, 2fs), y = fill(0.2f0, 2fs), z = vcat(fill(0.9f0, fs), fill(1.0f0, fs)))
    gaps = DataFrame(start = [t0 + Second(1)], stop = [t0 + Second(3)], n_missing = [2fs], kind = ["idle_sleep"])
    a0, g0 = impute_idle_sleep(acc, gaps, fs; method = :zeros)
    @test nrow(a0) == 4fs && issorted(a0.time) && nrow(g0) == 0
    k = (a0.time .>= t0 + Second(1)) .& (a0.time .< t0 + Second(3))
    @test count(k) == 2fs && all(a0.x[k] .== 0) && all(a0.z[k] .== 0)
    @test a0.time[fs + 1] == t0 + Second(1) && a0.time[fs + 2] == t0 + Second(1) + Millisecond(34)
    a1, _ = impute_idle_sleep(acc, gaps, fs; method = :last)
    @test all(a1.z[k] .== 0.9f0) && all(a1.x[k] .== 0.1f0)   # la última muestra antes del hueco
    @test nrow(acc) == 2fs                                     # el original no cambia
    g2 = DataFrame(start = [t0 + Second(1)], stop = [t0 + Second(3)], n_missing = [2fs], kind = ["missing"])
    @test nrow(impute_idle_sleep(acc, g2, fs; method = :zeros)[1]) == 2fs   # otros huecos: intactos
    @test_throws ArgumentError impute_idle_sleep(acc, gaps, fs; method = :media)
end
