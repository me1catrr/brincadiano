# Tests lentos: registros de 7 días frente a sus RAW.csv de ActiLife (~900 MB cada uno).
# Se activan con BRINCADIANO_SLOW_TESTS=1.
using Brincadiano.GT3X: read_gt3x, expected_samples

@testset "GT3X · registros de 7 días frente a ActiLife (lento)" begin
    for id in ("6005_V00", "6021_V00")
        f = raw(id * ".gt3x"); fcsv = raw(id * "RAW.csv")
        (isfile(f) && isfile(fcsv)) || (@info "$id no disponible; se omite"; continue)
        @testset "$id" begin
            r = read_gt3x(f)
            @test r.diagnostics["checksum_errors"] == 0
            @test nrow(r.acc) + sum(r.gaps.n_missing) == expected_samples(r) == 7 * 86400 * 30
            # Cada hueco de idle sleep contiene al menos una entrada; puede haber entradas
            # breves sin pérdida de muestras y huecos con varios ciclos.
            @test count(==("idle_sleep"), r.gaps.kind) <= r.diagnostics["idle_sleep_enter_events"]
            a = read_actilife_raw(fcsv)
            @test nrow(a) == expected_samples(r)
            c = compare_actilife(r, a)
            @info "$id frente a ActiLife" c.n c.n_recorded c.maxdiff_recorded c.maxdiff_filled
            @test c.times_ok && c.idx_ok && c.n_recorded == nrow(r.acc)
            @test c.maxdiff_recorded == 0
            @test c.maxdiff_filled == 0
            a = nothing; GC.gc()
        end
    end
end
