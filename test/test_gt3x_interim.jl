using Brincadiano.GT3X: read_gt3x, write_interim, read_interim, expected_samples
using DataFrames, Dates

@testset "GT3X · exportación Arrow (ida y vuelta)" begin
    files = [joinpath(FIXTURES, "legacy_sample.gt3x")]
    isfile(raw("VTL-02-008.gt3x")) && push!(files, raw("VTL-02-008.gt3x"))
    for f in files
        r = read_gt3x(f)
        dir = mktempdir()                      # nunca en data/: los tests no dejan derivados
        write_interim(r, dir)
        @test sort(readdir(dir)) == sort(["acc.arrow", "gaps.arrow", "lux.arrow", "battery.arrow",
                                          "capsense.arrow", "meta.arrow", "diagnostics.arrow"])
        s = read_interim(dir)
        @test s.id == r.id && s.format == r.format && s.scale == r.scale
        for fld in (:serial, :firmware, :sample_rate, :start, :stop, :last_sample, :tz_offset)
            @test getfield(s.info, fld) == getfield(r.info, fld)
        end
        for t in (:acc, :gaps, :lux, :battery, :capsense)
            @test isequal(getfield(s, t), getfield(r, t))
        end
        @test expected_samples(s) == expected_samples(r)
        rm(dir; recursive = true)
    end
end
