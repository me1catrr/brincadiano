using Brincadiano.GT3X: read_gt3x, nhanes_lux, expected_samples
using ZipArchives: ZipWriter, zip_newfile
using Dates, DataFrames

# ZIP NHANES sintético en memoria → archivo temporal
function nhanes_zip(activity::Vector{UInt8}, lux::Vector{UInt8}; serial = "NEO1C16110020")
    io = IOBuffer()
    ZipWriter(io) do w
        zip_newfile(w, "info.txt"); write(w, "Serial Number: $serial\nFirmware: 2.2.0\nSample Rate: 30\nStart Date: 634556532600000000\nStop Date: 0\nSubject Name: NO-DEBE-APARECER\n")
        zip_newfile(w, "activity.bin"); write(w, activity)
        zip_newfile(w, "lux.bin"); write(w, lux)
    end
    path = tempname() * ".gt3x"
    write(path, take!(io))
    return path
end

@testset "GT3X · formato NHANES (activity.bin + lux.bin)" begin
    @testset "lux.bin: ejemplo de la especificación" begin
        raw = UInt16[0x0000, 0x0013, 0xFFFF, 0x0898, 0x044C]      # 00 00 13 00 FF FF 98 08 4C 04
        @test nhanes_lux.(raw, "NEO1C16110020") == [0, 0, 0, 2500, 1375]
        @test nhanes_lux(UInt16(2000), "MRA0000") == 6000          # ActiSleep+: ×3,25, máx. 6000
        @test nhanes_lux(UInt16(21), "NEO") == 26                  # 26,25 → 26
    end

    @testset "ZIP sintético con el ejemplo de activity.bin" begin
        act = hexbytes("00 60 08 EB D0 07 00 9E BF 00 70 08 EB F0")
        lux = hexbytes("00 00 13 00 FF FF 98 08 4C 04")
        f = nhanes_zip(act, lux)
        r = read_gt3x(f; id = "sintetico")
        @test r.format == :nhanes && r.scale == 341.0
        # Valores de la especificación y de LibraryTests.cs (redondeo lejos de cero)
        @test actilife_round.(r.acc.x) == [0.023, 0.026, 0.023]
        @test actilife_round.(r.acc.y) == [0.018, 0.021, 0.021]
        @test actilife_round.(r.acc.z) == [-0.947, -0.941, -0.941]
        @test r.acc.time == DateTime(2011, 10, 31, 10, 21) .+ Millisecond.([0, 34, 67])
        @test r.lux.lux == [0, 0, 0, 2500, 1375]                   # los 5 pares, incluido el último
        @test r.gaps.kind == ["partial"] && r.gaps.n_missing == [27]
        @test nrow(r.acc) + sum(r.gaps.n_missing) == expected_samples(r)
        rm(f)
    end

    @testset "fixture público de ActiGraph (legacy_sample.gt3x)" begin
        r = read_gt3x(joinpath(FIXTURES, "legacy_sample.gt3x"))
        @test r.format == :nhanes && r.scale == 341.0 && r.scale_source == "número de serie"
        @test r.diagnostics["samples"] == (47104 * 8) ÷ 36 == 10467
        # ~84,5 s de señal real y el resto, relleno de ceros hasta completar la página
        @test nrow(r.acc) == 2536 && r.diagnostics["trailing_zero_samples"] == 10467 - 2536
        @test r.gaps.kind == ["zero_padding"] && (r.gaps.start[1] - r.acc.time[end]) in Millisecond.((33, 34))
        @test nrow(r.lux) == 704 ÷ 2
        # Rango válido de la especificación: |cuenta| ≤ 2046 → |a| ≤ 6 g
        @test maximum(abs, Matrix(r.acc[:, [:x, :y, :z]])) <= 2046 / 341
        @test issorted(r.acc.time) && first(r.acc.time) == r.info.start
        @test r.acc.time[31] == r.info.start + Second(1)
        @test nrow(r.acc) + sum(r.gaps.n_missing) == expected_samples(r)
        # Plausibilidad física: la norma mediana debe estar cerca de 1 g
        nrm = sqrt.(r.acc.x .^ 2 .+ r.acc.y .^ 2 .+ r.acc.z .^ 2)
        @info "legacy_sample: norma mediana (g)" median = sort(nrm)[end ÷ 2]
        @test 0.8 < sort(nrm)[end ÷ 2] < 1.2
    end
end
