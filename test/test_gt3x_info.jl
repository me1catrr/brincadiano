using Brincadiano.GT3X: parse_info, read_info, ticks2datetime, default_accel_scale, INFO_WHITELIST
using Dates

@testset "GT3X · info.txt" begin
    @testset "ticks .NET" begin
        # Ejemplo de la especificación NHANES (verificado de forma independiente
        # con datetime de Python): 634556532600000000 ticks = 2011-10-31T10:21:00
        @test ticks2datetime(634556532600000000) == DateTime(2011, 10, 31, 10, 21)
        @test ticks2datetime(634236969566163879) == DateTime(2010, 10, 26, 13, 35, 56, 616)  # trunca a ms
        @test ticks2datetime(0) === nothing
        @test ticks2datetime(10_000) == DateTime(1, 1, 1, 0, 0, 0, 1)
    end

    @testset "lista blanca y formatos" begin
        txt = """
        Serial Number: MOS2D00000000
        Device Type: wGT3XBT
        Firmware: 1.9.2
        Battery Voltage: 4,03
        Sample Rate: 30
        Start Date: 636734916000000000
        Stop Date: 0
        TimeZone: -05:30:00
        Subject Name: NO-DEBE-APARECER
        DateOfBirth: 635610667800000000
        Acceleration Scale: 256.0
        """
        info = parse_info(txt)
        @test info.serial == "MOS2D00000000"
        @test info.battery_voltage == 4.03                 # coma decimal
        @test info.sample_rate == 30
        @test info.stop === nothing                        # 0 ticks = sin fecha
        @test info.tz_offset == Minute(-330)
        @test info.accel_scale == 256.0
        # Ningún campo personal llega a la estructura
        @test !any(v -> v isa AbstractString && occursin("NO-DEBE", v),
                   getfield.(Ref(info), fieldnames(typeof(info))))
        @test !("Subject Name" in INFO_WHITELIST) && !("DateOfBirth" in INFO_WHITELIST)
        @test_throws ErrorException parse_info("Serial Number: X\n")
    end

    @testset "escala por defecto" begin
        @test default_accel_scale("MOS2E46240378") == 256.0
        @test default_accel_scale("NEO1B34100019") == 341.0
        @test_throws ErrorException default_accel_scale("XYZ0")
    end

    @testset "fixture formato antiguo (ActiGraph, MIT)" begin
        info = read_info(joinpath(FIXTURES, "legacy_sample.gt3x"))
        @test info.serial == "NEO1B34100019"
        @test info.firmware == "1.4.0"
        @test info.sample_rate == 30
        @test info.stop === nothing
    end

    f = raw("VTL-02-008.gt3x")
    if isfile(f)
        @testset "archivo real (VTL-02-008)" begin
            info = read_info(f)
            @test info.device_type == "wGT3XBT"
            @test info.sample_rate == 30
            @test info.start == DateTime(2026, 7, 22, 8)
            @test info.stop == DateTime(2026, 7, 22, 17)
            @test info.tz_offset == Minute(120)
            @test info.accel_scale == 256.0
        end
    else
        @info "data/raw/VTL-02-008.gt3x no disponible; se omiten tests con archivo real"
    end
end
