using Brincadiano.Viz
using DataFrames, Dates

@testset "Viz (funciones puras)" begin
    @testset "vector_norm" begin
        @test vector_norm([3.0], [4.0], [0.0]) == Float32[5]
        @test vector_norm(DataFrame(x = Float32[0, 1], y = Float32[0, 0], z = Float32[1, 0])) == Float32[1, 1]
    end

    @testset "binned_summary" begin
        o = DateTime(2020, 1, 1)
        # Dos muestras en el primer minuto, ninguna en el segundo, una en el tercero.
        t = [o + Second(5), o + Second(50), o + Minute(2) + Second(1)]
        s = binned_summary(t, [1.0, 3.0, 2.0], Minute(1); origin = o)
        @test nrow(s) == 3
        @test s.start == [o, o + Minute(1), o + Minute(2)]
        @test s.n == [2, 0, 1]
        @test (s.vmin[1], s.vmax[1], s.vmean[1]) == (1.0, 3.0, 2.0)
        @test isnan(s.vmin[2]) && isnan(s.vmax[2]) && isnan(s.vmean[2])   # hueco: NaN, no relleno
        @test s.vmean[3] == 2.0
        # stop alarga la serie con intervalos vacíos
        s2 = binned_summary(t, [1.0, 3.0, 2.0], Minute(1); origin = o, stop = o + Minute(4))
        @test nrow(s2) == 5 && s2.n[4:5] == [0, 0]
        @test_throws DimensionMismatch binned_summary(t, [1.0], Minute(1))
    end

    @testset "days_since y day_labels" begin
        o = DateTime(2024, 1, 6)            # sábado
        @test days_since(o + Hour(36), o) == 1.5
        @test day_labels(o, 3) == ["día 1 (sáb)", "día 2 (dom)", "día 3 (lun)"]
    end

    @testset "tamaños de figura en puntos" begin
        @test figure_size(:texto, 200) == (ANCHO_TEXTO, 200.0)
        @test ANCHO_TEXTO ≈ 453.54 atol = 0.01                 # 160 mm
        @test ANCHO_COMPLETO ≈ 731.34 atol = 0.01              # 258 mm
        @test figure_size(:completo, 900)[2] == ALTO_MAX       # nunca más alta que la página
        @test_throws ArgumentError figure_size(:otro, 100)
    end

    @testset "procedencia_tex" begin
        base = (; tipo = :datos, registros = ["6005_V00"], funcion = "plot_days24h",
                script = "scripts/04_qc.jl", entradas = ["data/interim/6005_V00/acc.arrow"],
                salida = "outputs/qc/figures/6005_V00_dias24h_qc", commit = "abc1234",
                fecha = Date(2026, 10, 5))
        t = procedencia_tex(; base...)
        @test occursin("\\BRprocedencia{datos}{6005\\_V00}{plot\\_days24h}{scripts/04\\_qc.jl}" *
                       "{data/interim/6005\\_V00/acc.arrow}{outputs/qc/figures/6005\\_V00\\_dias24h\\_qc}" *
                       "{abc1234}{2026-10-05}", t)
        @test count(==('\n'), t) == 2                                   # comentario + una llamada
        @test occursin("{abc1234+}", procedencia_tex(; base..., commit = "abc1234+"))
        s = procedencia_tex(; base..., tipo = :sintetico, registros = String[], entradas = String[])
        @test occursin("{sintetico}{}{plot", s)
        # Solo IDs y rutas relativas: nada que pueda arrastrar datos del sujeto o del equipo.
        @test_throws ArgumentError procedencia_tex(; base..., registros = ["Juan Pérez"])
        @test_throws ArgumentError procedencia_tex(; base..., entradas = ["/Users/x/data/raw/a.gt3x"])
        @test_throws ArgumentError procedencia_tex(; base..., entradas = ["../fuera.csv"])
        @test_throws ArgumentError procedencia_tex(; base..., script = "DateOfBirth: 2018-01-01")
        @test_throws ArgumentError procedencia_tex(; base..., tipo = :otro)
        @test_throws ArgumentError procedencia_tex(; base..., registros = String[])  # datos sin registros
        @test_throws ArgumentError procedencia_tex(; base..., commit = "rm -rf")
    end

    @testset "break_at_gaps" begin
        t, y = break_at_gaps([0.0, 1.0, 5.0, 6.0], [1.0, 2.0, 3.0, 4.0]; maxstep = 1.5)
        @test length(t) == 5 && isnan(y[3]) && t[3] == 3.0
        @test y[[1, 2, 4, 5]] == [1.0, 2.0, 3.0, 4.0]
        t2, y2 = break_at_gaps([0.0, 1.0], [1.0, 2.0]; maxstep = 1.5)
        @test t2 == [0.0, 1.0] && y2 == [1.0, 2.0]
    end
end
