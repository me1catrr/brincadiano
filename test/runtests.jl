using Brincadiano
using Test

const FIXTURES = joinpath(@__DIR__, "fixtures")

# Archivo original por nombre, en cualquier carpeta de data/raw/ (Brincadiano.Rutas). Si no está,
# devuelve una ruta inexistente y el test que lo usa se omite con un aviso; con
# BRINCADIANO_REQUIRE_DATA=1, falta un original = error (para comprobar una migración).
function raw(nombre)
    p = Brincadiano.Rutas.find_raw(nombre)
    p === nothing && get(ENV, "BRINCADIANO_REQUIRE_DATA", "0") == "1" &&
        error("BRINCADIANO_REQUIRE_DATA=1 y no se encuentra $nombre en data/raw/")
    return something(p, joinpath(Brincadiano.Rutas.raw_dir(), nombre))
end
include("actilife_csv.jl")

@testset "Brincadiano" begin
    @testset "Submódulos" begin
        for m in (:Rutas, :GT3X, :AGD, :QC, :Metrics, :Circadian, :Sleep, :Model, :Spectral, :Light, :Viz)
            @test isdefined(Brincadiano, m)
        end
    end
    include("test_rutas.jl")
    include("test_gt3x_info.jl")
    include("test_agd.jl")
    include("test_gt3x_logbin.jl")
    include("test_gt3x_legacy.jl")
    include("test_gt3x_interim.jl")
    include("test_impute.jl")
    include("test_viz.jl")
    include("test_qc.jl")
    include("test_metrics.jl")
    include("test_circadian.jl")
    include("test_sleep.jl")
    include("test_light.jl")
    include("test_model.jl")
    include("test_spectral.jl")
    if get(ENV, "BRINCADIANO_SLOW_TESTS", "0") == "1"
        include("test_gt3x_slow.jl")
    else
        @info "Tests lentos omitidos (activar con BRINCADIANO_SLOW_TESTS=1)"
    end
end
