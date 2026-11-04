using Brincadiano.Rutas
using Brincadiano.Rutas: formato

@testset "Rutas y clasificación" begin
    @testset "Clasificación por nombre" begin
        c = clasificar("C6005_V03.gt3x")
        @test c.proyecto == "COHORTE_A" && c.sujeto == "C6005" && c.visita == "V03" && c.destino == joinpath("cohorte_a", "gt3x", "C6005_V03.gt3x")
        c = clasificar("6005_V00.gt3x")
        @test c.proyecto == "MELIPOP" && c.sujeto == "6005" && c.destino == joinpath("melipop", "gt3x", "6005_V00.gt3x")
        @test clasificar("6005_V00RAW.csv").destino == joinpath("melipop", "actilife_raw_csv", "6005_V00RAW.csv")
        # Los C60xx_V00 no se atribuyen: ni a COHORTE_A (por el prefijo) ni a MELIPOP (por el código).
        c = clasificar("C6005_V00.agd")
        @test c.proyecto == "pendiente" && startswith(c.estado, "hipótesis")
        @test c.destino == joinpath("pendientes", "agd_C60xx_V00", "C6005_V00.agd")
        c = clasificar("H-01-0011sec.agd")
        @test c.proyecto == "MiopiaUNAV" && c.sujeto == "H-01-001" && c.destino == joinpath("miopia-unav", "agd", "H-01-0011sec.agd")
        @test clasificar("H-01-001.gt3x").sujeto == clasificar("H-01-0011sec.agd").sujeto
        @test id_registro("H-01-0011sec.agd") == "H-01-001" && id_registro("H-01-00110sec.agd") == "H-01-001"
        @test id_registro("C6005_V00.agd") == "C6005_V00" && id_registro("6005_V00RAW.csv") == "6005_V00"
        @test clasificar("VTL-02-008.gt3x").proyecto == "pruebas"
        @test clasificar("config_cohorte_a.csv").destino == joinpath("cohorte_a", "config", "config_cohorte_a.csv")
        @test clasificar("config_UNAV.csv").proyecto == "MiopiaUNAV"
        @test clasificar("EJEMPLO GGIR_COHORTE_AV03.xlsx").destino == joinpath("cohorte_a", "exportaciones", "2026-10-30_JMJC", "EJEMPLO GGIR_COHORTE_AV03.xlsx")
        @test clasificar("EJEMPLO_BBDD_Datos acel.xlsx").proyecto == "MiopiaUNAV"
        @test clasificar("desconocido.bin").proyecto == "pendiente"
        @test formato("x.GT3X") == "gt3x" && formato("ABCRAW.csv") == "actilife_raw_csv" && formato("config_x.csv") == "ggir_config"
    end

    @testset "Búsqueda en la estructura" begin
        mktempdir() do d
            raw = joinpath(d, "raw")
            mkpath(joinpath(raw, "cohorte_a", "gt3x")); mkpath(joinpath(raw, "pendientes", "agd_C60xx_V00")); mkpath(joinpath(raw, ".oculto"))
            touch(joinpath(raw, "cohorte_a", "gt3x", "C6001_V03.gt3x"))
            touch(joinpath(raw, "6005_V00.gt3x"))                     # raíz (transición)
            touch(joinpath(raw, "pendientes", "agd_C60xx_V00", "C6001_V00.agd"))
            touch(joinpath(raw, ".oculto", "X.gt3x"))
            withenv("BRINCADIANO_DATA" => d) do
                @test raw_dir() == raw && interim_dir() == joinpath(d, "interim")
                @test basename.(raw_files("gt3x")) == ["6005_V00.gt3x", "C6001_V03.gt3x"]
                @test basename.(raw_files("gt3x"; proyecto = "COHORTE_A")) == ["C6001_V03.gt3x"]
                @test basename.(raw_files("agd")) == ["C6001_V00.agd"]
                @test find_raw("C6001_V03.gt3x") == joinpath(raw, "cohorte_a", "gt3x", "C6001_V03.gt3x")
                @test find_raw("X.gt3x") === nothing && find_raw("nada.gt3x") === nothing
                touch(joinpath(raw, "cohorte_a", "C6001_V03.gt3x"))      # mismo nombre dos veces
                @test_throws ErrorException find_raw("C6001_V03.gt3x")
            end
        end
    end

    # Integridad de data/raw/ (solo si está en este equipo): nombres únicos, originales en solo
    # lectura y, en la estructura por proyecto, cada archivo en la carpeta que le corresponde.
    @testset "Integridad de data/raw/" begin
        raiz = raw_dir()
        if !isdir(raiz)
            @info "data/raw/ no está en este equipo: se omiten los tests de integridad"
        else
            todos = String[]
            for (dir, dirs, fs) in walkdir(raiz)
                filter!(x -> !startswith(x, "."), dirs)
                append!(todos, joinpath.(dir, filter(f -> !startswith(f, "."), fs)))
            end
            @test allunique(basename.(todos))
            @test all(f -> (filemode(f) & 0o222) == 0, todos)       # a-w
            for f in todos
                rel = relpath(f, raiz)
                dirname(rel) == "" && continue                        # raíz: aún sin migrar
                startswith(rel, "exportaciones") && continue          # envíos sin migrar
                @test rel == clasificar(basename(f)).destino
            end
        end
    end
end
