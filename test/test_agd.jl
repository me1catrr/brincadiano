using Brincadiano.AGD
using Brincadiano.GT3X: ticks2datetime
using Dates

# AGD sintético mínimo (mismo esquema que ActiLife 6.13): no es de ningún participante.
# Lleva campos del sujeto inventados para comprobar que la lista blanca no los deja salir.
function _agd_sintetico(path)
    t0 = 638000000000000000                       # ticks .NET
    filas = join(["($(t0 + k * 10_000_000), $(10k), $(2k), 0, $(k % 2), $(k == 2 ? "NULL" : 5k), 0, 1, 0, 0)" for k in 0:4], ", ")
    sql = """
    CREATE TABLE settings (settingID INTEGER PRIMARY KEY, settingName VARCHAR(64), settingValue VARCHAR(8192));
    INSERT INTO settings (settingName, settingValue) VALUES
      ('softwarename', 'ActiLife'), ('softwareversion', '6.13.3'), ('devicename', 'wGT3XBT'),
      ('deviceserial', 'MOS0TEST0001'), ('deviceversion', '1.9.2'), ('filter', 'Normal'),
      ('modenumber', '61'), ('epochlength', '1'), ('epochcount', '5'), ('startdatetime', '$t0'),
      ('stopdatetime', '$(t0 + 5 * 10_000_000)'), ('downloaddatetime', '0'), ('batteryvoltage', '4,07'),
      ('original sample rate', '30'), ('finished', 'true'),
      ('subjectname', 'NO_DEBE_SALIR'), ('dateOfBirth', '631000000000000000'), ('sex', 'X');
    CREATE TABLE data (dataTimestamp INTEGER, axis1 REAL, axis2 REAL, axis3 REAL, steps REAL, lux REAL,
      inclineOff REAL, inclineStanding REAL, inclineSitting REAL, inclineLying REAL);
    INSERT INTO data VALUES $filas;
    CREATE TABLE filters (filterID INTEGER PRIMARY KEY, filterStartTimestamp INTEGER, filterStopTimestamp INTEGER);
    INSERT INTO filters (filterStartTimestamp, filterStopTimestamp) VALUES ($t0, $(t0 + 2 * 10_000_000));
    CREATE TABLE logDiaryTimes (id INTEGER PRIMARY KEY AUTOINCREMENT, onDateTimestamp BIGINT NOT NULL,
      offDateTimestamp BIGINT NOT NULL, category VARCHAR(8192) NULL);
    INSERT INTO logDiaryTimes (onDateTimestamp, offDateTimestamp, category) VALUES ($t0, $t0, 'texto libre');
    """
    run(pipeline(`sqlite3 $path`; stdin = IOBuffer(sql)))
    return t0
end

@testset "AGD" begin
    @testset "Lista blanca de settings" begin
        @test !("subjectname" in AGD_SETTINGS_WHITELIST)
        @test !("dateOfBirth" in AGD_SETTINGS_WHITELIST)
        info = AGD.parse_agd_settings(Dict("deviceserial" => "MOS1", "epochlength" => "60",
            "startdatetime" => "638000000000000000", "subjectname" => "X", "batteryvoltage" => "4,1"))
        @test info.serial == "MOS1" && info.epoch_s == 60 && info.battery_voltage == 4.1
        @test info.download === nothing
        @test_throws ErrorException AGD.parse_agd_settings(Dict("deviceserial" => "MOS1"))
    end
    if Sys.which("sqlite3") === nothing
        @info "sqlite3 no está en el PATH: se omite la lectura de un AGD sintético"
    else
        mktempdir() do dir
            path = joinpath(dir, "sintetico con espacio.agd")
            t0 = _agd_sintetico(path)
            r = read_agd(path)
            @test r.info.serial == "MOS0TEST0001" && r.info.device == "wGT3XBT" && r.info.firmware == "1.9.2"
            @test r.info.filter == "Normal" && r.info.mode == 61 && r.info.epoch_s == 1 && r.info.sample_rate == 30
            @test r.info.start == ticks2datetime(t0) && r.info.download === nothing && r.info.finished === true
            @test r.info.battery_voltage == 4.07
            @test nrow(r.epochs) == 5 && r.epochs.time[2] - r.epochs.time[1] == Second(1)
            @test r.epochs.axis1 == [0, 10, 20, 30, 40] && r.epochs.steps == [0, 1, 0, 1, 0]
            @test ismissing(r.epochs.lux[3]) && r.epochs.lux[4] == 15
            @test nrow(r.filters) == 1 && r.filters.stop[1] - r.filters.start[1] == Second(2)
            @test r.rows["logDiaryTimes"] == 1 && ismissing(r.rows["proximity"])
            # Nada del sujeto en lo que devuelve el lector.
            @test !occursin("NO_DEBE_SALIR", sprint(show, MIME"text/plain"(), r))
            @test !any(occursin("NO_DEBE_SALIR", string(getfield(r.info, f))) for f in fieldnames(AGDInfo))
        end
    end
end
