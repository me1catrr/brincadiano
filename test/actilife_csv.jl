# Lectura de exportaciones RAW.csv de ActiLife (solo para validación).
# Cabecera de 10 líneas; columnas Timestamp, Accelerometer X/Y/Z; decimal con coma;
# fechas "dd/MM/yyyy HH:mm:ss.fff". No se lee ni se muestra la cabecera.
using CSV, DataFrames, Dates

function read_actilife_raw(path)
    df = CSV.read(path, DataFrame; header = 11, decimal = ',',
                  dateformat = dateformat"dd/mm/yyyy HH:MM:SS.sss",
                  types = Dict(1 => DateTime, 2 => Float64, 3 => Float64, 4 => Float64))
    rename!(df, [:time, :x, :y, :z])
    return df
end

# ActiLife redondea a 3 decimales con empates lejos de cero (convención .NET
# MidpointRounding.AwayFromZero), p. ej. 16/256 = 0,0625 → 0,063.
actilife_round(x) = round(Float64(x), RoundNearestTiesAway; digits = 3)

"""
    compare_actilife(r::GT3XRecording, a::DataFrame) -> NamedTuple

Alinea nuestro registro con una exportación RAW.csv de ActiLife por índice de muestra.
Devuelve la diferencia absoluta máxima en las muestras registradas (`maxdiff_recorded`)
y en las filas que ActiLife rellena durante los huecos (`maxdiff_filled`, comparando con
la última muestra registrada: ActiLife la repite), además de comprobaciones de tiempos.
El relleno hacia delante se hace **solo aquí**, para validar; el lector no imputa.
"""
function compare_actilife(r, a)
    fs = r.info.sample_rate
    N = nrow(a)
    ms = Dates.value.(r.acc.time .- r.info.start)                 # ms desde el inicio
    idx = (ms .÷ 1000) .* fs .+ (mod.(ms, 1000) .* fs) .÷ 1000 .+ 1  # válido con instantes «techo»
    expected_times = [r.info.start + Second((i - 1) ÷ fs) + Millisecond(cld(1000 * ((i - 1) % fs), fs)) for i in 1:N]
    recorded = falses(N); recorded[idx] .= true
    out = (; times_ok = a.time == expected_times, idx_ok = allunique(idx) && maximum(idx) <= N,
           n = N, n_recorded = count(recorded))
    maxrec = 0.0; maxfill = 0.0
    for c in (:x, :y, :z)
        full = fill(NaN, N)
        full[idx] .= actilife_round.(r.acc[!, c])
        for i in 2:N                                                # relleno hacia delante (solo test)
            isnan(full[i]) && (full[i] = full[i-1])
        end
        d = abs.(full .- a[!, c])
        maxrec = max(maxrec, maximum(d[recorded]))
        any(.!recorded) && (maxfill = max(maxfill, maximum(d[.!recorded])))
    end
    return merge(out, (; maxdiff_recorded = maxrec, maxdiff_filled = maxfill))
end
