# Imputación del idle sleep SOLO para análisis de sensibilidad (hito H14, scripts/95). El
# pipeline nunca imputa: los huecos se marcan y se respetan (cap. 5).

"""
    impute_idle_sleep(acc, gaps, fs; method) -> (acc2, gaps2)

Rellena los huecos de *idle sleep* como lo hacen otras herramientas, para medir cuánto
cambian los resultados si se imputan:
- `method = :zeros`: muestras (0, 0, 0), como `read.gt3x(imputeZeroes = TRUE)`;
- `method = :last`: la última muestra registrada antes del hueco, como ActiLife y pygt3x
  (comprobado frente a las exportaciones RAW.csv de ActiLife, `test/actilife_csv.jl`).
Los instantes siguen la convención del lector (segundo + ⌈1000k/fs⌉ ms). Los huecos de
otro tipo se dejan como están; `gaps2` ya no contiene los de *idle sleep*.
"""
function impute_idle_sleep(acc::AbstractDataFrame, gaps::AbstractDataFrame, fs::Integer; method::Symbol)
    method in (:zeros, :last) || throw(ArgumentError("method debe ser :zeros o :last"))
    idle = gaps[gaps.kind .== "idle_sleep", :]
    secs = Int[]; vals = NTuple{3,Float32}[]
    for g in eachrow(idle)
        s0 = floor(Int, datetime2unix(g.start)); s1 = floor(Int, datetime2unix(g.stop))
        v = (0f0, 0f0, 0f0)
        if method === :last
            i = searchsortedlast(acc.time, g.start - Millisecond(1))
            i >= 1 && (v = (acc.x[i], acc.y[i], acc.z[i]))
        end
        for s in s0:(s1 - 1)
            push!(secs, s); push!(vals, v)
        end
    end
    t = _sample_times(secs, fill(Int(fs), length(secs)), Int(fs))
    v = repeat(vals; inner = Int(fs))
    add = DataFrame(time = t, x = first.(v), y = getindex.(v, 2), z = last.(v))
    acc2 = sort!(vcat(DataFrame(acc; copycols = false), add), :time)
    return acc2, gaps[gaps.kind .!= "idle_sleep", :]
end
