"""
    Metrics

Métricas de aceleración desde la señal RAW (hito H5): ENMO y MAD agregadas en *epochs*
de 1, 5 y 60 s alineados al reloj, con la cobertura de cada *epoch* explícita, y una
comprobación de calibración en reposo.

Definiciones (r = ‖a‖ = √(x² + y² + z²), en g; resultados en mg):
- ENMO (van Hees et al., 2013): media en el *epoch* de max(r − 1, 0).
- MAD (Vähä-Ypyä et al., 2015): media en el *epoch* de |r − r̄|, con r̄ la media del *epoch*.

Huecos: los *epochs* solo usan muestras registradas. Nada se rellena en silencio:
- `enmo` y `mad` se calculan con las muestras registradas (`NaN` si no hay ninguna);
- `enmo_idle0` es la única imputación, explícita: trata el tiempo de *idle sleep* del
  *epoch* como movimiento nulo (ENMO = 0), porque el dispositivo no detectó movimiento
  suficiente para despertarse. Es `NaN` si falta tiempo por otras causas.
"""
module Metrics

using Arrow, DataFrames, Dates

export enmo, epoch_metrics, flag_nonwear!, stationary_points, calibration_check, autocalibrate,
       daily_summary, intensity_time, save_epochs, load_epochs

"""
    enmo(x, y, z) -> Vector{Float32}

ENMO muestra a muestra, en g: max(‖a‖ − 1, 0).
"""
enmo(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    Float32[max(sqrt(Float64(a)^2 + Float64(b)^2 + Float64(c)^2) - 1, 0) for (a, b, c) in zip(x, y, z)]

# Milisegundos de [a, b) dentro de los intervalos ordenados y disjuntos [s, e).
function _overlap_ms(a::DateTime, b::DateTime, s::AbstractVector{DateTime}, e::AbstractVector{DateTime})
    tot = 0
    i = searchsortedfirst(e, a + Millisecond(1))
    @inbounds while i <= length(s) && s[i] < b
        tot += max(0, Dates.value(min(b, e[i]) - max(a, s[i])))
        i += 1
    end
    return tot
end

"""
    epoch_metrics(acc, gaps, fs, epoch; origin = floor(first(acc.time), Day),
                  stop = last(acc.time)) -> DataFrame

*Epochs* consecutivos de duración `epoch`, alineados al reloj desde `origin` (medianoche
del primer día), hasta `stop`. Columnas:
- `start`; `n` (muestras registradas); `cov` (`n` / muestras esperadas);
- `idle_frac`, `missing_frac`: fracción del *epoch* en huecos de *idle sleep* y de otro tipo;
- `enmo`, `mad` (mg, solo muestras registradas; `NaN` si `n = 0`);
- `enmo_idle0` (mg): ENMO con el tiempo de *idle sleep* como movimiento nulo
  (= suma de ENMO / muestras esperadas); `NaN` si `missing_frac > 0` o si el *epoch*
  queda fuera del registro.

`offset` y `scale` (de `autocalibrate`) aplican la calibración, (a + offset)·scale, antes
de calcular la norma; por defecto no se corrige nada.
"""
function epoch_metrics(acc::AbstractDataFrame, gaps::AbstractDataFrame, fs::Integer, epoch::Period;
                       origin::DateTime = floor(first(acc.time), Day), stop::DateTime = last(acc.time),
                       offset = (0.0, 0.0, 0.0), scale = (1.0, 1.0, 1.0))
    w = Dates.value(convert(Millisecond, epoch))
    nb = fld(Dates.value(stop - origin), w) + 1
    expected = fs * w / 1000
    t = acc.time; x = acc.x; y = acc.y; z = acc.z
    n = zeros(Int, nb); se = zeros(nb); sr = zeros(nb)
    first_i = zeros(Int, nb); last_i = zeros(Int, nb)
    r = Vector{Float64}(undef, length(t))
    @inbounds for i in eachindex(t)
        k = fld(Dates.value(t[i] - origin), w) + 1
        (1 <= k <= nb) || continue
        ri = sqrt(((Float64(x[i]) + offset[1]) * scale[1])^2 + ((Float64(y[i]) + offset[2]) * scale[2])^2 +
                  ((Float64(z[i]) + offset[3]) * scale[3])^2)
        r[i] = ri
        n[k] += 1; se[k] += max(ri - 1, 0); sr[k] += ri
        first_i[k] == 0 && (first_i[k] = i); last_i[k] = i
    end
    mad = fill(NaN, nb)
    @inbounds for k in 1:nb                                  # segunda pasada: |r − r̄|
        n[k] == 0 && continue
        m = sr[k] / n[k]; s = 0.0
        for i in first_i[k]:last_i[k]
            s += abs(r[i] - m)
        end
        mad[k] = 1000 * s / n[k]
    end
    starts = [origin + Millisecond((k - 1) * w) for k in 1:nb]
    idle = gaps[gaps.kind .== "idle_sleep", :]; other = gaps[gaps.kind .!= "idle_sleep", :]
    fi = [_overlap_ms(a, a + Millisecond(w), idle.start, idle.stop) / w for a in starts]
    fm = [_overlap_ms(a, a + Millisecond(w), other.start, other.stop) / w for a in starts]
    t0 = first(t)
    inside = [a >= floor(t0, Second) for a in starts]       # epochs anteriores al inicio: fuera
    enmo_mg = [n[k] > 0 ? 1000 * se[k] / n[k] : NaN for k in 1:nb]
    idle0 = [inside[k] && fm[k] == 0 && (n[k] > 0 || fi[k] > 0) ? 1000 * se[k] / expected : NaN for k in 1:nb]
    return DataFrame(start = starts, n = n, cov = n ./ expected, idle_frac = fi, missing_frac = fm,
                     enmo = enmo_mg, mad = mad, enmo_idle0 = idle0)
end

"""
    flag_nonwear!(epochs, runs; label = "no uso probable") -> epochs

Añade la columna `nonwear::Bool`: el *epoch* cae (en más de la mitad) dentro de un tramo
inmóvil de `QC.qc_report` con esa etiqueta. Así los módulos siguientes pueden excluir la
retirada del dispositivo (p. ej., L5 en el hito H6).
"""
function flag_nonwear!(ep::AbstractDataFrame, runs::AbstractDataFrame; label::AbstractString = "no uso probable")
    rr = runs[runs.etiqueta .== label, :]
    w = nrow(ep) > 1 ? Dates.value(ep.start[2] - ep.start[1]) : 1
    ep.nonwear = [_overlap_ms(a, a + Millisecond(w), rr.start, rr.stop) > w / 2 for a in ep.start]
    return ep
end

"""
    stationary_points(acc, fs; window = Second(10), sd = 0.013) -> Matrix{Float64}

Medias (x, y, z) de las ventanas de `window` con muestreo completo, desviación estándar
< `sd` g en los tres ejes y |media| < 2 g: los puntos de reposo, en los que el sensor solo
mide la gravedad. Base de `calibration_check` y `autocalibrate`.
"""
function stationary_points(acc::AbstractDataFrame, fs::Integer; window::Period = Second(10), sd::Real = 0.013)
    w = Dates.value(convert(Millisecond, window)); need = fs * w ÷ 1000
    t = acc.time; o = floor(first(t), Second)
    pts = NTuple{3,Float64}[]
    i = 1
    while i <= length(t)
        k = fld(Dates.value(t[i] - o), w)
        j = i
        while j < length(t) && fld(Dates.value(t[j + 1] - o), w) == k; j += 1; end
        if j - i + 1 == need
            ms = map((acc.x, acc.y, acc.z)) do v
                m = sum(Float64, @view v[i:j]) / need
                (m, sqrt(sum(x -> (Float64(x) - m)^2, @view v[i:j]) / (need - 1)))
            end
            all(p -> p[2] < sd && abs(p[1]) < 2, ms) && push!(pts, (ms[1][1], ms[2][1], ms[3][1]))
        end
        i = j + 1
    end
    return reduce(vcat, ([p[1] p[2] p[3]] for p in pts); init = zeros(0, 3))
end

"""
    calibration_check(acc, fs; window = Second(10), sd = 0.013) -> NamedTuple

Comprueba la calibración en reposo: en los puntos de `stationary_points`, resume la
desviación de ‖a‖ respecto de 1 g (mediana y cuartiles, en mg). Un sesgo de unos pocos mg
se traslada íntegro a ENMO (van Hees et al., 2014). Es un diagnóstico: no corrige nada.
"""
function calibration_check(acc::AbstractDataFrame, fs::Integer; kw...)
    P = stationary_points(acc, fs; kw...)
    size(P, 1) == 0 && return (ventanas = 0, sesgo_mg = NaN, p25_mg = NaN, p75_mg = NaN)
    s = sort(1000 .* (sqrt.(sum(abs2, P; dims = 2))[:, 1] .- 1))
    q(p) = s[clamp(ceil(Int, p * length(s)), 1, length(s))]
    return (ventanas = length(s), sesgo_mg = q(0.5), p25_mg = q(0.25), p75_mg = q(0.75))
end

"""
    autocalibrate(acc, fs; spherecrit = 0.3, maxiter = 1000, tol = 1e-10) -> NamedTuple

Autocalibración por gravedad local (van Hees et al., 2014), como `g.calibrate` de GGIR
3.3.9 sin término de temperatura (los archivos no la registran): se buscan los factores
`offset` y `scale` por eje tales que los puntos de reposo corregidos, (a + offset)·scale,
queden sobre la esfera de radio 1 g. Iteración de GGIR: proyectar sobre la esfera,
regresión lineal ponderada por eje, actualizar, pesos = min(1/distancia, 100).

Devuelve `offset`, `scale` (3 componentes), `error_antes_mg` y `error_despues_mg` (media de
|‖a‖ − 1| en los puntos de reposo), `puntos`, `esfera_ok` (cada eje con puntos por debajo de
−`spherecrit` y por encima de +`spherecrit` g) y `aplicar::Bool`: la corrección solo se
aplica si la esfera está bien cubierta y el error baja y queda por debajo de 10 mg, como en
GGIR.
"""
function autocalibrate(acc::AbstractDataFrame, fs::Integer; spherecrit::Real = 0.3,
                       maxiter::Integer = 1000, tol::Real = 1e-10, kw...)
    P = stationary_points(acc, fs; kw...)
    N = size(P, 1)
    err(Q) = 1000 * sum(abs.(sqrt.(sum(abs2, Q; dims = 2)) .- 1)) / size(Q, 1)
    N < 10 && return (offset = (0.0, 0.0, 0.0), scale = (1.0, 1.0, 1.0), error_antes_mg = NaN,
                      error_despues_mg = NaN, puntos = N, esfera_ok = false, aplicar = false)
    esfera = all(k -> minimum(P[:, k]) < -spherecrit && maximum(P[:, k]) > spherecrit, 1:3)
    offset = zeros(3); scale = ones(3); w = ones(N); res = Inf
    for _ in 1:maxiter
        curr = (P .+ offset') .* scale'
        closest = curr ./ sqrt.(sum(abs2, curr; dims = 2))
        och = zeros(3); sch = ones(3)
        for k in 1:3                                     # mínimos cuadrados ponderados
            x = curr[:, k]; y = closest[:, k]
            sw = sum(w); mx = sum(w .* x) / sw; my = sum(w .* y) / sw
            b = sum(w .* (x .- mx) .* (y .- my)) / sum(w .* (x .- mx) .^ 2)
            a = my - b * mx
            och[k] = a; sch[k] = b
            curr[:, k] = a .+ b .* x
        end
        offset .+= och ./ (scale .* sch)
        scale .*= sch
        r = 3 * sum(w .* (curr .- closest) .^ 2) / size(curr, 1) / sum(w)
        w = min.(1 ./ sqrt.(sum(abs2, curr .- closest; dims = 2))[:, 1], 100.0)
        abs(r - res) < tol && break
        res = r
    end
    e0 = err(P); e1 = err((P .+ offset') .* scale')
    return (offset = Tuple(offset), scale = Tuple(scale), error_antes_mg = e0, error_despues_mg = e1,
            puntos = N, esfera_ok = esfera, aplicar = esfera && e1 < e0 && e1 < 10)
end

"""
    daily_summary(ep60) -> DataFrame

Por día natural, a partir de *epochs* de 60 s con `nonwear`: minutos de uso (sin no uso
probable ni huecos sin datos) y aceleración media en ellos (`enmo_idle0`, mg), la
métrica de volumen que recomienda el consenso GRANADA.
"""
function daily_summary(ep::AbstractDataFrame)
    d = transform(ep, :start => ByRow(Date) => :dia)
    use = .!d.nonwear .& .!isnan.(d.enmo_idle0)
    d.uso = use
    return combine(groupby(d, :dia), :uso => sum => :minutos_uso,
                   [:enmo_idle0, :uso] => ((e, u) -> any(u) ? sum(e[u]) / count(u) : NaN) => :enmo_medio_mg)
end

"""
    intensity_time(ep, despierto; umbrales = (63.3, 142.6, 464.6), col = :enmo_idle0) -> DataFrame

Minutos por día natural en cada intensidad durante el tiempo despierto (`despierto`, un
`Bool` por *epoch*, p. ej. fuera del periodo de reposo nocturno), con el dispositivo puesto
(`nonwear` falso) y ENMO conocido: sedentario (< `umbrales[1]`), ligera, moderada y
vigorosa (≥ `umbrales[3]`). Umbrales por defecto: Hildebrand en cadera (mg), como
`threshold.lig/mod/vig` en la parte 5 de GGIR.
"""
function intensity_time(ep::AbstractDataFrame, despierto::AbstractVector{Bool};
                        umbrales = (63.3, 142.6, 464.6), col::Symbol = :enmo_idle0)
    w = nrow(ep) > 1 ? Dates.value(ep.start[2] - ep.start[1]) / 60_000 : 1.0
    v = ep[!, col]
    ok = despierto .& .!ep.nonwear .& isfinite.(v)
    d = DataFrame(dia = Date.(ep.start), ok = ok, v = v)
    return combine(groupby(d, :dia),
        :ok => (o -> w * count(o)) => :min_despierto,
        [:ok, :v] => ((o, x) -> w * count(o .& (x .< umbrales[1]))) => :sedentario,
        [:ok, :v] => ((o, x) -> w * count(o .& (umbrales[1] .<= x .< umbrales[2]))) => :ligera,
        [:ok, :v] => ((o, x) -> w * count(o .& (umbrales[2] .<= x .< umbrales[3]))) => :moderada,
        [:ok, :v] => ((o, x) -> w * count(o .& (x .>= umbrales[3]))) => :vigorosa)
end

"""
    save_epochs(dir, epochs, len) / load_epochs(dir, len) -> DataFrame

Guarda o lee los *epochs* de `len` segundos en `dir/epochs_<len>s.arrow`
(`data/interim/<ID>/`, zstd).
"""
save_epochs(dir::AbstractString, ep::AbstractDataFrame, len::Integer) =
    Arrow.write(joinpath(dir, "epochs_$(len)s.arrow"), ep; compress = :zstd)
load_epochs(dir::AbstractString, len::Integer) =
    DataFrame(Arrow.Table(joinpath(dir, "epochs_$(len)s.arrow")); copycols = true)

end # module
