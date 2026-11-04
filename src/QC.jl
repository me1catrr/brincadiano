"""
    QC

Control de calidad por archivo (hito H4): duración y cobertura, frecuencia de muestreo
efectiva, huecos, saturación, no uso, luz y días válidos.

Todas las funciones son puras sobre tablas (`acc`, `gaps`, `lux`, `capsense` de un
`GT3XRecording`) y tienen además un método que recibe el registro completo.

**No uso** (`nonwear_windows`): criterio de van Hees et al. (2013) tal como lo implementa
GGIR 3.3.9 (`detect_nonwear_clipping`, método «2023», el predeterminado): ventanas de
15 min y ventanas largas de 60 min; un eje está quieto si su rango es < 50 mg **y** su DE
< 13 mg; la ventana es de no uso si lo están al menos 2 de los 3 ejes. Adaptación al *idle sleep*
(decisión propia): si la ventana larga casi no tiene muestras pero está cubierta por
*idle sleep* en al menos la mitad, se considera quieta en los tres ejes (el dispositivo no
detectó movimiento suficiente para despertarse); si lo que falta son datos sin *idle
sleep*, la ventana queda sin clasificar (`missing`). Los tramos inmóviles no se dan por
«no uso» sin más: `nonwear_runs` los anota con el contacto con la piel (CAPSENSE) y la
hora, y `label_runs` distingue «no uso probable» de «inmóvil con contacto».
"""
module QC

using Arrow, DataFrames, Dates
using ..GT3X: GT3XRecording

export sampling_summary, gap_summary, clipping_summary, lux_summary
export nonwear_windows, nonwear_runs, label_runs, wear_by_day, qc_report, save_runs, load_runs

const MS_H = 3_600_000

_ms(p::Period) = Dates.value(convert(Millisecond, p))

# Milisegundos de [a, b) que caen dentro de los intervalos [s, e) de `gaps` (sin solapes entre sí).
function _overlap_ms(a::DateTime, b::DateTime, s::AbstractVector{DateTime}, e::AbstractVector{DateTime})
    tot = 0
    i = searchsortedfirst(e, a + Millisecond(1))       # primer hueco que termina después de a
    @inbounds while i <= length(s) && s[i] < b
        tot += max(0, Dates.value(min(b, e[i]) - max(a, s[i])))
        i += 1
    end
    return tot
end

"""
    sampling_summary(acc, fs; start = first(acc.time), stop = last(acc.time)) -> NamedTuple

Duración del registro (de la primera a la última muestra), muestras registradas, número de
segundos con alguna muestra, fracción de esos segundos con exactamente `fs` muestras y
frecuencia efectiva media dentro de los segundos registrados (muestras/segundo).
"""
function sampling_summary(acc::AbstractDataFrame, fs::Integer;
                          start::DateTime = first(acc.time), stop::DateTime = last(acc.time))
    secs = Dates.value.(floor.(acc.time, Second)) .÷ 1000
    nsec = 0; full = 0; k = 1
    @inbounds while k <= length(secs)
        j = k
        while j < length(secs) && secs[j + 1] == secs[k]; j += 1; end
        nsec += 1; full += (j - k + 1) == fs
        k = j + 1
    end
    dur_s = Dates.value(stop - start) / 1000 + 1 / fs
    return (duracion_h = dur_s / 3600, muestras = nrow(acc), segundos_registrados = nsec,
            pct_segundos_completos = 100 * full / max(nsec, 1),
            fs_efectiva = nrow(acc) / max(nsec, 1),
            pct_tiempo_muestreado = 100 * nrow(acc) / (dur_s * fs))
end

"""
    gap_summary(gaps) -> DataFrame

Una fila por tipo de hueco: número, horas totales, mediana (s) y máximo (h).
"""
function gap_summary(gaps::AbstractDataFrame)
    d = transform(gaps, [:start, :stop] => ByRow((a, b) -> Dates.value(b - a) / 1000) => :dur_s)
    out = combine(groupby(d, :kind), nrow => :n, :dur_s => (x -> sum(x) / 3600) => :horas,
                  :dur_s => (x -> isempty(x) ? NaN : sort(x)[cld(length(x), 2)]) => :mediana_s,
                  :dur_s => (x -> isempty(x) ? NaN : maximum(x) / 3600) => :max_h)   # sin huecos: 0 filas
    return sort!(out, :kind)
end

"""
    clipping_summary(acc; threshold = 7.5) -> NamedTuple

Muestras con algún eje por encima de `threshold` g en valor absoluto (saturación; 7,5 g es
el umbral de GGIR para dispositivos de ±8 g).
"""
function clipping_summary(acc::AbstractDataFrame; threshold::Real = 7.5)
    n = count(i -> max(abs(acc.x[i]), abs(acc.y[i]), abs(acc.z[i])) > threshold, 1:nrow(acc))
    return (muestras_saturadas = n, pct_saturadas = 100 * n / max(nrow(acc), 1), umbral_g = threshold)
end

"""
    lux_summary(lux, gaps) -> NamedTuple

Registros de luz (el dispositivo solo escribe cuando cambia): número, % de ceros, valor
bruto máximo y % de registros en ese máximo (posible saturación), y % de registros que
caen dentro de un hueco de *idle sleep* (¿se registra la luz mientras duerme?).
"""
function lux_summary(lux::AbstractDataFrame, gaps::AbstractDataFrame)
    n = nrow(lux)
    n == 0 && return (registros = 0, pct_ceros = NaN, max_bruto = missing, pct_en_max = NaN,
                      pct_en_idle = NaN)
    mx = maximum(lux.lux_raw)
    idle = gaps[gaps.kind .== "idle_sleep", :]
    en_idle = count(t -> _overlap_ms(t, t + Millisecond(1), idle.start, idle.stop) > 0, lux.time)
    return (registros = n, pct_ceros = 100 * count(==(0), lux.lux_raw) / n, max_bruto = Int(mx),
            pct_en_max = 100 * count(==(mx), lux.lux_raw) / n, pct_en_idle = 100 * en_idle / n)
end

"""
    nonwear_windows(acc, gaps, fs; medium = Minute(15), long = Minute(60),
                    sd = 0.013, range = 0.05, min_axes = 2, idle_frac = 0.5,
                    min_samples = 10fs) -> DataFrame

Clasificación de no uso por ventanas: método «2023» de GGIR (el predeterminado en GGIR
3.3.9), con adaptación al *idle sleep* (ver la documentación del módulo). Para cada
ventana de 15 min se evalúa la ventana larga de 60 min que empieza en ella; si un eje está
quieto (rango < 50 mg y DE < 13 mg), se marca ese eje en **todas** las ventanas de 15 min
que cubre. El `score` de cada ventana es el número de ejes marcados (0–3), o `missing` si
ninguna ventana larga que la cubre pudo evaluarse (huecos sin *idle sleep*). Es de no uso
si `score ≥ min_axes`, con el suavizado de GGIR (1 eje entre dos ventanas de no uso pasa a
no uso). A diferencia de GGIR, que submuestrea a unos 5 Hz, se usan todas las muestras.

Columnas: `start`, `stop`, `n` (muestras en la ventana larga que empieza en ella),
`idle_frac` y `missing_frac` (fracción de esa ventana larga en huecos de cada tipo),
`score` y `nonwear::Bool`.
"""
function nonwear_windows(acc::AbstractDataFrame, gaps::AbstractDataFrame, fs::Integer;
                         medium::Period = Minute(15), long::Period = Minute(60),
                         sd::Real = 0.013, range::Real = 0.05, min_axes::Integer = 2,
                         idle_frac::Real = 0.5, min_samples::Integer = 10fs,
                         origin::DateTime = floor(first(acc.time), medium),
                         stop::DateTime = last(acc.time))
    t = acc.time
    M = _ms(medium); L = _ms(long)
    span_w = max(1, L ÷ M)                                     # ventanas medias por ventana larga
    nw = fld(Dates.value(stop - origin), M) + 1
    idle = gaps[gaps.kind .== "idle_sleep", :]
    other = gaps[gaps.kind .!= "idle_sleep", :]
    starts = [origin + Millisecond((k - 1) * M) for k in 1:nw]
    n = zeros(Int, nw); fi = zeros(nw); fm = zeros(nw)
    flags = falses(nw, 3); evaluated = falses(nw)
    axes = (acc.x, acc.y, acc.z)
    for k in 1:nw
        a = max(starts[k], first(t)); b = min(starts[k] + Millisecond(L), stop + Millisecond(1))
        span = max(Dates.value(b - a), 1)
        i0 = searchsortedfirst(t, a); i1 = searchsortedlast(t, b - Millisecond(1))
        n[k] = max(0, i1 - i0 + 1)
        fi[k] = _overlap_ms(a, b, idle.start, idle.stop) / span
        fm[k] = _overlap_ms(a, b, other.start, other.stop) / span
        cover = k:min(nw, k + span_w - 1)
        if n[k] >= min_samples
            evaluated[cover] .= true
            for (j, v) in enumerate(axes)
                lo = Inf; hi = -Inf; m = 0.0; q = 0.0; cnt = 0
                @inbounds for i in i0:i1                       # Welford: media y varianza
                    x = Float64(v[i]); cnt += 1
                    lo = min(lo, x); hi = max(hi, x)
                    d = x - m; m += d / cnt; q += d * (x - m)
                end
                sdv = cnt > 1 ? sqrt(q / (cnt - 1)) : 0.0
                (hi - lo < range) && (sdv < sd) && (flags[cover, j] .= true)
            end
        elseif fi[k] >= idle_frac                               # dormido: sin movimiento detectado
            evaluated[cover] .= true
            flags[cover, :] .= true
        end
    end
    score = Union{Missing,Int}[evaluated[k] ? count(flags[k, :]) : missing for k in 1:nw]
    flag = [!ismissing(s) && s >= min_axes for s in score]
    for k in 2:nw-1                                             # suavizado de GGIR
        if !ismissing(score[k]) && score[k] == min_axes - 1 && flag[k - 1] && flag[k + 1]
            flag[k] = true
        end
    end
    return DataFrame(start = starts, stop = starts .+ Millisecond(M), n = n, idle_frac = fi,
                     missing_frac = fm, score = score, nonwear = flag)
end

"""
    nonwear_runs(windows, capsense; night = (22, 8)) -> DataFrame

Une las ventanas de no uso consecutivas en tramos (`start`, `stop`, `horas`) y los anota
con la fracción de registros CAPSENSE con contacto (`contacto`), la fracción del tramo en
horario nocturno (`nocturno`, de `night[1]` a `night[2]` h) y la de *idle sleep*.
"""
function nonwear_runs(w::AbstractDataFrame, capsense::AbstractDataFrame; night = (22, 8))
    runs = DataFrame(start = DateTime[], stop = DateTime[], horas = Float64[],
                     idle_frac = Float64[], contacto = Float64[], nocturno = Float64[])
    k = 1
    while k <= nrow(w)
        if !w.nonwear[k]; k += 1; continue; end
        j = k
        while j < nrow(w) && w.nonwear[j + 1]; j += 1; end
        a, b = w.start[k], w.stop[j]
        sel = (capsense.time .>= a) .& (capsense.time .< b)
        cont = any(sel) ? count(capsense.worn[sel]) / count(sel) : NaN
        mins = a:Minute(1):(b - Minute(1))
        noc = count(m -> (hour(m) >= night[1]) || (hour(m) < night[2]), mins) / length(mins)
        push!(runs, (a, b, Dates.value(b - a) / MS_H, sum(w.idle_frac[k:j]) / (j - k + 1), cont, noc))
        k = j + 1
    end
    return runs
end

"""
    label_runs(runs, contacto_ref; ratio = 0.25) -> Vector{String}

Etiqueta cada tramo inmóvil: «no uso probable» si su fracción de contacto es menor que
`ratio` veces la del propio registro mientras se usa (`contacto_ref`), e «inmóvil con
contacto» en otro caso (posible sueño o reposo con el dispositivo puesto). Se compara con
el propio registro porque CAPSENSE no marca contacto de forma fiable ni siquiera en
movimiento (hito H3: 58 % en `6005_V00`).
"""
label_runs(runs::AbstractDataFrame, contacto_ref::Real; ratio::Real = 0.25) =
    [isnan(c) ? "sin CAPSENSE" : c < ratio * contacto_ref ? "no uso probable" : "inmóvil con contacto"
     for c in runs.contacto]

"""
    wear_by_day(acc, gaps, runs, labels; valid_hours = 10) -> DataFrame

Por día natural: horas cubiertas por el registro, horas de no uso probable, horas
inmóviles con contacto, horas sin datos (huecos que no son *idle sleep*), horas de uso
(cubiertas − no uso probable − sin datos) y si el día es válido (≥ `valid_hours` de uso;
10 h en COHORTE_A, 8 h en GENOBOX).
"""
function wear_by_day(acc::AbstractDataFrame, gaps::AbstractDataFrame, runs::AbstractDataFrame,
                     labels::AbstractVector{<:AbstractString}; valid_hours::Real = 10)
    t0, t1 = first(acc.time), last(acc.time)
    nu = runs[labels .== "no uso probable", :]
    ic = runs[labels .== "inmóvil con contacto", :]
    other = gaps[gaps.kind .!= "idle_sleep", :]
    rows = map(Date(t0):Day(1):Date(t1)) do d
        a = max(DateTime(d), t0); b = min(DateTime(d) + Day(1), t1)
        cub = Dates.value(b - a) / MS_H
        h_nu = _overlap_ms(a, b, nu.start, nu.stop) / MS_H
        h_ic = _overlap_ms(a, b, ic.start, ic.stop) / MS_H
        h_sd = _overlap_ms(a, b, other.start, other.stop) / MS_H
        uso = cub - h_nu - h_sd
        (dia = d, cubiertas_h = cub, no_uso_h = h_nu, inmovil_contacto_h = h_ic, sin_datos_h = h_sd,
         uso_h = uso, valido = uso >= valid_hours)
    end
    return DataFrame(rows)
end

"""
    qc_report(acc, gaps, lux, capsense, fs; kwargs...) -> NamedTuple

Ejecuta todo el control de calidad de un registro y devuelve sus piezas: `muestreo`,
`huecos`, `saturacion`, `luz`, `ventanas`, `tramos` (con `etiqueta`), `contacto_ref` y
`dias`. Los argumentos con nombre se pasan a `nonwear_windows`, `label_runs` y
`wear_by_day` (`valid_hours`, `ratio`, `clip`).
"""
function qc_report(acc::AbstractDataFrame, gaps::AbstractDataFrame, lux::AbstractDataFrame,
                   capsense::AbstractDataFrame, fs::Integer;
                   valid_hours::Real = 10, ratio::Real = 0.25, clip::Real = 7.5, kw...)
    w = nonwear_windows(acc, gaps, fs; kw...)
    runs = nonwear_runs(w, capsense)
    # Contacto de referencia: CAPSENSE fuera de cualquier tramo inmóvil.
    out = trues(nrow(capsense))
    for r in eachrow(runs)
        out .&= .!((capsense.time .>= r.start) .& (capsense.time .< r.stop))
    end
    cref = any(out) ? count(capsense.worn[out]) / count(out) : NaN
    labels = label_runs(runs, cref; ratio = ratio)
    runs.etiqueta = labels
    return (muestreo = sampling_summary(acc, fs), huecos = gap_summary(gaps),
            saturacion = clipping_summary(acc; threshold = clip), luz = lux_summary(lux, gaps),
            ventanas = w, tramos = runs, contacto_ref = cref,
            dias = wear_by_day(acc, gaps, runs, labels; valid_hours = valid_hours))
end

"""
    save_runs(dir, runs) / load_runs(dir) -> DataFrame

Guarda o lee los tramos inmóviles etiquetados en `dir/qc_tramos.arrow`
(`data/interim/<ID>/`), para que los hitos siguientes no tengan que repetir el QC.
"""
save_runs(dir::AbstractString, runs::AbstractDataFrame) =
    Arrow.write(joinpath(dir, "qc_tramos.arrow"), runs; compress = :zstd)
load_runs(dir::AbstractString) =
    DataFrame(Arrow.Table(joinpath(dir, "qc_tramos.arrow")); copycols = true)

"""
    qc_report(rec::GT3XRecording; kwargs...)

Igual que el método sobre tablas, con `fs` del registro y umbral de saturación
`accel_max − 0,5 g` (7,5 g en un dispositivo de ±8 g, como GGIR; 5,5 g en uno de ±6 g).
"""
qc_report(r::GT3XRecording; kw...) =
    qc_report(r.acc, r.gaps, r.lux, r.capsense, r.info.sample_rate;
              clip = something(r.info.accel_max, 8.0) - 0.5, kw...)

end # module
