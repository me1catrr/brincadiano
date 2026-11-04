"""
    Sleep

Periodo de reposo nocturno sin diario para un monitor en la **cadera** con *idle sleep* (cap. 9).

Sigue las reglas de GGIR 3.3.9 para la cadera (`HASPT.algo = "HorAngle"`, función
`HASPT`), reimplementadas a mano:
1. Ángulos de los ejes respecto al plano horizontal, atan(a_i / √(a_j² + a_k²)) en grados,
   sobre la señal calibrada, diezmada a 10 Hz y con mediana móvil de 5 s (como
   `g.applymetrics`), y promediados en *epochs* de 5 s.
2. Eje longitudinal: el que está más cerca de la vertical durante el uso (en ActiGraph en
   cadera, el Y). Tumbado, ese eje queda casi horizontal: |ángulo| < 60° = «tumbado».
3. Ventana de 12:00 a 12:00: se descartan los tramos tumbados de 30 min o menos, se unen
   los que distan menos de 60 min y el periodo de reposo nocturno (SPT) es el tramo más largo.
4. Bloques de inactividad sostenida (SIB, van Hees 2015): sin cambios del ángulo z de más
   de 5° durante más de 5 min. Dentro de la SPT, SIB ≈ sueño y el resto ≈ vigilia.

Adaptaciones propias, declaradas:
- *Idle sleep*: el dispositivo solo se duerme cuando está quieto, así que la postura no
  cambia; los *epochs* sin muestras cubiertos por *idle sleep* heredan el último ángulo
  (`fill_idle`, imputación explícita con su columna `imputado`). Otros huecos no se rellenan.
- No uso: nunca cuenta como «tumbado». Si una noche tiene ≥ 2 h de no uso probable, la
  postura no informa y la ventana se estima como el tramo de no uso más largo (SPT*,
  `metodo = "retirada"`), el mismo supuesto que L5* (cap. 10): el niño se quita el
  dispositivo para dormir.
"""
module Sleep

using Arrow, DataFrames, Dates
using ..Metrics: _overlap_ms

export epoch_angles, fill_idle, longitudinal_axis, sib_vanhees2015, spt_block, sleep_nights,
       night_l5, sleep_summary, sleep_regularity_index, night_subsample, nap_candidates,
       save_angles, load_angles

const MS_H = 3_600_000

# Mediana de v[a:b] usando buf como espacio de trabajo.
function _median!(buf, v, a, b)
    m = b - a + 1
    @inbounds for j in 0:m-1
        buf[j+1] = v[a+j]
    end
    w = view(buf, 1:m)
    if isodd(m)
        return partialsort!(w, (m + 1) ÷ 2)
    else
        s = partialsort!(w, (m ÷ 2):(m ÷ 2 + 1))
        return (s[1] + s[2]) / 2
    end
end

"""
    epoch_angles(acc, gaps, fs, epoch = Second(5); origin = floor(first(acc.time), Day),
                 stop = last(acc.time), offset = (0, 0, 0), scale = (1, 1, 1)) -> DataFrame

Ángulos medios de los tres ejes (grados) por *epoch*, como `anglex/y/z` de GGIR: la señal
calibrada se diezma a 10 Hz, se suaviza con una mediana móvil de 5 s (51 muestras) **dentro
de cada tramo continuo** (sin cruzar huecos) y se promedia el ángulo en el *epoch*.
Columnas: `start`, `n` (muestras diezmadas), `idle_frac`, `missing_frac`, `anglex`,
`angley`, `anglez` (`NaN` sin muestras).
"""
function epoch_angles(acc::AbstractDataFrame, gaps::AbstractDataFrame, fs::Integer, epoch::Period = Second(5);
                      origin::DateTime = floor(first(acc.time), Day), stop::DateTime = last(acc.time),
                      offset = (0.0, 0.0, 0.0), scale = (1.0, 1.0, 1.0))
    t = acc.time
    step = max(fs ÷ 10, 1)
    newfs = step > 1 ? 10 : fs
    win = round(Int, newfs * 5); iseven(win) && (win += 1)
    half = win ÷ 2
    # Tramos continuos: saltos de más de 1,5 periodos de muestreo.
    lim = Millisecond(ceil(Int, 1500 / fs))
    seg = [1; [i for i in 2:length(t) if t[i] - t[i-1] > lim]]
    segend = [seg[2:end] .- 1; length(t)]
    idx = Int[]; segd = Int[]                                 # índices diezmados y tramo
    for (s, (a, b)) in enumerate(zip(seg, segend))
        r = a:step:b
        append!(idx, r); append!(segd, fill(s, length(r)))
    end
    nd = length(idx)
    cal = [[(Float64(col[i]) + offset[c]) * scale[c] for i in idx] for (c, col) in enumerate((acc.x, acc.y, acc.z))]
    med = [Vector{Float64}(undef, nd) for _ in 1:3]
    buf = Vector{Float64}(undef, win)
    lo = 1
    while lo <= nd                                            # mediana móvil por tramo
        hi = lo
        while hi < nd && segd[hi+1] == segd[lo]; hi += 1; end
        for c in 1:3, i in lo:hi
            med[c][i] = _median!(buf, cal[c], max(lo, i - half), min(hi, i + half))
        end
        lo = hi + 1
    end
    w = Dates.value(convert(Millisecond, epoch))
    nb = fld(Dates.value(stop - origin), w) + 1
    n = zeros(Int, nb); s = zeros(3, nb)
    @inbounds for j in 1:nd
        k = fld(Dates.value(t[idx[j]] - origin), w) + 1
        (1 <= k <= nb) || continue
        x, y, z = med[1][j], med[2][j], med[3][j]
        n[k] += 1
        s[1, k] += atand(x, sqrt(y^2 + z^2))
        s[2, k] += atand(y, sqrt(x^2 + z^2))
        s[3, k] += atand(z, sqrt(x^2 + y^2))
    end
    starts = [origin + Millisecond((k - 1) * w) for k in 1:nb]
    idle = gaps[gaps.kind .== "idle_sleep", :]; other = gaps[gaps.kind .!= "idle_sleep", :]
    fi = [_overlap_ms(a, a + Millisecond(w), idle.start, idle.stop) / w for a in starts]
    fm = [_overlap_ms(a, a + Millisecond(w), other.start, other.stop) / w for a in starts]
    ang(c) = [n[k] > 0 ? s[c, k] / n[k] : NaN for k in 1:nb]
    return DataFrame(start = starts, n = n, idle_frac = fi, missing_frac = fm,
                     anglex = ang(1), angley = ang(2), anglez = ang(3))
end

"""
    fill_idle(ep) -> DataFrame

Copia de `ep` en la que los *epochs* sin muestras cubiertos solo por *idle sleep* heredan
el último ángulo observado (el dispositivo solo se duerme si está quieto). Un hueco de otro
tipo corta la herencia. Añade `imputado::Bool`.
"""
function fill_idle(ep::AbstractDataFrame)
    out = copy(ep)
    imp = falses(nrow(ep))
    last = [NaN, NaN, NaN]
    cols = (:anglex, :angley, :anglez)
    for k in 1:nrow(ep)
        if ep.n[k] > 0
            last .= (ep.anglex[k], ep.angley[k], ep.anglez[k])
        elseif ep.missing_frac[k] > 0
            last .= NaN
        elseif ep.idle_frac[k] > 0 && !isnan(last[1])
            for (c, col) in enumerate(cols); out[k, col] = last[c]; end
            imp[k] = true
        end
    end
    out.imputado = imp
    return out
end

"""
    longitudinal_axis(ep; mask = trues(nrow(ep))) -> Int

Eje (1 = x, 2 = y, 3 = z) con mayor mediana de |ángulo| en los *epochs* de `mask` con
muestras: el más cercano a la vertical en el uso (en cadera y de pie, el longitudinal).
"""
function longitudinal_axis(ep::AbstractDataFrame; mask = trues(nrow(ep)))
    meds = map((:anglex, :angley, :anglez)) do c
        v = sort!(abs.(filter(isfinite, ep[mask .& (ep.n .> 0), c])))
        isempty(v) ? -Inf : v[cld(length(v), 2)]
    end
    return argmax(collect(meds))
end

"""
    sib_vanhees2015(angle; epoch_s = 5, angle_th = 5, time_min = 5) -> BitVector

Bloques de inactividad sostenida (van Hees et al., 2015), como `HASIB` de GGIR: se marcan
los cambios de ángulo > `angle_th` grados entre *epochs* consecutivos; dos cambios separados
por más de `time_min` minutos delimitan un SIB. Un `NaN` cuenta como cambio (no se cruza un
hueco sin datos). Con menos de 2 cambios: todo SIB si hay < 10 cambios (regla de GGIR).
"""
function sib_vanhees2015(angle::AbstractVector{<:Real}; epoch_s::Real = 5, angle_th::Real = 5, time_min::Real = 5)
    N = length(angle)
    sib = falses(N)
    d = diff(angle)
    postch = findall(x -> isnan(x) || abs(x) > angle_th, d)
    q1 = length(postch) > 1 ? findall(>(time_min * 60 / epoch_s), diff(postch)) : Int[]
    if isempty(q1)
        length(postch) < 10 && (sib .= true)
    else
        for q in q1
            sib[postch[q]:postch[q+1]] .= true
        end
    end
    sib[isnan.(angle)] .= false
    return sib
end

# Tramos (inicio, fin) de los `true` consecutivos.
function _runs(v::AbstractVector{Bool})
    r = Tuple{Int,Int}[]; i = 1; N = length(v)
    while i <= N
        if v[i]
            j = i
            while j < N && v[j+1]; j += 1; end
            push!(r, (i, j)); i = j + 1
        else
            i += 1
        end
    end
    return r
end

"""
    spt_block(nomov; epoch_s = 5, min_block_min = 30, max_gap_min = 60) -> (i, j) o nothing

Reglas finales de `HASPT` de GGIR: se descartan los tramos de `nomov` de `min_block_min`
minutos o menos, se rellenan los huecos interiores de menos de `max_gap_min` minutos y se
devuelve el tramo más largo (el primero si hay empate).
"""
function spt_block(nomov::AbstractVector{Bool}; epoch_s::Real = 5, min_block_min::Real = 30, max_gap_min::Real = 60)
    v = collect(Bool, nomov)
    for (i, j) in _runs(v)
        j - i + 1 <= min_block_min * 60 / epoch_s && (v[i:j] .= false)
    end
    r = _runs(v)
    for k in 1:length(r)-1
        g0, g1 = r[k][2] + 1, r[k+1][1] - 1
        g1 - g0 + 1 < max_gap_min * 60 / epoch_s && (v[g0:g1] .= true)
    end
    r = _runs(v)
    isempty(r) && return nothing
    return r[argmax([j - i + 1 for (i, j) in r])]
end

"""
    sleep_nights(ep; axis, nonwear, threshold = 60, core = (18, 10), max_nonwear = Hour(2),
                 rec_start, rec_stop, kw...) -> DataFrame

Una fila por noche (ventana de 12:00 del día `noche` a 12:00 del siguiente) con *epochs* de
5 s ya pasados por `fill_idle`. Solo se evalúan las noches cuyo núcleo (`core`, de 18:00 a
10:00) cae dentro del registro.

- `metodo = "postura"`: SPT = tramo «tumbado» (|ángulo del eje `axis`| < `threshold`) más
  largo según `spt_block`, sin contar el no uso como tumbado; SIB (`sib_vanhees2015` sobre
  `anglez`) dentro de la SPT = sueño.
- `metodo = "retirada"`: si el núcleo tiene ≥ `max_nonwear` de no uso, SPT* = tramo de no
  uso más largo con las mismas reglas; el sueño dentro no es medible (`NaN`).
- `metodo = "sin ventana"`: no hay ningún tramo que cumpla las reglas.

Horas (`inicio_h`, `fin_h`, `punto_medio_h`) desde la medianoche del día `noche`
(p. ej., 22,5 = 22:30; 31 = 07:00 del día siguiente).
"""
function sleep_nights(ep::AbstractDataFrame; axis::Integer, nonwear::AbstractVector{Bool},
                      threshold::Real = 60, core = (18, 10), max_nonwear::Period = Hour(2),
                      rec_start::DateTime = ep.start[findfirst(>(0), ep.n)],
                      rec_stop::DateTime = ep.start[findlast(>(0), ep.n)], kw...)
    w = Dates.value(ep.start[2] - ep.start[1]); es = w / 1000
    ang = ep[!, (:anglex, :angley, :anglez)[axis]]
    rows = NamedTuple[]
    d0 = Date(ep.start[1]); d1 = Date(ep.start[end])
    for d in d0:Day(1):d1
        w0 = DateTime(d) + Hour(12); w1 = w0 + Day(1)
        c0 = DateTime(d) + Hour(core[1]); c1 = DateTime(d) + Day(1) + Hour(core[2])
        (rec_start <= c0 && rec_stop >= c1) || continue
        k = findall(t -> w0 <= t < w1, ep.start)
        incore = [c0 <= t < c1 for t in ep.start[k]]
        nw = nonwear[k]
        nw_h = count(nw .& incore) * es / 3600
        valid = .!isnan.(ang[k]) .& .!nw
        metodo = nw_h >= Dates.value(convert(Millisecond, max_nonwear)) / MS_H ? "retirada" : "postura"
        mask = metodo == "retirada" ? nw : valid .& (abs.(ang[k]) .< threshold)
        b = spt_block(mask; epoch_s = es, kw...)
        sueno = waso = efic = NaN; imp_pct = NaN
        if b === nothing
            metodo = "sin ventana"; ini = fin = DateTime(0); dur = NaN
        else
            i, j = b
            ini = ep.start[k[i]]; fin = ep.start[k[j]] + Millisecond(w)
            dur = (j - i + 1) * es / 3600
            if metodo == "postura"
                sib = sib_vanhees2015(ep.anglez[k]; epoch_s = es)
                ins = i:j
                okk = valid[ins]
                sueno = count(sib[ins] .& okk) * es / 3600
                waso = count(.!sib[ins] .& okk) * es / 3600
                efic = sueno / max(sueno + waso, eps())
                imp_pct = 100 * count(ep.imputado[k[ins]]) / length(ins)
            end
        end
        h(t) = Dates.value(t - DateTime(d)) / MS_H
        push!(rows, (noche = d, dia_semana = dayofweek(d), libre = dayofweek(d) in (5, 6),
                     metodo = metodo, inicio = ini, fin = fin,
                     inicio_h = b === nothing ? NaN : h(ini), fin_h = b === nothing ? NaN : h(fin),
                     punto_medio_h = b === nothing ? NaN : (h(ini) + h(fin)) / 2,
                     duracion_h = dur, sueno_h = sueno, waso_h = waso, eficiencia = efic,
                     nonwear_nucleo_h = nw_h, pct_spt_idle = imp_pct))
    end
    return DataFrame(rows)
end

"""
    night_l5(ep60, noche; col = :enmo_idle0, nonwear_as = :missing, hours = 5,
             min_cov = 2/3) -> Float64

Punto medio (horas desde la medianoche de `noche`) de las `hours` horas consecutivas de
menor actividad media en la ventana de 12:00 a 12:00, con *epochs* de 60 s. Ventanas con
menos de `min_cov` de minutos válidos se ignoran; `nonwear_as = :rest` cuenta el no uso
como actividad nula (L5*). Equivale al L5 nocturno de GGIR (cap. 10).
"""
function night_l5(ep::AbstractDataFrame, d::Date; col::Symbol = :enmo_idle0, nonwear_as::Symbol = :missing,
                  hours::Real = 5, min_cov::Real = 2 / 3)
    w0 = DateTime(d) + Hour(12)
    k = findall(t -> w0 <= t < w0 + Day(1), ep.start)
    v = Float64.(ep[k, col])
    nw = ep.nonwear[k]
    v[nw] .= nonwear_as == :rest ? 0.0 : NaN
    L = round(Int, 60hours)
    nw_ = length(v) - L + 1
    means = fill(Inf, max(nw_, 0))
    for i in 1:nw_
        seg = view(v, i:i+L-1)
        ok = count(isfinite, seg)
        ok >= min_cov * L || continue
        means[i] = sum(x for x in seg if isfinite(x)) / ok
    end
    (isempty(means) || !isfinite(minimum(means))) && return NaN
    # Empates (p. ej., meseta de ceros del no uso como reposo): centro de la meseta, como
    # `Circadian.lx_window`.
    best = minimum(means); bi = findfirst(<=(best + 1e-9), means)
    bj = bi
    while bj < nw_ && means[bj+1] <= best + 1e-9; bj += 1; end
    c = (bi + bj) / 2
    return (Dates.value(ep.start[k[1]] - DateTime(d)) / 60_000 + (c - 1)) / 60 + hours / 2
end

_sd(v) = length(v) < 2 ? NaN : sqrt(sum((v .- sum(v) / length(v)) .^ 2) / (length(v) - 1))
_mean(v) = isempty(v) ? NaN : sum(v) / length(v)

"""
    sleep_summary(nights) -> NamedTuple

Medias y regularidad (desviación estándar entre noches) del inicio, fin, punto medio y
duración de la ventana; sueño y eficiencia (solo `postura`); y *jet lag* social = punto
medio en noches libres (viernes y sábado) − en noches de diario.
"""
function sleep_summary(n::AbstractDataFrame)
    ok = n[n.metodo .!= "sin ventana", :]
    po = ok[ok.metodo .== "postura", :]
    mid_l = ok.punto_medio_h[ok.libre]; mid_e = ok.punto_medio_h[.!ok.libre]
    return (noches = nrow(n), postura = nrow(po), retirada = count(==("retirada"), n.metodo),
            sin_ventana = count(==("sin ventana"), n.metodo),
            inicio_h = _mean(ok.inicio_h), inicio_de_h = _sd(ok.inicio_h),
            fin_h = _mean(ok.fin_h), fin_de_h = _sd(ok.fin_h),
            punto_medio_h = _mean(ok.punto_medio_h), punto_medio_de_h = _sd(ok.punto_medio_h),
            duracion_h = _mean(ok.duracion_h), duracion_de_h = _sd(ok.duracion_h),
            sueno_h = _mean(po.sueno_h), eficiencia = _mean(po.eficiencia),
            jetlag_social_h = isempty(mid_l) || isempty(mid_e) ? NaN : _mean(mid_l) - _mean(mid_e))
end

"""
    sleep_regularity_index(start, sleep; epoch_s = 5) -> DataFrame

Índice de regularidad del sueño (SRI, Phillips et al., 2017) por pareja de días naturales
consecutivos, como `CalcSleepRegularityIndex` de GGIR 3.3.9 (parte 3):
1. `sleep` (por *epoch*: `true` sueño, `false` vigilia, `missing` no válido) se agrega a
   *epochs* de 30 s con la media redondeada (un `missing` deja el de 30 s sin valor);
2. cada día se compara con el siguiente a la misma hora:
   SRI = −100 + 200 · (coincidencias / comparaciones válidas).
Devuelve `dia` (el primero de la pareja), `SRI` y `frac_valida`. 100 = mismo patrón cada
día; 0 = al azar.
"""
function sleep_regularity_index(start::AbstractVector{DateTime}, sleep::AbstractVector; epoch_s::Real = 5)
    k = round(Int, 30 / epoch_s)
    k >= 1 || throw(ArgumentError("epoch_s debe ser ≤ 30 s"))
    # Agregación a 30 s alineada al reloj.
    b30 = Dict{DateTime,Vector{Any}}()
    for (t, v) in zip(start, sleep)
        push!(get!(b30, floor(t, Second(30)), Any[]), v)
    end
    estado = Dict{DateTime,Union{Missing,Bool}}()
    for (t, vs) in b30
        estado[t] = any(ismissing, vs) ? missing : round(sum(Float64.(vs)) / length(vs)) == 1
    end
    dias = sort(unique(Date.(keys(estado))))
    M = 2880                                                     # epochs de 30 s por día
    filas = NamedTuple[]
    for i in 1:length(dias)-1
        d0 = DateTime(dias[i]); d1 = DateTime(dias[i+1])
        iguales = 0; validas = 0
        for j in 0:M-1
            a = get(estado, d0 + Second(30j), missing); b = get(estado, d1 + Second(30j), missing)
            (ismissing(a) || ismissing(b)) && continue
            validas += 1; iguales += (a == b)
        end
        push!(filas, (dia = dias[i], SRI = validas > 0 ? -100 + 200 * iguales / validas : NaN,
                      frac_valida = validas / M))
    end
    return DataFrame(filas)
end

"""
    night_subsample(noches; cols = (:inicio_h, :fin_h, :punto_medio_h, :duracion_h)) -> DataFrame

¿Cuántas noches hacen falta para estimar el reposo nocturno de un niño? Para cada número de
noches `k` de 1 a `n − 1` y **todos** los subconjuntos de `k` noches con ventana (`metodo` distinto
de `"sin ventana"` y valores finitos), compara la media de cada columna de `cols` con la de las `n`
noches. Una fila por `k` y variable: `subconjuntos`, `referencia` (media de las `n`), y el error
absoluto `mediana_abs`, `p90_abs` y `max_abs` (mismas unidades que la columna, horas). No imputa
noches ni usa los días como fecha: solo la numeración. Lanza `ArgumentError` si hay más de 20 noches
(los subconjuntos se enumeran con máscaras de bits) o menos de 2.
"""
function night_subsample(noches::AbstractDataFrame; cols = (:inicio_h, :fin_h, :punto_medio_h, :duracion_h))
    ok = noches[noches.metodo .!= "sin ventana", :]
    n = nrow(ok)
    n >= 2 || throw(ArgumentError("hacen falta al menos 2 noches con ventana"))
    n <= 20 || throw(ArgumentError("más de 20 noches: demasiados subconjuntos"))
    filas = NamedTuple[]
    for c in cols
        v = Float64.(ok[!, c])
        all(isfinite, v) || throw(ArgumentError("la columna $(c) tiene valores no finitos"))
        ref = _mean(v)
        err = [Float64[] for _ in 1:n-1]
        for mask in 1:(2^n - 2)
            k = count_ones(mask)
            s = 0.0
            for i in 1:n
                (mask >> (i - 1)) & 1 == 1 && (s += v[i])
            end
            push!(err[k], abs(s / k - ref))
        end
        for k in 1:n-1
            e = sort(err[k])
            push!(filas, (k = k, variable = String(c), subconjuntos = length(e), referencia = ref,
                          mediana_abs = e[cld(length(e), 2)], p90_abs = e[clamp(ceil(Int, 0.9 * length(e)), 1, length(e))],
                          max_abs = e[end]))
        end
    end
    return DataFrame(filas)
end

"""
    nap_candidates(ep, noches; window = (9, 18), dur_min = (15, 240), gap_min = 0) -> DataFrame

Siestas candidatas: bloques de inactividad sostenida diurna (`sib_vanhees2015` sobre `ep.anglez`,
*epochs* de 5 s ya pasados por `fill_idle`) fuera del reposo nocturno estimado (`noches` de
`sleep_nights`), con el dispositivo puesto, que empiezan dentro de `window` (horas del reloj, `[a, b)`)
y duran entre `dur_min[1]` y `dur_min[2]` minutos. `gap_min` une dos bloques separados por menos
de esos minutos (0 = no unir). Es la misma idea que `possible_nap_window` y `possible_nap_dur` de
GGIR 3.3.9 (ayuda, sin leer su código): **no está validada** frente a siestas declaradas ni frente a
GGIR, y en la cadera una siesta con el dispositivo puesto y quieto se parece a una retirada que el
control de calidad no haya marcado (`junto_a_no_uso` avisa si hay no uso justo antes o después).

Una fila por candidata: `dia` (número de día desde el primer epoch, sin fecha), `dia_semana`,
`inicio_h`, `fin_h`, `duracion_min`, `pct_idle` (parte rellenada con el último ángulo del *idle
sleep*, que el lector no rellena) y `junto_a_no_uso`.
"""
function nap_candidates(ep::AbstractDataFrame, noches::AbstractDataFrame; window = (9, 18),
                        dur_min = (15, 240), gap_min::Real = 0)
    N = nrow(ep)
    N >= 2 || return DataFrame(dia = Int[], dia_semana = Int[], inicio_h = Float64[], fin_h = Float64[],
                               duracion_min = Float64[], pct_idle = Float64[], junto_a_no_uso = Bool[])
    es = Dates.value(ep.start[2] - ep.start[1]) / 1000
    sib = sib_vanhees2015(ep.anglez; epoch_s = es)
    iv = [(n.inicio, n.fin) for n in eachrow(noches) if n.metodo != "sin ventana"]
    nocturna = [any(a <= t < b for (a, b) in iv) for t in ep.start]
    hora = [hour(t) + minute(t) / 60 + second(t) / 3600 for t in ep.start]
    cand = sib .& .!ep.nonwear .& .!isnan.(ep.anglez) .& .!nocturna
    r = _runs(cand)
    # Une tramos separados por menos de `gap_min` minutos (los huecos no son de no uso ni de noche).
    if gap_min > 0 && length(r) > 1
        v = copy(cand)
        for k in 1:length(r)-1
            g0, g1 = r[k][2] + 1, r[k+1][1] - 1
            gap = g0:g1
            (g1 - g0 + 1) * es / 60 < gap_min && !any(ep.nonwear[gap]) && !any(nocturna[gap]) && (v[gap] .= true)
        end
        r = _runs(v)
    end
    d0 = Date(ep.start[1]); k_ = round(Int, 600 / es)       # 10 min de margen para «junto a no uso»
    filas = NamedTuple[]
    for (i, j) in r
        window[1] <= hora[i] < window[2] || continue
        dur = (j - i + 1) * es / 60
        dur_min[1] <= dur <= dur_min[2] || continue
        a = max(i - k_, 1); b = min(j + k_, N)
        push!(filas, (dia = Dates.value(Date(ep.start[i]) - d0) + 1, dia_semana = dayofweek(ep.start[i]),
                      inicio_h = hora[i], fin_h = hora[i] + dur / 60, duracion_min = dur,
                      pct_idle = 100 * count(ep.imputado[i:j]) / (j - i + 1),
                      junto_a_no_uso = any(ep.nonwear[a:b])))
    end
    isempty(filas) && return DataFrame(dia = Int[], dia_semana = Int[], inicio_h = Float64[], fin_h = Float64[],
                                       duracion_min = Float64[], pct_idle = Float64[], junto_a_no_uso = Bool[])
    return DataFrame(filas)
end

"""
    save_angles(dir, ep) / load_angles(dir) -> DataFrame

Guarda o lee los ángulos por *epoch* de 5 s en `dir/angulos_5s.arrow` (`data/interim/<ID>/`).
"""
save_angles(dir::AbstractString, ep::AbstractDataFrame) =
    Arrow.write(joinpath(dir, "angulos_5s.arrow"), ep; compress = :zstd)
load_angles(dir::AbstractString) =
    DataFrame(Arrow.Table(joinpath(dir, "angulos_5s.arrow")); copycols = true)

end # module
