"""
    Circadian

Ritmo de reposo-actividad y *timing* circular de la actividad vigorosa (hito H6).
Estadística circular implementada a mano (regla del proyecto).

- **Perfil de 24 h** (`profile24h`): media de la actividad por minuto del día sobre los
  días válidos, con el no uso probable como dato ausente (`nonwear_as = :missing`) o como
  reposo (`:rest`, ENMO = 0; supuesto: retirada = sueño).
- **L5 / M10 / RA** (`lx_window`, `rar_metrics`): ventana circular de 5 h menos activa y de
  10 h más activa del perfil; RA = (M10 − L5)/(M10 + L5) (Van Someren et al., 1999).
- **IS / IV** (`is_iv`): sobre medias horarias, con la generalización a horas ausentes
  documentada en la función.
- **Episodios vigorosos** (`vigorous_episodes`): *epochs* de 1 s por encima de un umbral
  de ENMO, unidos si los separan pausas breves.
- **Timing circular** (`circular_mean`, `rayleigh`, `timing`): fase de cada episodio
  respecto al punto medio de L5, media circular ponderada por la «dosis» del episodio
  → ángulo θ y concentración R.
- **Fiabilidad** (`timing_reliability`): ICC de R·cos θ y R·sin θ entre niños y entre
  días, y su fiabilidad para la media de k días (Spearman–Brown).
"""
module Circadian

using DataFrames, Dates

export profile24h, lx_window, rar_metrics, is_iv, vigorous_episodes,
       circular_mean, rayleigh, hours_to_phase, phase_to_hours, episode_phases, timing, timing_by_day,
       bootstrap_days, clock_coverage, coverage_weights, timing_reliability

const HILDEBRAND_VIGOROUS_HIP = 464.6   # mg, ActiGraph en cadera, 7–11 años (Hildebrand 2014)

# Epoch utilizable para el ritmo: día válido, no es no uso y tiene valor.
_usable(ep, valid_days, col) = [(Date(t) in valid_days) && !nw && isfinite(v)
                                for (t, nw, v) in zip(ep.start, ep.nonwear, ep[!, col])]

"""
    profile24h(ep60, valid_days; col = :enmo_idle0, nonwear_as = :missing) -> DataFrame

Perfil medio de 24 h a resolución de minuto (1440 filas: `minuto`, `media`, `dias`) con
los *epochs* de 60 s de los días válidos. `nonwear_as = :missing` excluye el no uso
probable; `:rest` lo cuenta como ENMO = 0 (supuesto explícito: el dispositivo se retira
para dormir). Minutos sin datos: `media = NaN`, `dias = 0`.
"""
function profile24h(ep::AbstractDataFrame, valid_days; col::Symbol = :enmo_idle0, nonwear_as::Symbol = :missing)
    nonwear_as in (:missing, :rest) || throw(ArgumentError("nonwear_as debe ser :missing o :rest"))
    s = zeros(1440); c = zeros(Int, 1440)
    for (t, nw, v) in zip(ep.start, ep.nonwear, ep[!, col])
        Date(t) in valid_days || continue
        m = hour(t) * 60 + minute(t) + 1
        if nw
            nonwear_as === :rest && (c[m] += 1)            # suma 0
        elseif isfinite(v)
            s[m] += v; c[m] += 1
        end
    end
    return DataFrame(minuto = 0:1439, media = [c[m] > 0 ? s[m] / c[m] : NaN for m in 1:1440], dias = c)
end

"""
    lx_window(profile, hours; mode = :min, min_cov = 2/3) -> NamedTuple

Ventana circular (cruza la medianoche) de `hours` horas con la menor (`:min`, L5) o mayor
(`:max`, M10) actividad media del perfil. Solo se aceptan ventanas con datos en al menos
`min_cov` de sus minutos. Si varias ventanas empatan (una meseta, p. ej. de ceros), se
toma la del centro del tramo de ventanas empatadas. Devuelve `valor` (mg), `inicio_h`,
`punto_medio_h` (horas del reloj, 0–24) y `cobertura` de la ventana elegida.
"""
function lx_window(p::AbstractDataFrame, hours::Real; mode::Symbol = :min, min_cov::Real = 2 / 3,
                   tol::Real = 1e-9)
    L = round(Int, 60hours); v = p.media
    vals = fill(NaN, 1440); covs = zeros(1440)
    for s in 0:1439
        idx = mod.(s:(s + L - 1), 1440) .+ 1
        ok = filter(isfinite, v[idx])
        covs[s + 1] = length(ok) / L
        covs[s + 1] >= min_cov && (vals[s + 1] = sum(ok) / length(ok))
    end
    cand = findall(isfinite, vals)
    isempty(cand) && return (valor = NaN, inicio_h = NaN, punto_medio_h = NaN, cobertura = NaN)
    best = mode === :min ? minimum(vals[cand]) : maximum(vals[cand])
    tied = [k for k in cand if abs(vals[k] - best) <= tol]
    # Empates (p. ej., una meseta de ceros): se toma el centro del tramo circular de inicios
    # empatados que contiene al primero, no el primero, para no sesgar la fase.
    k0 = first(tied); run = [k0]
    k = k0
    while (nk = mod(k, 1440) + 1) in tied && nk != k0; push!(run, nk); k = nk; end
    k = k0
    while (pk = mod(k - 2, 1440) + 1) in tied && pk != k0 && !(pk in run); pushfirst!(run, pk); k = pk; end
    sel = run[cld(length(run), 2)]
    s = sel - 1
    return (valor = vals[sel], inicio_h = s / 60, punto_medio_h = mod(s + L / 2, 1440) / 60, cobertura = covs[sel])
end

"""
    rar_metrics(profile; min_profile_cov = 0.9) -> NamedTuple

L5, M10 (valor e inicio y punto medio en horas del reloj), amplitud relativa
RA = (M10 − L5)/(M10 + L5) y la cobertura del perfil. `estimable` es falso si el perfil
tiene datos en menos de `min_profile_cov` de los minutos del día: entonces el reposo
probablemente cae en el no uso y L5 no describiría el sueño.
"""
function rar_metrics(p::AbstractDataFrame; min_profile_cov::Real = 0.9)
    cov = count(isfinite, p.media) / 1440
    l5 = lx_window(p, 5; mode = :min); m10 = lx_window(p, 10; mode = :max)
    ra = (m10.valor - l5.valor) / (m10.valor + l5.valor)
    return (L5 = l5.valor, L5_inicio_h = l5.inicio_h, L5_medio_h = l5.punto_medio_h,
            M10 = m10.valor, M10_inicio_h = m10.inicio_h, M10_medio_h = m10.punto_medio_h,
            RA = ra, cobertura_perfil = cov, estimable = cov >= min_profile_cov)
end

# Medias horarias de los días válidos (`days` ordenados), 24 por día seguidas; `NaN` si la hora
# tiene menos de `min_minutes` minutos utilizables. La usan `is_iv` y `Spectral.harmonic_is`.
function _hourly_series(ep::AbstractDataFrame, days::AbstractVector{Date}; col::Symbol = :enmo_idle0,
                        nonwear_as::Symbol = :missing, min_minutes::Integer = 30)
    pos = Dict(d => i for (i, d) in enumerate(days))     # día válido → posición
    nh = 24 * length(days)
    s = zeros(nh); c = zeros(Int, nh)
    for (t, nw, v) in zip(ep.start, ep.nonwear, ep[!, col])
        i = get(pos, Date(t), 0)
        i == 0 && continue
        k = (i - 1) * 24 + hour(t) + 1
        if nw
            nonwear_as === :rest && (c[k] += 1)
        elseif isfinite(v)
            s[k] += v; c[k] += 1
        end
    end
    return [c[k] >= min_minutes ? s[k] / c[k] : NaN for k in 1:nh]
end

"""
    is_iv(ep60, valid_days; col = :enmo_idle0, nonwear_as = :missing, min_minutes = 30) -> NamedTuple

Estabilidad interdiaria (IS) y variabilidad intradiaria (IV) sobre medias horarias de los
días válidos (hora con al menos `min_minutes` minutos utilizables; si no, ausente).

Con datos completos son las definiciones de Van Someren et al. (1999):
IS = n Σₕ (x̄ₕ − x̄)² / (p Σᵢ (xᵢ − x̄)²) y IV = n Σᵢ (xᵢ − xᵢ₋₁)² / ((n − 1) Σᵢ (xᵢ − x̄)²).
Con horas ausentes se generalizan así (decisión propia; coinciden con las anteriores si no
falta nada): IS = Σₕ nₕ (x̄ₕ − x̄)² / Σᵢ (xᵢ − x̄)², con nₕ las horas disponibles de la hora
del día h; IV usa solo las parejas de horas consecutivas disponibles,
IV = n Σ (Δx)² / (n_pares Σ (xᵢ − x̄)²).
"""
function is_iv(ep::AbstractDataFrame, valid_days; col::Symbol = :enmo_idle0, nonwear_as::Symbol = :missing,
               min_minutes::Integer = 30)
    days = sort(collect(valid_days))
    isempty(days) && return (IS = NaN, IV = NaN, horas = 0)
    x = _hourly_series(ep, days; col = col, nonwear_as = nonwear_as, min_minutes = min_minutes)
    nh = length(x)
    ok = isfinite.(x); n = count(ok)
    n < 2 && return (IS = NaN, IV = NaN, horas = n)
    xm = sum(x[ok]) / n
    tot = sum((x[ok] .- xm) .^ 2)
    hm = [(v = filter(isfinite, x[h:24:end]); (length(v), isempty(v) ? 0.0 : sum(v) / length(v))) for h in 1:24]
    IS = sum(nh_ * (m - xm)^2 for (nh_, m) in hm if nh_ > 0) / tot
    dd = Float64[]
    for k in 2:nh
        # Entre días válidos no consecutivos no hay pareja (la hora anterior no es contigua).
        contiguo = (k - 1) % 24 != 0 || days[(k - 1) ÷ 24 + 1] - days[(k - 1) ÷ 24] == Day(1)
        (ok[k] && ok[k - 1] && contiguo) && push!(dd, (x[k] - x[k - 1])^2)
    end
    IV = isempty(dd) ? NaN : n * sum(dd) / (length(dd) * tot)
    return (IS = IS, IV = IV, horas = n)
end

"""
    vigorous_episodes(ep1, valid_days; threshold = 464.6, tolerance = 2, min_seconds = 3,
                      col = :enmo) -> DataFrame

Episodios de actividad vigorosa en *epochs* de 1 s de los días válidos: segundos con
`col` ≥ `threshold` mg (por defecto, el umbral de Hildebrand et al. (2014) para
ActiGraph en cadera, 7–11 años), unidos si los separan como mucho `tolerance` segundos por
debajo, y con al menos `min_seconds` segundos por encima. Columnas: `inicio`, `fin`,
`medio`, `duracion_s`, `segundos_vigorosos`, `enmo_medio` (mg, en los segundos vigorosos)
y `dosis` (Σ ENMO en esos segundos, mg·s), que pondera el episodio en la media circular.
Con *epochs* de otra duración (p. ej. 5 s, sensibilidad del cap. 12), `tolerance` y
`min_seconds` siguen en segundos y la dosis se multiplica por la duración del *epoch*; con
1 s el resultado es el de siempre.
"""
function vigorous_episodes(ep::AbstractDataFrame, valid_days; threshold::Real = HILDEBRAND_VIGOROUS_HIP,
                           tolerance::Integer = 2, min_seconds::Integer = 3, col::Symbol = :enmo)
    out = DataFrame(inicio = DateTime[], fin = DateTime[], medio = DateTime[], duracion_s = Int[],
                    segundos_vigorosos = Int[], enmo_medio = Float64[], dosis = Float64[])
    v = ep[!, col]; t = ep.start
    above = [(Date(ti) in valid_days) && !nw && isfinite(x) && x >= threshold for (ti, nw, x) in zip(t, ep.nonwear, v)]
    k = 1; N = length(t)
    step = N >= 2 ? Dates.value(t[2] - t[1]) : 1000        # ms por epoch (1000 con epochs de 1 s)
    es = step / 1000; tol = floor(Int, tolerance / es)     # tolerancia en epochs
    while k <= N
        if !above[k]; k += 1; continue; end
        i = k; j = k; gap = 0; m = k + 1
        while m <= N && gap <= tol && Dates.value(t[m] - t[m - 1]) == step
            if above[m]; j = m; gap = 0; else gap += 1; end
            m += 1
        end
        idx = [q for q in i:j if above[q]]
        if length(idx) * es >= min_seconds
            a, b = t[i], t[j] + Millisecond(step)
            dosis = sum(v[idx])
            push!(out, (a, b, a + (b - a) ÷ 2, round(Int, (j - i + 1) * es), round(Int, length(idx) * es),
                        dosis / length(idx), dosis * es))
        end
        k = j + 1
    end
    return out
end

"""
    clock_coverage(ep1, valid_days; nonwear_as = :missing) -> Vector{Int}

Cobertura del reloj: para cada segundo del día (índice 1 = 00:00:00), en cuántos días
válidos se observó con el dispositivo en uso (con `nonwear_as = :rest`, el no uso también
cuenta como observado). Con días completos y sin no uso, vale el número de días en todo el
reloj. Sirve para corregir medias circulares por la cobertura desigual
(`coverage_weights` para episodios; `Spectral.harmonic_timing` para segundos).
"""
function clock_coverage(ep::AbstractDataFrame, valid_days; nonwear_as::Symbol = :missing)
    nobs = zeros(Int, 86400)
    for (t, nw) in zip(ep.start, ep.nonwear)
        (Date(t) in valid_days && (!nw || nonwear_as === :rest)) && (nobs[_sod(t)] += 1)
    end
    return nobs
end
_sod(t) = hour(t) * 3600 + minute(t) * 60 + second(t) + 1

"""
    coverage_weights(episodes, cobertura) -> Vector{Float64}

Peso de cada episodio corregido por la cobertura del reloj: su `dosis` dividida por el
número de días en que se observó el segundo del reloj de su instante medio
(`cobertura = clock_coverage(...)`). Con estos pesos, la media circular es la del perfil
medio de dosis y no depende de qué horas se observaron más días (cap. 11).
"""
coverage_weights(ep::AbstractDataFrame, cob::AbstractVector{<:Integer}) =
    [d / max(cob[_sod(t)], 1) for (d, t) in zip(ep.dosis, ep.medio)]

"""
    hours_to_phase(h, anchor_h) / phase_to_hours(φ)

Hora del reloj (decimal) → fase en radianes respecto a `anchor_h` (0 en el ancla, sentido
de las agujas del reloj = hacia delante en el tiempo), y la inversa, en horas tras el ancla
(0–24).
"""
hours_to_phase(h::Real, anchor_h::Real) = 2π * mod(h - anchor_h, 24) / 24
phase_to_hours(φ::Real) = mod(φ * 24 / (2π), 24)

"""
    circular_mean(φ, w = ones) -> NamedTuple

Media circular ponderada: C = Σ w cos φ / Σ w, S = Σ w sin φ / Σ w, R = √(C² + S²),
θ = atan2(S, C) ∈ [0, 2π). Implementación propia (Fisher, 1993).
"""
function circular_mean(φ::AbstractVector{<:Real}, w::AbstractVector{<:Real} = ones(length(φ)))
    isempty(φ) && return (θ = NaN, R = NaN, C = NaN, S = NaN)
    sw = sum(w)
    C = sum(w .* cos.(φ)) / sw; S = sum(w .* sin.(φ)) / sw
    return (θ = mod(atan(S, C), 2π), R = sqrt(C^2 + S^2), C = C, S = S)
end

"""
    rayleigh(n, R) -> NamedTuple

Test de Rayleigh de uniformidad circular para `n` direcciones con longitud media `R`:
Z = nR² y la aproximación de Zar (1999) p = exp(√(1 + 4n + 4(n² − (nR)²)) − (1 + 2n)).
Supone direcciones independientes: con episodios agrupados en días, el p es optimista.
"""
function rayleigh(n::Integer, R::Real)
    n < 2 && return (Z = NaN, p = NaN)
    Rn = n * R
    p = exp(sqrt(1 + 4n + 4(n^2 - Rn^2)) - (1 + 2n))
    return (Z = n * R^2, p = clamp(p, 0, 1))
end

"""
    episode_phases(episodes, anchor; max_gap = Hour(36)) -> Vector{Float64}

Fase (radianes, 0 en el ancla, hacia delante en el tiempo) del instante medio de cada
episodio. `anchor` puede ser una hora fija del reloj (`Real`, decimal: el mismo ancla
todos los días) o un ancla por noche (`AbstractVector{DateTime}`, p. ej. el punto medio de
cada periodo de reposo nocturno). En el segundo caso, cada episodio se refiere al ancla anterior más
próxima (la de la noche previa) si está a menos de `max_gap`; si no la hay (el primer día
del registro, o tras una noche sin ancla), a la siguiente, con el mismo criterio; si
tampoco, la fase es `NaN` y el episodio no entra en la media.
"""
function episode_phases(ep::AbstractDataFrame, anchor_h::Real; kw...)
    return [hours_to_phase(hour(t) + minute(t) / 60 + second(t) / 3600, anchor_h)   # resolución de 1 s
            for t in ep.medio]
end
function episode_phases(ep::AbstractDataFrame, anclas::AbstractVector{DateTime}; max_gap::Period = Hour(36))
    a = sort(anclas); lim = Dates.value(convert(Millisecond, max_gap))
    return map(ep.medio) do t
        i = searchsortedlast(a, t)
        j = i >= 1 && Dates.value(t - a[i]) < lim ? i :
            i < length(a) && Dates.value(a[i+1] - t) < lim ? i + 1 : 0
        j == 0 ? NaN : 2π * mod(Dates.value(t - a[j]) / 3_600_000, 24) / 24
    end
end

"""
    timing(episodes, anchor; peso = :dosis) -> NamedTuple

*Timing* circular de los episodios respecto al ancla (una hora fija, `anchor_h`, p. ej. el
punto medio de L5, o un ancla por noche; ver `episode_phases`): fase de cada episodio (su
instante medio; los episodios sin ancla se descartan), media circular ponderada por la
dosis (θ, R y sus componentes `R_cos` = R·cos θ y `R_sin` = R·sin θ, los predictores del
modelo metabólico),
θ en horas tras el ancla, la R sin ponderar y el test de Rayleigh. `R2` es la longitud del
segundo armónico (media de e^{2iφ}): un `R2` alto con `R` moderado indica dos grupos de
episodios separados unas 12 h, que la media circular simple resume mal; `θ2_h` es la
orientación de ese eje (y su opuesto, +12 h). `peso` es la columna que pondera cada
episodio: `:dosis`, o una columna con la dosis corregida por cobertura (`coverage_weights`).
"""
function timing(ep::AbstractDataFrame, anchor; peso::Symbol = :dosis)
    φ = episode_phases(ep, anchor)
    ok = .!isnan.(φ); φ = φ[ok]; w = ep[ok, peso]
    cw = circular_mean(φ, w); cu = circular_mean(φ)
    c2 = circular_mean(2 .* φ, w)                        # segundo armónico (bimodalidad)
    ry = rayleigh(length(φ), cu.R)
    return (n = length(φ), θ = cw.θ, θ_h = phase_to_hours(cw.θ), R = cw.R,
            R2 = c2.R, θ2_h = phase_to_hours(c2.θ / 2),
            R_cos = cw.R * cos(cw.θ), R_sin = cw.R * sin(cw.θ),
            R_sin_ponderar = cu.R, rayleigh_Z = ry.Z, rayleigh_p = ry.p, fase = φ)
end

"""
    timing_by_day(episodes, anchor) -> DataFrame

θ (horas tras el ancla), R y número de episodios de cada día, para ver la estabilidad del
*timing* entre días.
"""
function timing_by_day(ep::AbstractDataFrame, anchor)
    out = DataFrame(dia = Date[], episodios = Int[], θ_h = Float64[], R = Float64[])
    for g in groupby(transform(ep, :medio => ByRow(Date) => :dia), :dia; sort = true)
        tm = timing(DataFrame(g), anchor)
        push!(out, (first(g.dia), tm.n, tm.θ_h, tm.R))
    end
    return out
end

"""
    timing_reliability(por_dia; id = :id, k = nothing) -> NamedTuple

Fiabilidad del *timing* de una semana: qué parte de la variación de los predictores
C = R·cos θ y S = R·sin θ es diferencia estable entre niños y qué parte es variación de un
día a otro. `por_dia` tiene una fila por niño y día con episodios y las columnas `id`,
`θ_h` (horas tras el ancla) y `R`, como las de `timing_by_day` apiladas con su ID.

Para cada componente se ajusta un ANOVA de efectos aleatorios de una vía con grupos
desiguales: σ²_entre = (MSB − MSW)/n₀ (truncada en 0), σ²_dentro = MSW y
n₀ = (N − Σnᵢ²/N)/(g − 1). El ICC de un día es σ²_entre/(σ²_entre + σ²_dentro); el de la
media de `k` días, k·ICC/(1 + (k − 1)·ICC) (Spearman–Brown). Por defecto `k` = número medio
de días con episodios por niño. `icc_vector` usa la suma de las varianzas de C y S (la
fiabilidad del vector en el plano). Es la fiabilidad que atenúa el efecto estimado en el
modelo metabólico (error de medida clásico): el coeficiente se multiplica, en promedio,
por algo del orden de esa fiabilidad.

Aproximación: la exposición del modelo es la media ponderada por dosis de todos los
episodios de la semana, no la media sin ponderar de los vectores diarios, y los días sin
episodios no entran. Con un solo niño (`n_ninos = 1`) solo se estima la varianza dentro del
niño; los ICC salen `NaN`.
"""
function timing_reliability(por_dia::AbstractDataFrame; id::Symbol = :id, k = nothing)
    d = por_dia[isfinite.(por_dia.θ_h) .& isfinite.(por_dia.R), :]
    φ = 2π .* d.θ_h ./ 24
    C = d.R .* cos.(φ); S = d.R .* sin.(φ)
    ids = d[!, id]; grupos = unique(ids)
    g = length(grupos); N = length(ids)
    nᵢ = [count(==(u), ids) for u in grupos]
    function componentes(y)
        medias = [sum(y[ids .== u]) / n for (u, n) in zip(grupos, nᵢ)]
        ȳ = sum(y) / N
        ssw = sum((y[i] - medias[findfirst(==(ids[i]), grupos)])^2 for i in 1:N; init = 0.0)
        msw = N > g ? ssw / (N - g) : NaN
        g < 2 && return (entre = NaN, dentro = msw)
        msb = sum(n * (m - ȳ)^2 for (n, m) in zip(nᵢ, medias)) / (g - 1)
        n0 = (N - sum(nᵢ .^ 2) / N) / (g - 1)
        return (entre = max((msb - msw) / n0, 0.0), dentro = msw)
    end
    vc = componentes(C); vs = componentes(S)
    icc(v) = v.entre / (v.entre + v.dentro)
    kk = k === nothing ? N / g : float(k)
    sb(r) = kk * r / (1 + (kk - 1) * r)
    iC = icc(vc); iS = icc(vs)
    iV = (vc.entre + vs.entre) / (vc.entre + vs.entre + vc.dentro + vs.dentro)
    return (n_ninos = g, n_dias = N, k = kk,
            var_entre_C = vc.entre, var_dentro_C = vc.dentro,
            var_entre_S = vs.entre, var_dentro_S = vs.dentro,
            icc_C = iC, icc_S = iS, icc_vector = iV,
            icc_k_C = sb(iC), icc_k_S = sb(iS), icc_k_vector = sb(iV))
end

# Generador xorshift64* mínimo y reproducible para el bootstrap (evita añadir Random como
# dependencia del paquete).
mutable struct _XorShift; s::UInt64; end
function _rand_index(g::_XorShift, n::Integer)
    g.s ⊻= g.s >> 12; g.s ⊻= g.s << 25; g.s ⊻= g.s >> 27
    return Int(((g.s * 0x2545f4914f6cdd1d) >> 11) % UInt64(n)) + 1
end

"""
    bootstrap_days(episodes, anchor; B = 2000, seed = 1, level = 0.95, peso = :dosis) -> NamedTuple

Incertidumbre de θ y R por *bootstrap* de días: se remuestrean los días (con reposición) y
con ellos todos sus episodios, porque los episodios de un mismo día no son independientes.
Devuelve el intervalo de θ (horas tras el ancla) como el arco circular más corto que
contiene `level` de las réplicas (`θ_ic_h = (desde, hasta)`, puede cruzar las 24 h) y su
amplitud en horas, y el intervalo percentil de R. Con pesos corregidos por cobertura, los
pesos se mantienen fijos en las réplicas (aproximación: no se recalcula la cobertura de
cada remuestreo).
"""
function bootstrap_days(ep::AbstractDataFrame, anchor; B::Integer = 2000, seed::Integer = 1,
                        level::Real = 0.95, peso::Symbol = :dosis)
    d = transform(ep, :medio => ByRow(Date) => :dia)
    grupos = [DataFrame(g) for g in groupby(d, :dia)]
    nd = length(grupos)
    nd < 2 && return (θ_ic_h = (NaN, NaN), amplitud_h = NaN, R_ic = (NaN, NaN), B = 0)
    g = _XorShift(UInt64(seed) * 0x9e3779b97f4a7c15 + 1)
    θs = Float64[]; Rs = Float64[]
    for _ in 1:B
        m = reduce(vcat, (grupos[_rand_index(g, nd)] for _ in 1:nd))
        tm = timing(m, anchor; peso = peso)
        push!(θs, tm.θ_h); push!(Rs, tm.R)
    end
    # Arco más corto que contiene `level` de las réplicas: se ordenan y se busca la ventana
    # circular de k puntos consecutivos de menor amplitud.
    sθ = sort(θs); n = length(sθ); k = ceil(Int, level * n)
    best = (Inf, 0)
    for i in 1:n
        j = i + k - 1
        span = j <= n ? sθ[j] - sθ[i] : sθ[j - n] + 24 - sθ[i]
        span < best[1] && (best = (span, i))
    end
    i = best[2]; j = i + k - 1
    desde, hasta = sθ[i], j <= n ? sθ[j] : sθ[j - n]
    sR = sort(Rs); q(p) = sR[clamp(ceil(Int, p * n), 1, n)]
    return (θ_ic_h = (desde, hasta), amplitud_h = best[1], R_ic = (q((1 - level) / 2), q((1 + level) / 2)), B = B)
end

end # module
