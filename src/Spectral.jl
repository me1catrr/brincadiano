"""
    Spectral

Análisis espectral del ritmo de reposo-actividad (escala circadiana y ultradiana: *epochs*
de minutos a lo largo de días). Todo implementado a mano, sin rellenar huecos:

- [`gls_periodogram`](@ref): periodograma de Lomb–Scargle generalizado (Zechmeister y
  Kürster, 2009): en cada frecuencia, ajuste por mínimos cuadrados de media + sinusoide
  sobre los instantes con dato. Admite huecos y no uso sin imputarlos.
- [`spectral_window`](@ref): ventana espectral del patrón de muestreo (qué picos introduce
  la propia forma de los huecos, p. ej. las noches sin dispositivo).
- [`chi2_periodogram`](@ref): periodograma χ² de Sokolove y Bushell (1978).
- [`harmonic_is`](@ref): IS descompuesto en los armónicos de 24 h (Parseval).
- [`harmonic_timing`](@ref): momentos trigonométricos de una dosis de actividad respecto a
  un ancla; el primero es el θ y R del *timing* (cap. 11) sin definir episodios.

En la escala del ritmo, las frecuencias van en ciclos por día (cpd) y los tiempos en horas.
En la escala del movimiento (señal a la frecuencia de muestreo, 0,25–15 Hz a 30 Hz):

- [`window_spectra`](@ref): densidad espectral de la norma del vector en ventanas cortas
  alineadas al reloj, solo si están completas (sin huecos), y sus rasgos: frecuencia
  dominante (cadencia), potencia por bandas, entropía espectral.
"""
module Spectral

using DataFrames, Dates
using FFTW: rfft
using DSP: hanning
using ..Circadian: hours_to_phase, phase_to_hours, _hourly_series, clock_coverage, _sod
using ..Model: normcdf

export freq_grid, gls_periodogram, spectral_window, chi2_quantile, chi2_periodogram,
       harmonic_is, harmonic_timing, clock_coverage, window_spectra, BANDAS_MOVIMIENTO

"""
    freq_grid(T_h; fmin = 0.25, fmax = 24, ofac = 10) -> Vector{Float64}

Rejilla de frecuencias (cpd) para un registro de `T_h` horas: paso 1/(T·ofac), la
resolución natural 1/T sobremuestreada `ofac` veces.
"""
function freq_grid(T_h::Real; fmin::Real = 0.25, fmax::Real = 24, ofac::Real = 10)
    df = 24 / (T_h * ofac)
    return collect(fmin:df:fmax)
end

# Resuelve un sistema 3×3 por eliminación con pivoteo parcial (sin LinearAlgebra).
function _solve3(M::NTuple{9,Float64}, r::NTuple{3,Float64})
    A = [M[1] M[2] M[3] r[1]; M[4] M[5] M[6] r[2]; M[7] M[8] M[9] r[3]]
    for k in 1:3
        p = k - 1 + argmax(abs.(A[k:3, k]))
        p != k && (A[[k, p], :] = A[[p, k], :])
        abs(A[k, k]) < 1e-300 && return (NaN, NaN, NaN)
        for i in k+1:3
            f = A[i, k] / A[k, k]
            A[i, k:4] .-= f .* A[k, k:4]
        end
    end
    x3 = A[3, 4] / A[3, 3]
    x2 = (A[2, 4] - A[2, 3] * x3) / A[2, 2]
    x1 = (A[1, 4] - A[1, 2] * x2 - A[1, 3] * x3) / A[1, 1]
    return (x1, x2, x3)
end

"""
    gls_periodogram(t_h, y, freqs) -> DataFrame

Lomb–Scargle generalizado con media flotante: para cada frecuencia `f` (cpd) ajusta
y ≈ c + a cos ωt + b sin ωt (ω = 2πf/24, t en horas) por mínimos cuadrados sobre los
puntos con `y` finito. Columnas:

- `frecuencia` (cpd) y `periodo_h`;
- `potencia` = 1 − RSS/TSS: fracción de la varianza que explica la sinusoide (0–1);
- `amplitud` = √(a² + b²), en las unidades de `y`;
- `acrofase_h`: instante del máximo, en horas desde t = 0, dentro de un periodo.

Con muestreo regular y completo coincide con el periodograma clásico de Lomb–Scargle.
"""
function gls_periodogram(t_h::AbstractVector{<:Real}, y::AbstractVector{<:Real}, freqs::AbstractVector{<:Real})
    ok = isfinite.(y) .& isfinite.(t_h)
    t = Float64.(t_h[ok]); v = Float64.(y[ok]); N = length(v)
    out = DataFrame(frecuencia = Float64[], periodo_h = Float64[], potencia = Float64[],
                    amplitud = Float64[], acrofase_h = Float64[])
    N < 4 && return out
    Sy = sum(v); Syy = sum(abs2, v); TSS = Syy - Sy^2 / N
    for f in freqs
        ω = 2π * f / 24
        Sc = Ss = Scc = Sss = Scs = Syc = Sys = 0.0
        @inbounds for i in 1:N
            s, c = sincos(ω * t[i])
            Sc += c; Ss += s; Scc += c * c; Sss += s * s; Scs += c * s
            Syc += v[i] * c; Sys += v[i] * s
        end
        β = _solve3((Float64(N), Sc, Ss, Sc, Scc, Scs, Ss, Scs, Sss), (Sy, Syc, Sys))
        RSS = Syy - (β[1] * Sy + β[2] * Syc + β[3] * Sys)
        a, b = β[2], β[3]
        P = 24 / f
        push!(out, (f, P, TSS > 0 ? 1 - RSS / TSS : NaN, hypot(a, b), mod(atan(b, a) / ω, P)))
    end
    return out
end

"""
    spectral_window(t_h, freqs) -> Vector{Float64}

Ventana espectral del muestreo: W(f) = |Σ exp(iωtⱼ)|² / N², con los instantes `t_h` (horas)
en que hay dato. Vale 1 en f = 0; si el patrón de huecos se repite cada 24 h (p. ej. el
dispositivo se retira cada noche), W tiene picos en 1, 2, 3… cpd, y cualquier pico del
periodograma en esas frecuencias puede venir en parte de los huecos y no de la señal.
"""
function spectral_window(t_h::AbstractVector{<:Real}, freqs::AbstractVector{<:Real})
    N = length(t_h)
    return [begin
                ω = 2π * f / 24; C = S = 0.0
                for t in t_h
                    s, c = sincos(ω * t); C += c; S += s
                end
                (C^2 + S^2) / N^2
            end for f in freqs]
end

# Cuantil de la normal estándar por bisección sobre `Model.normcdf`.
function _norm_quantile(p::Real)
    lo, hi = -10.0, 10.0
    for _ in 1:100
        m = (lo + hi) / 2
        normcdf(m) < p ? (lo = m) : (hi = m)
    end
    return (lo + hi) / 2
end

"""
    chi2_quantile(k, p) -> Float64

Cuantil `p` de la χ² con `k` grados de libertad por la aproximación de Wilson–Hilferty
(error relativo < 1 % para k ≥ 5; basta para el umbral del periodograma).
"""
function chi2_quantile(k::Real, p::Real)
    z = _norm_quantile(p); h = 2 / (9k)
    return k * (1 - h + z * sqrt(h))^3
end

"""
    chi2_periodogram(y, bin_h; periods_h = 2:bin_h:36, α = 0.01) -> DataFrame

Periodograma χ² de Sokolove y Bushell (1978) para una serie equiespaciada `y` (un valor
cada `bin_h` horas; `NaN` = ausente). Para cada periodo de P intervalos se pliega la serie
y se calcula Qp = N Σₕ nₕ (M̄ₕ − M̄)² / Σᵢ (xᵢ − M̄)², con nₕ los valores disponibles en la
posición h del ciclo (con datos completos, la fórmula original). Bajo ausencia de ritmo,
Qp ~ χ² con (posiciones con dato − 1) grados de libertad; `Qp_crit` es su cuantil 1 − α.
Los *epochs* consecutivos están autocorrelacionados, así que el umbral es optimista: sirve
para ordenar periodos, no como contraste. Relación con IS: con medias horarias y P = 24,
Qp = N · IS exactamente (también con horas ausentes, con la IS generalizada de `is_iv`).
"""
function chi2_periodogram(y::AbstractVector{<:Real}, bin_h::Real; periods_h = 2:bin_h:36, α::Real = 0.01)
    ok = isfinite.(y)
    N = count(ok); M = sum(y[ok]) / max(N, 1)
    tot = sum((y[ok] .- M) .^ 2)
    out = DataFrame(periodo_h = Float64[], Qp = Float64[], gl = Int[], Qp_crit = Float64[], Qp_rel = Float64[])
    for P in unique(round.(Int, periods_h ./ bin_h))
        P < 2 && continue
        s = zeros(P); c = zeros(Int, P)
        for (i, v) in enumerate(y)
            ok[i] || continue
            h = mod(i - 1, P) + 1
            s[h] += v; c[h] += 1
        end
        Qp = N * sum(c[h] * (s[h] / c[h] - M)^2 for h in 1:P if c[h] > 0) / tot
        gl = count(>(0), c) - 1
        crit = chi2_quantile(gl, 1 - α)
        push!(out, (P * bin_h, Qp, gl, crit, Qp / crit))
    end
    return out
end

"""
    harmonic_is(ep60, valid_days; col = :enmo_idle0, nonwear_as = :missing, min_minutes = 30)
        -> NamedTuple

IS en los armónicos de 24 h. Con las medias horarias de los días válidos (las mismas que
`Circadian.is_iv`), el perfil medio d(h) = M̄ₕ − M̄ se descompone por Fourier en sus 24
puntos, y por Parseval la varianza del perfil se reparte entre los armónicos k = 1…12
(periodos de 24, 12, 8, 6… h). IS_k = n̄ · (varianza del armónico k) / Σ(xᵢ − M̄)², con n̄
el número medio de días por hora. Con datos completos, Σₖ IS_k = IS exactamente; con horas
ausentes la suma se aparta un poco (`IS` es el valor de `is_iv`). Devuelve `IS`,
`IS_armonicos` (la suma), `media` (M̄), `perfil` (las 24 medias M̄ₕ, de 0 a 23 h) y
`armonicos`: `k`, `periodo_h`, `IS_k`, `amplitud` (unidades de `col`) y `acrofase_h` (hora
del reloj del primer máximo del armónico). Perfil reconstruido con los K primeros:
M̄ + Σₖ amplitudₖ cos(2πk(h − acrofaseₖ)/24), con h en el centro de cada hora.
"""
function harmonic_is(ep::AbstractDataFrame, valid_days; col::Symbol = :enmo_idle0, nonwear_as::Symbol = :missing,
                     min_minutes::Integer = 30)
    days = sort(collect(valid_days))
    vacio = DataFrame(k = Int[], periodo_h = Float64[], IS_k = Float64[], amplitud = Float64[], acrofase_h = Float64[])
    isempty(days) && return (IS = NaN, IS_armonicos = NaN, media = NaN, perfil = fill(NaN, 24), armonicos = vacio)
    x = _hourly_series(ep, days; col = col, nonwear_as = nonwear_as, min_minutes = min_minutes)
    ok = isfinite.(x); n = count(ok)
    n < 2 && return (IS = NaN, IS_armonicos = NaN, media = NaN, perfil = fill(NaN, 24), armonicos = vacio)
    xm = sum(x[ok]) / n; tot = sum((x[ok] .- xm) .^ 2)
    nh = [count(isfinite, x[h:24:end]) for h in 1:24]
    Mh = [nh[h] > 0 ? sum(filter(isfinite, x[h:24:end])) / nh[h] : NaN for h in 1:24]
    IS = sum(nh[h] * (Mh[h] - xm)^2 for h in 1:24 if nh[h] > 0) / tot
    any(==(0), nh) && return (IS = IS, IS_armonicos = NaN, media = xm, perfil = Mh, armonicos = vacio)   # perfil incompleto
    d = Mh .- xm; nbar = sum(nh) / 24
    arm = copy(vacio)
    for k in 1:12
        P = sum(d[h+1] * cis(-2π * k * h / 24) for h in 0:23)
        var_k = (k < 12 ? 2 : 1) * abs2(P) / 24              # Σₕ d² = Σₖ var_k (Parseval)
        A = (k < 12 ? 2 : 1) * abs(P) / 24
        ψ = -angle(P)                                       # d(h) ≈ A cos(2πkh/24 − ψ)
        acro = mod(ψ * 24 / (2π * k) + 0.5, 24 / k)          # +0,5: la media horaria representa el centro de la hora
        push!(arm, (k, 24 / k, nbar * var_k / tot, A, acro))
    end
    return (IS = IS, IS_armonicos = sum(arm.IS_k), media = xm, perfil = Mh, armonicos = arm)
end


"""
    harmonic_timing(ep1, valid_days, anchor_h; threshold = 0, col = :enmo, K = 4,
                    cobertura = false, nonwear_as = :missing) -> NamedTuple

Momentos trigonométricos de la «dosis» de actividad respecto al ancla: cada segundo de los
días válidos, en uso y con `col` ≥ `threshold` (m*g*), pesa su valor; con fases φ respecto
al ancla, Rₖ e^{iθₖ} = Σ w e^{ikφ} / Σ w, k = 1…K. Es el coeficiente k de Fourier de la
dosis plegada en 24 h, normalizado por la media: R₁ y θ₁ son el θ y R del *timing*
(cap. 11) sin pasar por episodios (sin pausas ni duración mínima), y R₂ es su segundo
armónico. Con `threshold = 0`, θ₁ es la acrofase de toda la actividad respecto al ancla.
`θ_h` es la posición del primer máximo del armónico k, en horas tras el ancla (< 24/k).

La media circular supone que todas las horas del reloj se observaron el mismo tiempo. Si el
primer o el último día válido están incompletos, o hay no uso, unas horas pesan más que
otras. Con `cobertura = true`, cada segundo se divide por el número de días en que se
observó ese segundo del reloj (en uso; con `nonwear_as = :rest`, también el no uso cuenta
como observado): son los momentos del perfil medio de dosis, no de la dosis acumulada.
"""
function harmonic_timing(ep::AbstractDataFrame, valid_days, anchor_h::Real; threshold::Real = 0,
                         col::Symbol = :enmo, K::Integer = 4, cobertura::Bool = false,
                         nonwear_as::Symbol = :missing)
    nobs = cobertura ? clock_coverage(ep, valid_days; nonwear_as = nonwear_as) : ones(Int, 86400)
    Z = zeros(ComplexF64, K); W = 0.0; n = 0
    for (t, nw, v) in zip(ep.start, ep.nonwear, ep[!, col])
        (Date(t) in valid_days && !nw && isfinite(v) && v >= threshold) || continue
        w = v / nobs[_sod(t)]
        φ = hours_to_phase(hour(t) + minute(t) / 60 + second(t) / 3600, anchor_h)
        for k in 1:K
            Z[k] += w * cis(k * φ)
        end
        W += w; n += 1
    end
    arm = DataFrame(k = 1:K, R = W > 0 ? abs.(Z) ./ W : fill(NaN, K),
                    θ_h = [W > 0 ? phase_to_hours(mod(angle(Z[k]), 2π)) / k : NaN for k in 1:K])
    return (segundos = n, dosis = W, armonicos = arm)
end

# ─── Escala del movimiento ────────────────────────────────────────────────────

"""
Bandas de la PSD del movimiento (Hz): movimiento voluntario (marcha y carrera, con su
fundamental hasta unos 4 Hz en niños), armónicos e impactos, y vibración o ruido.
"""
const BANDAS_MOVIMIENTO = ((0.25, 5.0), (5.0, 10.0), (10.0, 15.0))

"""
    window_spectra(acc, fs, window = Second(5); origin = floor(first(acc.time), Day),
                   offset = (0, 0, 0), scale = (1, 1, 1), nfft = nextpow(2, muestras),
                   bandas = BANDAS_MOVIMIENTO, rango_dominante = (0.5, 5.0), psd = false)
        -> NamedTuple

Espectro del movimiento en ventanas consecutivas de duración `window`, alineadas al reloj
desde `origin` (como `Metrics.epoch_metrics`). Solo entran las ventanas **completas** (todas
sus muestras; ninguna con huecos de *idle sleep* ni de otro tipo): no se rellena nada. En
cada una, la norma del vector calibrada (`offset`, `scale` de `autocalibrate`) menos su
media, con ventana de Hann y relleno de ceros hasta `nfft`, da la densidad espectral
unilateral en m*g*²/Hz, normalizada para que su integral sea la varianza de la ventana.

Devuelve `rasgos` (un `DataFrame` por ventana completa: `start`; `potencia`, la varianza en
0,25 Hz–Nyquist, m*g*²; `f_dom`, la frecuencia del máximo dentro de `rango_dominante`, la
cadencia si hay locomoción; `frac_1`…`frac_3`, la fracción de potencia en cada banda de
`bandas`; `entropia`, la entropía espectral normalizada (0: un solo pico; 1: ruido
blanco)), `frecuencias` y, con `psd = true`, la matriz de PSD (ventanas × frecuencias,
`Float32`).
"""
function window_spectra(acc::AbstractDataFrame, fs::Integer, window::Period = Second(5);
                        origin::DateTime = floor(first(acc.time), Day), offset = (0.0, 0.0, 0.0),
                        scale = (1.0, 1.0, 1.0), nfft::Union{Nothing,Integer} = nothing,
                        bandas = BANDAS_MOVIMIENTO, rango_dominante = (0.5, 5.0), psd::Bool = false)
    w = Dates.value(convert(Millisecond, window))
    L = round(Int, fs * w / 1000)
    nf = something(nfft, nextpow(2, L))
    win = hanning(L); U = sum(abs2, win)
    fr = collect(0:nf÷2) .* (fs / nf); df = fs / nf
    util = fr .>= 0.25
    kb = [findall(f -> a <= f < b || (b == fs / 2 && f == b), fr) for (a, b) in bandas]
    kd = findall(f -> rango_dominante[1] <= f <= rango_dominante[2], fr)
    rasgos = DataFrame(start = DateTime[], potencia = Float64[], f_dom = Float64[],
                       frac_1 = Float64[], frac_2 = Float64[], frac_3 = Float64[], entropia = Float64[])
    P = Vector{Float32}[]
    _ventanas!(rasgos, P, acc.time, acc.x, acc.y, acc.z, origin, w, L, nf, fs, win, U, fr, df, util, kb, kd,
               offset, scale, psd)
    M = psd ? (isempty(P) ? zeros(Float32, 0, length(fr)) : permutedims(reduce(hcat, P))) : nothing
    return (rasgos = rasgos, frecuencias = fr, psd = M)
end

# Bucle de `window_spectra` (función barrera: columnas con tipo concreto).
function _ventanas!(rasgos, P, t, x, y, z, origin, w, L, nf, fs, win, U, fr, df, util, kb, kd, offset, scale, psd)
    buf = zeros(L); xpad = zeros(nf)
    N = length(t); i = 1
    while i <= N
        k = fld(Dates.value(t[i] - origin), w)
        j = i                                              # muestras i:j en la misma ventana
        while j < N && fld(Dates.value(t[j+1] - origin), w) == k
            j += 1
        end
        if k >= 0 && j - i + 1 == L                        # ventana completa
            m = 0.0
            @inbounds for q in 1:L
                ii = i + q - 1
                buf[q] = sqrt(((Float64(x[ii]) + offset[1]) * scale[1])^2 +
                              ((Float64(y[ii]) + offset[2]) * scale[2])^2 +
                              ((Float64(z[ii]) + offset[3]) * scale[3])^2)
                m += buf[q]
            end
            m /= L
            @inbounds for q in 1:L; xpad[q] = (buf[q] - m) * win[q]; end
            X = rfft(xpad)
            p = abs2.(X) ./ (fs * U) .* 1e6               # m𝑔²/Hz (la norma va en 𝑔)
            p[2:end-(iseven(nf) ? 1 : 0)] .*= 2           # unilateral
            pu = view(p, util); su = sum(pu); tot = su * df
            H = 0.0
            @inbounds for q in pu
                q > 0 && (H -= (q / su) * log(q / su))
            end
            push!(rasgos, (origin + Millisecond(k * w), tot, fr[kd[argmax(view(p, kd))]],
                           (sum(view(p, kk)) * df / tot for kk in kb)..., H / log(length(pu))))
            psd && push!(P, Float32.(p))
        end
        i = j + 1
    end
    return nothing
end

end # module
