"""
    Model

Modelo metabólico del primer objetivo y su potencia (hito H7), probado con cohortes
**sintéticas** mientras no se disponga de la tabla clínica.

Modelo principal (cap. 14): log(HOMA-IR) = β₀ + β₁ R cos θ + β₂ R sin θ + γᵀz + ε, con un
único contraste confirmatorio, H₀: β₁ = β₂ = 0 (F con 2 g.l.). El término circular se
reescribe como A·R·cos(θ − θ₀), A = √(β₁² + β₂²), θ₀ = atan2(β₂, β₁): fase de máximo
HOMA-IR predicho; la más favorable es θ₀ + π.

Sin dependencias nuevas: la regresión usa una QR de Householder propia; el valor p de la F,
una beta incompleta regularizada (fracción continua) y una gamma logarítmica de Lanczos; el
tobit, una función de distribución normal propia; y el azar, un generador xorshift64*
reproducible. Todo se verifica en `test/test_model.jl` frente a valores conocidos.
"""
module Model

using DataFrames

export simulate_cohort, ols, fdist_sf, normcdf, fit_timing_model, bootstrap_theta0,
       tobit_em, power_grid, Rng

# ─── Azar reproducible ────────────────────────────────────────────────────────

"""
    Rng(seed)

Generador xorshift64* reproducible (evita depender de `Random`). `rand(g)` ∈ (0, 1),
`randn(g)` por Box–Muller.
"""
mutable struct Rng; s::UInt64; end
Rng(seed::Integer) = Rng(UInt64(seed) * 0x9e3779b97f4a7c15 + 0x2545f4914f6cdd1d)
function Base.rand(g::Rng)
    g.s ⊻= g.s >> 12; g.s ⊻= g.s << 25; g.s ⊻= g.s >> 27
    return ((g.s * 0x2545f4914f6cdd1d) >> 11) / 2.0^53 + 2.0^-54
end
Base.randn(g::Rng) = sqrt(-2log(rand(g))) * cos(2π * rand(g))
_randidx(g::Rng, n::Integer) = min(n, floor(Int, rand(g) * n) + 1)

# ─── Funciones especiales ─────────────────────────────────────────────────────

# Gamma logarítmica (Lanczos, g = 7, 9 coeficientes; error relativo < 1e-13 para x > 0).
const _LANCZOS = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                  771.32342877765313, -176.61502916214059, 12.507343278686905,
                  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
function _lgamma(x::Real)
    x < 0.5 && return log(π / abs(sin(π * x))) - _lgamma(1 - x)
    x -= 1
    a = _LANCZOS[1]; t = x + 7.5
    for i in 1:8
        a += _LANCZOS[i + 1] / (x + i)
    end
    return 0.5log(2π) + (x + 0.5) * log(t) - t + log(a)
end

# Fracción continua de la beta incompleta (Lentz modificado).
function _betacf(a, b, x; maxit = 300, eps = 3e-14)
    tiny = 1e-300
    qab = a + b; qap = a + 1; qam = a - 1
    c = 1.0; d = 1 - qab * x / qap
    abs(d) < tiny && (d = tiny); d = 1 / d; h = d
    for m in 1:maxit
        m2 = 2m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; abs(d) < tiny && (d = tiny)
        c = 1 + aa / c; abs(c) < tiny && (c = tiny)
        d = 1 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; abs(d) < tiny && (d = tiny)
        c = 1 + aa / c; abs(c) < tiny && (c = tiny)
        d = 1 / d; del = d * c; h *= del
        abs(del - 1) < eps && break
    end
    return h
end

# Beta incompleta regularizada I_x(a, b).
function _betainc(a, b, x)
    (x <= 0) && return 0.0
    (x >= 1) && return 1.0
    lbt = _lgamma(a + b) - _lgamma(a) - _lgamma(b) + a * log(x) + b * log(1 - x)
    return x < (a + 1) / (a + b + 2) ? exp(lbt) * _betacf(a, b, x) / a :
                                       1 - exp(lbt) * _betacf(b, a, 1 - x) / b
end

"""
    fdist_sf(F, d1, d2) -> Float64

P(F_{d1,d2} > F): valor p de un contraste F, = I_{d2/(d2 + d1 F)}(d2/2, d1/2).
"""
fdist_sf(F::Real, d1::Real, d2::Real) = F <= 0 ? 1.0 : _betainc(d2 / 2, d1 / 2, d2 / (d2 + d1 * F))

"""
    normcdf(x) -> Float64

Función de distribución normal estándar (vía erfc de Numerical Recipes, error < 1,2e-7).
"""
function normcdf(x::Real)
    z = abs(x) / sqrt(2); t = 1 / (1 + 0.5z)
    erfc = t * exp(-z^2 - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
           t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
           t * (-0.82215223 + t * 0.17087277)))))))))
    return x >= 0 ? 1 - erfc / 2 : erfc / 2
end
_normpdf(x) = exp(-x^2 / 2) / sqrt(2π)

# ─── Mínimos cuadrados ────────────────────────────────────────────────────────

# Producto matriz-vector propio (el operador * de matrices vive en LinearAlgebra).
_mul(X::AbstractMatrix, b::AbstractVector) = [sum(X[i, j] * b[j] for j in axes(X, 2)) for i in axes(X, 1)]

"""
    ols(X, y) -> NamedTuple

Mínimos cuadrados ordinarios por QR de Householder (estable sin formar XᵀX). Devuelve
`β`, la suma de cuadrados residual `rss`, `n`, `p` y los residuos.
"""
function ols(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real})
    A = Matrix{Float64}(X); b = Vector{Float64}(y)
    n, p = size(A)
    n > p || throw(ArgumentError("hace falta n > p"))
    for k in 1:p
        v = A[k:n, k]
        α = -sign(v[1] == 0 ? 1.0 : v[1]) * sqrt(sum(abs2, v))
        v[1] -= α
        nv = sum(abs2, v)
        nv == 0 && continue
        for j in k:p
            A[k:n, j] .-= (2 * sum(v .* A[k:n, j]) / nv) .* v
        end
        b[k:n] .-= (2 * sum(v .* b[k:n]) / nv) .* v
    end
    β = zeros(p)
    for i in p:-1:1                                   # sustitución regresiva con R
        β[i] = (b[i] - sum(A[i, j] * β[j] for j in i+1:p; init = 0.0)) / A[i, i]
    end
    r = Vector{Float64}(y) .- _mul(X, β)
    return (β = β, rss = sum(abs2, r), n = n, p = p, residuos = r)
end

# ─── Cohortes sintéticas ──────────────────────────────────────────────────────

"""
    simulate_cohort(n; A = 0.2, θ0_h = 18, θ_media_h = 11.3, θ_de_h = 1.5, R_rango = (0.4, 0.8),
                    error_θ_h = 0, error_R = 0, error_cs = 0, homa_mediana = 0.65, de_log = 0.62,
                    suelo = nothing, seed = 1) -> DataFrame

Cohorte sintética de `n` niños para probar el modelo y estimar la potencia. **Todos los
valores por defecto son supuestos de simulación**, documentados en el cap. 14:
- covariables: edad ~ U(3, 6) años, sexo ~ Bernoulli(0,5), IMC z ~ N(0,14; 1,27) y
  AFMV ~ N(72,1; 28,4) min/día (descriptivos de COHORTE_A);
- timing verdadero: θ ~ N(`θ_media_h`, `θ_de_h`) horas tras el punto medio de L5 y
  R ~ U(`R_rango`) (en el rango de los dos registros de H6);
- log(HOMA-IR) = log(`homa_mediana`) + 0,15·IMC z − 0,05·AFMV z + A·R·cos(θ − θ₀) + ε,
  ε ~ N(0, `de_log`) (mediana y dispersión del HOMA-IR preescolar de MELI-POP);
- θ y R observados = verdaderos + error de medida (`error_θ_h`, `error_R`), para imitar la
  incertidumbre de una semana de registro;
- con `error_cs` > 0, además, error clásico en el plano: C = R·cos θ y S = R·sin θ
  observados = verdaderos + N(0, `error_cs`) en cada componente (la R observada puede
  pasar de 1: no se recorta, para que el error siga siendo clásico). Así la fiabilidad de
  la exposición se fija directamente (ver `Circadian.timing_reliability`);
- con `suelo`, el HOMA-IR observado no baja de ese valor (censura por el límite de
  detección de la insulina); `censurado` lo marca.
"""
function simulate_cohort(n::Integer; A::Real = 0.2, θ0_h::Real = 18.0, θ_media_h::Real = 11.3,
                         θ_de_h::Real = 1.5, R_rango = (0.4, 0.8), error_θ_h::Real = 0.0,
                         error_R::Real = 0.0, error_cs::Real = 0.0, homa_mediana::Real = 0.65,
                         de_log::Real = 0.62,
                         suelo = nothing, seed::Integer = 1)
    g = Rng(seed)
    edad = [3 + 3rand(g) for _ in 1:n]
    sexo = [rand(g) < 0.5 ? 1.0 : 0.0 for _ in 1:n]
    imc_z = [0.14 + 1.27randn(g) for _ in 1:n]
    afmv = [72.1 + 28.4randn(g) for _ in 1:n]
    afmv_z = (afmv .- 72.1) ./ 28.4
    θh = [mod(θ_media_h + θ_de_h * randn(g), 24) for _ in 1:n]
    R = [R_rango[1] + (R_rango[2] - R_rango[1]) * rand(g) for _ in 1:n]
    θ0 = 2π * θ0_h / 24
    loghoma = [log(homa_mediana) + 0.15imc_z[i] - 0.05afmv_z[i] + A * R[i] * cos(2π * θh[i] / 24 - θ0) +
               de_log * randn(g) for i in 1:n]
    θobs = [mod(θh[i] + error_θ_h * randn(g), 24) for i in 1:n]
    Robs = [clamp(R[i] + error_R * randn(g), 0, 1) for i in 1:n]
    if error_cs > 0          # solo entonces se consumen más números: con 0 la cohorte no cambia
        for i in 1:n
            c = Robs[i] * cos(2π * θobs[i] / 24) + error_cs * randn(g)
            s_ = Robs[i] * sin(2π * θobs[i] / 24) + error_cs * randn(g)
            Robs[i] = hypot(c, s_); θobs[i] = mod(24 * atan(s_, c) / (2π), 24)
        end
    end
    homa = exp.(loghoma)
    cens = suelo === nothing ? falses(n) : homa .<= suelo
    suelo === nothing || (homa = max.(homa, suelo))
    return DataFrame(edad = edad, sexo = sexo, imc_z = imc_z, afmv = afmv, θ_h = θobs, R = Robs,
                     θ_verdadero_h = θh, R_verdadero = R, homa = homa, censurado = cens)
end

# ─── Modelo principal ─────────────────────────────────────────────────────────

const COVARIABLES = (:edad, :sexo, :imc_z, :afmv)

function _design(df::AbstractDataFrame, covs)
    n = nrow(df)
    φ = 2π .* df.θ_h ./ 24
    Xc = hcat(df.R .* cos.(φ), df.R .* sin.(φ))
    Z = isempty(covs) ? zeros(n, 0) : hcat((Float64.(df[!, c]) for c in covs)...)
    return hcat(ones(n), Xc, Z), hcat(ones(n), Z)
end

"""
    fit_timing_model(df; covs = (:edad, :sexo, :imc_z, :afmv), y = log.(df.homa)) -> NamedTuple

Ajusta el modelo principal y el reducido (sin el término circular) y devuelve los
coeficientes circulares `β_cos`, `β_sin`, `A`, `θ0_h` (fase de máximo HOMA-IR predicho, en
horas tras el ancla), `fase_favorable_h` (= θ₀ + 12 h), el contraste `F` con 2 y n − p g.l.,
su valor `p` y `n`.
"""
function fit_timing_model(df::AbstractDataFrame; covs = COVARIABLES, y = log.(df.homa))
    Xf, Xr = _design(df, covs)
    f = ols(Xf, y); r = ols(Xr, y)
    d2 = f.n - f.p
    F = ((r.rss - f.rss) / 2) / (f.rss / d2)
    b1, b2 = f.β[2], f.β[3]
    θ0 = mod(atan(b2, b1), 2π)
    return (β_cos = b1, β_sin = b2, A = sqrt(b1^2 + b2^2), θ0_h = 24θ0 / (2π),
            fase_favorable_h = mod(24θ0 / (2π) + 12, 24), F = F, p = fdist_sf(F, 2, d2), n = f.n)
end

"""
    bootstrap_theta0(df; B = 1000, seed = 1, level = 0.95, kw...) -> NamedTuple

Intervalo de θ₀ por *bootstrap* de niños: arco circular más corto que contiene `level` de
las réplicas (puede cruzar las 24 h), su amplitud en horas, y el intervalo percentil de A.
Si el contraste conjunto no es significativo, el arco abarca casi toda la circunferencia y
θ₀ no se interpreta (cap. 14).
"""
function bootstrap_theta0(df::AbstractDataFrame; B::Integer = 1000, seed::Integer = 1,
                          level::Real = 0.95, kw...)
    g = Rng(seed); n = nrow(df)
    θs = Float64[]; As = Float64[]
    for _ in 1:B
        idx = [_randidx(g, n) for _ in 1:n]
        m = fit_timing_model(df[idx, :]; kw...)
        push!(θs, m.θ0_h); push!(As, m.A)
    end
    s = sort(θs); k = ceil(Int, level * B)
    best = (Inf, 1)
    for i in 1:B
        j = i + k - 1
        span = j <= B ? s[j] - s[i] : s[j - B] + 24 - s[i]
        span < best[1] && (best = (span, i))
    end
    i = best[2]; j = i + k - 1
    sa = sort(As); q(p) = sa[clamp(ceil(Int, p * B), 1, B)]
    return (θ0_ic_h = (s[i], j <= B ? s[j] : s[j - B]), amplitud_h = best[1],
            A_ic = (q((1 - level) / 2), q((1 + level) / 2)))
end

"""
    tobit_em(X, y, censored; iters = 200, tol = 1e-8) -> NamedTuple

Regresión lineal con censura por la izquierda (tobit) por el algoritmo EM: en cada paso,
los valores censurados se sustituyen por su esperanza condicionada a estar por debajo del
límite (`y` ya contiene el límite en esos casos) y se reestima por mínimos cuadrados, con
la corrección de la varianza. Devuelve `β` y `σ`.
"""
function tobit_em(X::AbstractMatrix, y::AbstractVector, censored::AbstractVector{Bool};
                  iters::Integer = 200, tol::Real = 1e-8)
    f = ols(X, y); β = f.β; σ = sqrt(f.rss / (f.n - f.p))
    c = Vector{Float64}(y)                         # límites (iguales a y en los censurados)
    for _ in 1:iters
        μ = _mul(X, β)
        ey = copy(c); vy = zeros(length(c))
        for i in findall(censored)
            a = (c[i] - μ[i]) / σ
            λ = _normpdf(a) / max(normcdf(a), 1e-300)
            ey[i] = μ[i] - σ * λ                   # E[y | y < c]
            vy[i] = σ^2 * (1 - a * λ - λ^2)         # Var[y | y < c]
        end
        fn = ols(X, ey)
        σn = sqrt((sum(abs2, ey .- _mul(X, fn.β)) + sum(vy)) / length(c))
        done = maximum(abs.(fn.β .- β)) < tol && abs(σn - σ) < tol
        β, σ = fn.β, σn
        done && break
    end
    return (β = β, σ = σ)
end

"""
    power_grid(ns, As; errores = (0.0,), nsim = 500, α = 0.05, kw...) -> DataFrame

Potencia del contraste conjunto por simulación: para cada tamaño `n`, efecto `A` y error
de medida de θ (`error_θ_h`), proporción de las `nsim` cohortes sintéticas con p < α. Con
A = 0 da el error de tipo I.
"""
function power_grid(ns, As; errores = (0.0,), nsim::Integer = 500, α::Real = 0.05, kw...)
    out = DataFrame(n = Int[], A = Float64[], error_θ_h = Float64[], potencia = Float64[])
    for n in ns, A in As, e in errores
        hits = 0
        for s in 1:nsim
            df = simulate_cohort(n; A = A, error_θ_h = e, error_R = e > 0 ? 0.1 : 0.0,
                                 seed = 1_000_003 * s + 7919n + round(Int, 1000A) + round(Int, 100e), kw...)
            hits += fit_timing_model(df).p < α
        end
        push!(out, (n, A, e, hits / nsim))
    end
    return out
end

end # module
