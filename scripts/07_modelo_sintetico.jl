# 07 · Modelo metabólico con cohortes sintéticas: potencia, ejemplo y censura.
#
# Qué hace
#   Simula cohortes con Brincadiano.Model y estima la potencia del contraste de 2 g.l.
#   sobre R·cos θ y R·sin θ, un ajuste de ejemplo y el efecto de censurar el HOMA-IR
#   (mínimos cuadrados frente a tobit). No abre ningún .gt3x ni data/interim/: los
#   supuestos están en el capítulo 15, y por eso las tablas y figuras sí se versionan.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/07_modelo_sintetico.jl              # 500 simulaciones por casilla
#   julia --project=scripts scripts/07_modelo_sintetico.jl --nsim 500
#   julia --project=scripts scripts/07_modelo_sintetico.jl --solo-fiabilidad   # solo la sección 4
#
# Depende de
#   De ningún script del pipeline. Solo de Brincadiano.Model y de la extensión de
#   figuras (CairoMakie).
#
# Lo usan
#   10, que lee outputs/metabolic/tables/ejemplo_modelo.csv (escenario «con efecto»)
#   para colocar los dos registros reales en el plano circular. El resto de las
#   salidas entra directo en el capítulo 15; no pasa por `make sync-assets`.
#
# Entrada
#   Ningún archivo. Parámetro --nsim (entero; por defecto 500).
#
# Salida versionada (material editorial del capítulo 15)
#   report/assets/14_salud_metabolica/report/tables/potencia.tex
#   report/assets/14_salud_metabolica/report/tables/ejemplo.tex
#   report/assets/14_salud_metabolica/report/tables/censura.tex
#   report/assets/14_salud_metabolica/report/figures/potencia.pdf
#   report/assets/14_salud_metabolica/report/figures/efecto-circular.pdf
#   report/assets/14_salud_metabolica/report/figures/fiabilidad.pdf
# Copia de trabajo (no hace falta versionarla; el 10 lee la segunda)
#   outputs/metabolic/tables/potencia.csv
#   outputs/metabolic/tables/ejemplo_modelo.csv
#   outputs/metabolic/tables/fiabilidad.csv, fiabilidad_semana.csv
#
# Rejilla
#   n ∈ {100, 200, 300, 420, 600, 900}. 420 ≈ cruce estimado de COHORTE_A (anexo de fichas).
#   A ∈ {0, 0,1, 0,2, 0,3, 0,4}. Con R ≈ 0,6, la fase mejor y la peor difieren en un
#   factor exp(2·A·0,6) del HOMA-IR: 13, 27, 43 y 62 %.
#   Error de medida de θ: 0 o 1,1 h (el intervalo de una semana del cap. 11); de R: 0,1.
#   Dispersión de θ entre niños: 1,5 h o 4 h.
#   Censura: suelo de HOMA-IR 0,39 (insulina 2 µU/mL × glucosa 4,4 mmol/L / 22,5).
#
# Fiabilidad (sección 4): error clásico en C = R·cos θ y S = R·sin θ (`error_cs`) que fija
#   la fiabilidad del vector λ = var(C,S verdaderos) / var(C,S observados) en
#   {1; 0,9; …; 0,3}; n ∈ {420, 900}, A = 0,3, dispersión de θ de 1,5 o 4 h.
#
# Consola: potencias y el n que alcanza el 80 %. No hay identificadores de participantes.

using Brincadiano.Model
using Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REP = joinpath(ROOT, "report", "assets", "14_salud_metabolica", "report")
const OUT = joinpath(ROOT, "outputs", "metabolic", "tables")
const NSIM = (i = findfirst(==("--nsim"), ARGS); i === nothing ? 500 : parse(Int, ARGS[i + 1]))
const SOLO_FIABILIDAD = "--solo-fiabilidad" in ARGS

const AZUL = RGBf(13 / 255, 45 / 255, 75 / 255); const TURQUESA = RGBf(20 / 255, 150 / 255, 145 / 255)
const AMBAR = RGBf(224 / 255, 121 / 255, 34 / 255); const MORADO = RGBf(116 / 255, 90 / 255, 170 / 255)

# Coma decimal (convención del informe).
fmt(x; d = 2) = replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
pct(x) = @sprintf("%.0f", 100x)

mkpath(joinpath(REP, "tables")); mkpath(joinpath(REP, "figures")); mkpath(OUT)

set_theme!(tema_informe())
coma(x) = replace(string(x), "." => ",")        # coma decimal en los textos de las figuras

if !SOLO_FIABILIDAD
# ─── 1. Potencia ──────────────────────────────────────────────────────────────
ns = [100, 200, 300, 420, 600, 900]; As = [0.0, 0.1, 0.2, 0.3, 0.4]
pot = DataFrame()
for de in (1.5, 4.0)
    g = power_grid(ns, As; errores = (0.0, 1.1), nsim = NSIM, θ_de_h = de)
    g.θ_de_h .= de
    append!(pot, g)
    println("potencia calculada con dispersión de θ entre niños = $de h")
end
CSV.write(joinpath(OUT, "potencia.csv"), pot)

fig = Figure(size = figure_size(:texto, 250))
for (j, de) in enumerate((1.5, 4.0))
    ax = Axis(fig[1, j]; xlabel = "niños (n)", ylabel = j == 1 ? "potencia" : "",
              title = "Dispersión de θ entre niños: $(coma(de)) h", xticks = ns, yticks = 0:0.2:1)
    hlines!(ax, [0.8]; color = (:black, 0.4), linestyle = :dot)
    hlines!(ax, [0.05]; color = (:black, 0.25), linestyle = :dot)
    for (A, c) in zip(As[2:end], (TURQUESA, AZUL, MORADO, AMBAR))
        for (e, ls) in ((0.0, :solid), (1.1, :dash))
            d = pot[(pot.θ_de_h .== de) .& (pot.A .== A) .& (pot.error_θ_h .== e), :]
            lines!(ax, d.n, d.potencia; color = c, linestyle = ls, linewidth = 1,
                   label = e == 0 ? "A = $(coma(A)) ($(pct(exp(2 * A * 0.6) - 1)) %)" : nothing)
        end
    end
    d0 = pot[(pot.θ_de_h .== de) .& (pot.A .== 0.0) .& (pot.error_θ_h .== 0.0), :]
    scatter!(ax, d0.n, d0.potencia; color = :gray, markersize = 4, label = j == 1 ? "A = 0 (error de tipo I)" : nothing)
    ylims!(ax, 0, 1)
    j == 1 && Legend(fig[2, 1:2], ax; orientation = :horizontal, nbanks = 2)
end
Label(fig[3, 1:2], "Continua: θ y R sin error de medida. Discontinua: con el error de una semana de registro (1,1 h en θ y 0,1 en R).\n" *
      "Entre paréntesis, diferencia de HOMA-IR entre la fase mejor y la peor con R = 0,6.";
      fontsize = 7, color = (:black, 0.6), tellwidth = false)
guardar(joinpath(REP, "figures", "potencia"), fig; png = false,
        procedencia = (; tipo = :sintetico, funcion = "Makie"))

# Tabla: n necesario para 80 % de potencia (primer n de la rejilla que lo alcanza).
function nmin(de, A, e)
    d = sort(pot[(pot.θ_de_h .== de) .& (pot.A .== A) .& (pot.error_θ_h .== e), :], :n)
    i = findfirst(>=(0.8), d.potencia)
    return i === nothing ? "\\,$(">")$(maximum(ns))" : string(d.n[i])
end
open(joinpath(REP, "tables", "potencia.tex"), "w") do io
    println(io, "% Generado por scripts/07_modelo_sintetico.jl (datos sintéticos; se versiona).")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X r r r r@{}}")
    println(io, "\\toprule")
    println(io, " & \\multicolumn{2}{c}{Dispersión de \$\\theta\$: 1,5\\,h} & \\multicolumn{2}{c}{Dispersión de \$\\theta\$: 4\\,h} \\\\")
    println(io, "Efecto & sin error & con error & sin error & con error \\\\")
    println(io, "\\midrule")
    for A in As[2:end]
        println(io, "\$A = ", replace(string(A), "." => "{,}", "-" => "−"), "\$ (", pct(exp(2A * 0.6) - 1), "\\,\\%) & ",
                nmin(1.5, A, 0.0), " & ", nmin(1.5, A, 1.1), " & ", nmin(4.0, A, 0.0), " & ", nmin(4.0, A, 1.1), " \\\\")
    end
    t1 = pot[pot.A .== 0.0, :potencia]
    println(io, "\\midrule")
    println(io, "\\multicolumn{5}{@{}l}{Error de tipo I con \$A = 0\$: entre ", pct(minimum(t1)), " y ", pct(maximum(t1)),
            "\\,\\% (nominal 5\\,\\%).} \\\\")
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
end

# ─── 2. Ejemplo de ajuste completo ────────────────────────────────────────────
ej = DataFrame()
for (A, lab) in ((0.3, "con efecto"), (0.0, "sin efecto"))
    df = simulate_cohort(420; A = A, θ0_h = 18.0, θ_de_h = 4.0, error_θ_h = 1.1, error_R = 0.1, seed = 2026)
    m = fit_timing_model(df); b = bootstrap_theta0(df; B = 1000)
    push!(ej, (escenario = lab, A_verdadero = A, A = m.A, A_lo = b.A_ic[1], A_hi = b.A_ic[2],
               θ0_h = m.θ0_h, θ0_lo = b.θ0_ic_h[1], θ0_hi = b.θ0_ic_h[2], amplitud_h = b.amplitud_h,
               F = m.F, p = m.p, favorable_h = m.fase_favorable_h))
    @printf("ejemplo %s: A = %.2f [%.2f, %.2f]; θ0 = %.1f h, IC %.1f–%.1f h (%.1f h); F = %.2f, p = %.3g\n",
            lab, m.A, b.A_ic..., m.θ0_h, b.θ0_ic_h..., b.amplitud_h, m.F, m.p)
end
CSV.write(joinpath(OUT, "ejemplo_modelo.csv"), ej)   # para el plano circular (scripts/10_figuras.jl)
open(joinpath(REP, "tables", "ejemplo.tex"), "w") do io
    println(io, "% Generado por scripts/07_modelo_sintetico.jl (datos sintéticos; se versiona).")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X r r r r r@{}}")
    println(io, "\\toprule")
    println(io, "Cohorte sintética (n = 420) & \$A\$ [IC 95\\,\\%] & \$\\theta_0\$ (h) [IC 95\\,\\%] & Amplitud del IC (h) & \$F_{2,n-7}\$ & \$p\$ \\\\")
    println(io, "\\midrule")
    for r in eachrow(ej)
        println(io, "Verdadero \$A = ", replace(string(r.A_verdadero), "." => "{,}", "-" => "−"), "\$, \$\\theta_0 = 18\$\\,h & ",
                fmt(r.A), " [", fmt(r.A_lo), "; ", fmt(r.A_hi), "] & ", fmt(r.θ0_h; d = 1), " [", fmt(r.θ0_lo; d = 1), "--",
                fmt(r.θ0_hi; d = 1), "] & ", fmt(r.amplitud_h; d = 1), " & ", fmt(r.F), " & ",
                r.p < 0.001 ? "\$<0{,}001\$" : fmt(r.p; d = 3), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
end

# Figura del efecto en el plano circular (escenario con efecto).
r1 = ej[1, :]
fig = Figure(size = figure_size(:texto, 210))
ax = Axis(fig[1, 1]; xlabel = "θ: horas tras el punto medio de L5", ylabel = "HOMA-IR relativo predicho",
          xticks = 0:3:24)
h = range(0, 24; length = 241)
for (R, c) in ((0.3, TURQUESA), (0.6, AZUL), (0.9, MORADO))
    lines!(ax, h, exp.(r1.A .* R .* cos.(2π .* (h .- r1.θ0_h) ./ 24)); color = c, linewidth = 1.2, label = "R = $(coma(R))")
end
vlines!(ax, [r1.θ0_h]; color = AMBAR, linestyle = :dash)
vlines!(ax, [r1.favorable_h]; color = TURQUESA, linestyle = :dash)
lo, hi = r1.θ0_lo, r1.θ0_hi
for (a, b) in (hi >= lo ? ((lo, hi),) : ((lo, 24.0), (0.0, hi)))
    vspan!(ax, [a], [b]; color = (AMBAR, 0.12))
end
hlines!(ax, [1.0]; color = (:black, 0.3))
text!(ax, r1.θ0_h + 0.3, exp(r1.A * 0.9) + 0.02; text = "θ₀ (máximo)", align = (:left, :bottom), fontsize = 7.5, color = AMBAR)
text!(ax, r1.favorable_h + 0.3, exp(-r1.A * 0.9) - 0.02; text = "θ₀ + 12 h (más favorable)", align = (:left, :top), fontsize = 7.5, color = TURQUESA)
xlims!(ax, 0, 24)
ylims!(ax, exp(-r1.A * 0.9) - 0.08, exp(r1.A * 0.9) + 0.08)
axislegend(ax; position = :rb)
Label(fig[2, 1], "HOMA-IR predicho respecto a un niño con R = 0, para tres concentraciones. Banda ámbar: IC 95 % de θ₀ por bootstrap de niños.";
      fontsize = 7, color = (:black, 0.6), tellwidth = false)
guardar(joinpath(REP, "figures", "efecto-circular"), fig; png = false,
        procedencia = (; tipo = :sintetico, funcion = "Makie"))

# ─── 3. Censura del HOMA-IR ───────────────────────────────────────────────────
# Suelo de 0,39 ≈ insulina 2 µU/mL × glucosa 4,4 mmol/L / 22,5 (percentil 25 de MELI-POP).
bias = DataFrame(metodo = String[], A_medio = Float64[], θ0_medio_h = Float64[])
cens_frac = Float64[]; Ao = Float64[]; At = Float64[]; θo = Float64[]; θt = Float64[]
for s in 1:min(NSIM, 200)
    df = simulate_cohort(420; A = 0.3, θ0_h = 18.0, θ_de_h = 4.0, suelo = 0.39, seed = 77_000 + s)
    push!(cens_frac, count(df.censurado) / nrow(df))
    m = fit_timing_model(df); push!(Ao, m.A); push!(θo, m.θ0_h)
    φ = 2π .* df.θ_h ./ 24
    X = hcat(ones(nrow(df)), df.R .* cos.(φ), df.R .* sin.(φ), df.edad, df.sexo, df.imc_z, df.afmv)
    t = tobit_em(X, log.(df.homa), df.censurado)
    push!(At, sqrt(t.β[2]^2 + t.β[3]^2)); push!(θt, mod(atan(t.β[3], t.β[2]) * 24 / (2π), 24))
end
cm(v) = sum(v) / length(v)
# Error de tipo I del F de mínimos cuadrados con censura (cohortes sin efecto).
tipo1 = count(s -> fit_timing_model(simulate_cohort(420; A = 0.0, θ_de_h = 4.0, suelo = 0.39,
                                                    seed = 500_000 + s)).p < 0.05, 1:1000) / 1000
open(joinpath(REP, "tables", "censura.tex"), "w") do io
    println(io, "% Generado por scripts/07_modelo_sintetico.jl (datos sintéticos; se versiona).")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X r r@{}}")
    println(io, "\\toprule")
    println(io, "Método (n = 420, ", pct(cm(cens_frac)), "\\,\\% de valores en el suelo; verdadero \$A = 0{,}3\$, \$\\theta_0 = 18\$\\,h) & \$A\$ medio & \$\\theta_0\$ medio (h) \\\\")
    println(io, "\\midrule")
    println(io, "Mínimos cuadrados sobre el logaritmo & ", fmt(cm(Ao)), " & ", fmt(cm(θo); d = 1), " \\\\")
    println(io, "Tobit (EM) & ", fmt(cm(At)), " & ", fmt(cm(θt); d = 1), " \\\\")
    println(io, "\\midrule")
    println(io, "\\multicolumn{3}{@{}l}{Error de tipo I del contraste \$F\$ con censura y \$A = 0\$: ", fmt(100tipo1; d = 1),
            "\\,\\% (1000 cohortes; nominal 5\\,\\%).} \\\\")
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
end
@printf("censura: %.0f %% en el suelo; A medio MCO %.3f vs tobit %.3f (verdadero 0,3); tipo I %.3f\n",
        100cm(cens_frac), cm(Ao), cm(At), tipo1)
end # !SOLO_FIABILIDAD

# ─── 4. Fiabilidad de la exposición y atenuación del efecto ───────────────────
# Con una semana, C y S se miden con error: el efecto estimado se atenúa y la potencia baja.
# La fiabilidad λ es la del vector (la que estima Circadian.timing_reliability como
# `icc_k_vector` con la cohorte real). Para cada λ se busca el error clásico por componente
# que la produce: e² = (var C + var S)·(1/λ − 1)/2, con las varianzas verdaderas de una
# cohorte grande. Se marca también la λ equivalente al error de una semana supuesto en la
# potencia (1,1 h en θ y 0,1 en R). La atenuación se mide con la proyección del efecto
# estimado sobre la dirección verdadera (β_cos·cos θ₀ + β_sin·sin θ₀), no con A, que es un
# módulo y está sesgado hacia arriba cuando la potencia es baja.
varCS(df, θ, R) = (φ = 2π .* df[!, θ] ./ 24; c = df[!, R] .* cos.(φ); s_ = df[!, R] .* sin.(φ);
                   sum(abs2, c .- cm2(c)) / length(c) + sum(abs2, s_ .- cm2(s_)) / length(s_))
cm2(v) = sum(v) / length(v)
λs = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
fia = DataFrame(θ_de_h = Float64[], n = Int[], λ = Float64[], error_cs = Float64[],
                A_proyectado = Float64[], atenuacion = Float64[], potencia = Float64[])
λ_semana = Dict{Float64, Float64}()
for de in (1.5, 4.0)
    grande = simulate_cohort(100_000; A = 0.3, θ_de_h = de, seed = 31)
    vt = varCS(grande, :θ_verdadero_h, :R_verdadero)
    sem = simulate_cohort(100_000; A = 0.3, θ_de_h = de, error_θ_h = 1.1, error_R = 0.1, seed = 31)
    λ_semana[de] = vt / varCS(sem, :θ_h, :R)
    for n in (420, 900), λ in λs
        e = sqrt(vt * (1 / λ - 1) / 2)
        As_ = Float64[]; hits = 0
        for s in 1:NSIM
            df = simulate_cohort(n; A = 0.3, θ0_h = 18.0, θ_de_h = de, error_cs = e,
                                 seed = 880_000 + 1009s + n + round(Int, 100λ) + round(Int, 10de))
            m = fit_timing_model(df); hits += m.p < 0.05
            push!(As_, m.β_cos * cos(2π * 18 / 24) + m.β_sin * sin(2π * 18 / 24))
        end
        push!(fia, (de, n, λ, e, cm2(As_), cm2(As_) / 0.3, hits / NSIM))
    end
    @printf("fiabilidad, dispersión %.1f h: λ equivalente al error de una semana ≈ %.2f\n", de, λ_semana[de])
end
CSV.write(joinpath(OUT, "fiabilidad.csv"), fia)
CSV.write(joinpath(OUT, "fiabilidad_semana.csv"),
          DataFrame(θ_de_h = [1.5, 4.0], λ_equivalente = [λ_semana[1.5], λ_semana[4.0]]))

fig = Figure(size = figure_size(:texto, 230))
axa = Axis(fig[1, 1]; xlabel = "fiabilidad de la exposición (λ)", ylabel = "efecto estimado / verdadero",
           xticks = 0.3:0.1:1, yticks = 0:0.2:1.2, xreversed = true)
axb = Axis(fig[1, 2]; xlabel = "fiabilidad de la exposición (λ)", ylabel = "potencia",
           xticks = 0.3:0.1:1, yticks = 0:0.2:1, xreversed = true)
lines!(axa, [0.3, 1.0], [0.3, 1.0]; color = (:black, 0.35), linestyle = :dot, label = "diagonal (= λ)")
for (de, c) in ((1.5, AZUL), (4.0, AMBAR))
    d = fia[(fia.θ_de_h .== de) .& (fia.n .== 420), :]
    scatterlines!(axa, d.λ, d.atenuacion; color = c, markersize = 4, linewidth = 1,
                  label = "dispersión de θ: $(coma(de)) h")
    for (n, ls) in ((420, :solid), (900, :dash))
        d = fia[(fia.θ_de_h .== de) .& (fia.n .== n), :]
        lines!(axb, d.λ, d.potencia; color = c, linestyle = ls, linewidth = 1)
    end
    vlines!(axa, [λ_semana[de]]; color = (c, 0.5), linestyle = :dash, linewidth = 0.8)
    vlines!(axb, [λ_semana[de]]; color = (c, 0.5), linestyle = :dash, linewidth = 0.8)
end
hlines!(axb, [0.8]; color = (:black, 0.4), linestyle = :dot)
xlims!(axa, 1.03, 0.27); xlims!(axb, 1.03, 0.27); ylims!(axa, 0, 1.1); ylims!(axb, 0, 1)
axislegend(axa; position = :rt)
Label(fig[2, 1:2], "A = 0,3 (43 % de HOMA-IR entre la fase peor y la mejor con R = 0,6). Izquierda: n = 420, efecto proyectado sobre la fase verdadera.\n" *
      "Derecha: continua n = 420, discontinua n = 900. Verticales discontinuas: λ equivalente al error de una semana\n" *
      "supuesto en la potencia (1,1 h en θ y 0,1 en R), con su color de dispersión.";
      fontsize = 7, color = (:black, 0.6), tellwidth = false)
guardar(joinpath(REP, "figures", "fiabilidad"), fig; png = false,
        procedencia = (; tipo = :sintetico, funcion = "Makie"))
for r in eachrow(fia[fia.n .== 420, :])
    @printf("  dispersión %.1f h, λ = %.1f: A/A_verd = %.2f, potencia(n=420) = %.2f\n", r.θ_de_h, r.λ, r.atenuacion, r.potencia)
end
println("Tablas y figuras en ", relpath(REP, ROOT))
