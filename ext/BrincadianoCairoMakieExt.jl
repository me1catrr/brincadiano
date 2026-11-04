"""
Extensión que implementa las funciones de `Brincadiano.Viz` con CairoMakie.
Se carga automáticamente al hacer `using Brincadiano, CairoMakie`.

Convenciones de todas las figuras:
- Los huecos (`gaps`) se **sombrean**, nunca se rellenan: naranja = *idle sleep*,
  rojo = sin datos o segundo incompleto, gris = relleno de ceros (formato antiguo).
- Los ejes de varios días se rotulan «día k (día de la semana)», sin fecha de calendario.
- Paleta del informe (styles/tikz-styles.tex): azul, turquesa, ámbar, morado.
- Tamaño final en puntos (`Viz.figure_size`): `ancho = :texto` (160 mm) o `:completo`
  (258 mm, para `figure*`); tema común (`Viz.tema_informe`); sin títulos internos: el pie de
  LaTeX titula la figura. Se guardan con `Viz.guardar` (1 unidad = 1 pt).
"""
module BrincadianoCairoMakieExt

using Brincadiano
using Brincadiano.Viz
using Brincadiano.GT3X: GT3XRecording
using Brincadiano.Metrics: epoch_metrics
using Brincadiano.Viz: figure_size
using CairoMakie
using DataFrames, Dates

# Los mismos RGB que report/styles/paleta.tex (BRBlue, BRTeal, BROrange, BRPurple).
const AZUL = RGBf(13 / 255, 45 / 255, 75 / 255)
const TURQUESA = RGBf(20 / 255, 150 / 255, 145 / 255)
const AMBAR = RGBf(224 / 255, 121 / 255, 34 / 255)
const MORADO = RGBf(116 / 255, 90 / 255, 170 / 255)
const EJES = (RGBf(0.80, 0.15, 0.20), RGBf(0.15, 0.55, 0.20), RGBf(0.15, 0.35, 0.80))

const HUECOS = [("idle_sleep", "idle sleep", (AMBAR, 0.18)),
                ("missing", "sin datos", (RGBf(0.85, 0.1, 0.1), 0.22)),
                ("partial", "segundo incompleto", (RGBf(0.85, 0.1, 0.1), 0.12)),
                ("zero_padding", "relleno de ceros", (RGBf(0.5, 0.5, 0.5), 0.2))]

const GRIS = RGBf(0.25, 0.25, 0.25)

# Marcas de los ejes con coma decimal y el mismo número de decimales en todas (convención
# del informe): 0 / 0,25 / 0,50, nunca 0.0 / 0.25 / 0.5. Makie no admite un formateador
# función junto a marcas dadas como (posiciones, etiquetas): esos ejes vuelven a
# `Makie.automatic`, que respeta las etiquetas explícitas.
function _coma(vals)
    k = findfirst(d -> all(v -> abs(v * 10.0^d - round(v * 10.0^d)) < 1e-6 * max(1, abs(v * 10.0^d)), vals), 0:4)
    d = k === nothing ? 4 : k - 1
    return map(vals) do v
        d == 0 && return replace(string(round(Int, v)), "-" => "−")
        s = string(round(abs(v); digits = d))
        i = findfirst('.', s)
        s = i === nothing ? s * "," * "0"^d : replace(rpad(s, i + d, '0'), "." => ",")
        return v < 0 ? "−" * s : s
    end
end

function Viz.tema_informe()
    Theme(fontsize = 8.5, figure_padding = (2, 6, 2, 4),
          Axis = (xlabelsize = 8.5, ylabelsize = 8.5, xticklabelsize = 7.5, yticklabelsize = 7.5,
                  titlesize = 8.5, titlefont = :regular, spinewidth = 0.6, xtickwidth = 0.6, ytickwidth = 0.6,
                  xticksize = 3, yticksize = 3, xgridcolor = (:black, 0.07), ygridcolor = (:black, 0.07),
                  xgridwidth = 0.5, ygridwidth = 0.5, xlabelpadding = 2, ylabelpadding = 3,
                  topspinevisible = false, rightspinevisible = false,
                  xtickformat = _coma, ytickformat = _coma),
          Legend = (labelsize = 7.5, framevisible = false, patchsize = (10, 7), padding = (0, 0, 0, 0),
                    rowgap = 0, colgap = 8, patchlabelgap = 3),
          Label = (fontsize = 8.5,), Lines = (linewidth = 0.8,))
end

function Viz.guardar(base::AbstractString, fig; png::Bool = true, procedencia = nothing)
    save(base * ".pdf", fig; pt_per_unit = 1)
    png && save(base * ".png", fig; px_per_unit = 3)
    procedencia === nothing || _escribir_procedencia(base, procedencia)
    return nothing
end

# Ficha de procedencia junto a la figura. El script es el que se está ejecutando
# (PROGRAM_FILE); el commit, el del código que dibuja (src/, ext/, scripts/), con «+» si
# tiene cambios sin commit. Todo relativo a la raíz del repositorio.
function _escribir_procedencia(base::AbstractString, p)
    raiz = _raiz_git(dirname(abspath(base)))
    rel(x) = raiz === nothing ? basename(x) : relpath(abspath(x), raiz)
    commit = "sin git"
    if raiz !== nothing
        h = strip(read(`git -C $raiz rev-parse --short HEAD`, String))
        sucio = !isempty(strip(read(`git -C $raiz status --porcelain --untracked-files=no -- src ext scripts Project.toml Manifest.toml`, String)))
        commit = h * (sucio ? "+" : "")
    end
    tex = Viz.procedencia_tex(; tipo = p.tipo, registros = collect(String, get(p, :registros, String[])),
                              funcion = (startswith(String(p.funcion), "plot_") ? "Viz." : "") * String(p.funcion),
                              script = rel(PROGRAM_FILE),
                              entradas = collect(String, get(p, :entradas, String[])),
                              salida = rel(base), commit, fecha = today())
    write(base * ".procedencia.tex", tex)
    return nothing
end

function _raiz_git(dir)
    try
        return strip(read(pipeline(`git -C $dir rev-parse --show-toplevel`; stderr = devnull), String))
    catch
        return nothing
    end
end

# Nota explicativa al pie de una figura (gris, pequeña); el pie de LaTeX da el detalle.
_nota!(pos, txt) = Label(pos, txt; fontsize = 7, color = (:black, 0.6), tellwidth = false)

# Sombrea los huecos de `gaps` que caen en [lo, hi] (en las unidades de `tox`).
function _shade_gaps!(ax, gaps::AbstractDataFrame, tox; lo = -Inf, hi = Inf)
    for (kind, _, color) in HUECOS
        g = gaps[gaps.kind .== kind, :]
        isempty(g) && continue
        a = clamp.(tox.(g.start), lo, hi)
        b = clamp.(tox.(g.stop), lo, hi)
        keep = b .> a
        any(keep) && vspan!(ax, a[keep], b[keep]; color = color)
    end
end

function _gap_legend!(pos, gaps::AbstractDataFrame)
    present = [(lab, col) for (k, lab, col) in HUECOS if k in gaps.kind]
    isempty(present) && return nothing
    Legend(pos, [PolyElement(color = c) for (_, c) in present], first.(present);
           orientation = :horizontal, tellheight = true)
end

# Percentil sin dependencias (para fijar límites de ejes robustos a picos aislados).
function _pct(v, p)
    w = sort!(filter(isfinite, collect(Float64, v)))
    isempty(w) && return NaN
    return w[clamp(ceil(Int, p * length(w)), 1, length(w))]
end

# Duración en castellano para títulos y rótulos («1 min», «20 s», «1 h»).
function _dur(p::Period)
    ms = Dates.value(convert(Millisecond, p))
    ms % 3_600_000 == 0 && return "$(ms ÷ 3_600_000) h"
    ms % 60_000 == 0 && return "$(ms ÷ 60_000) min"
    ms % 1000 == 0 && return "$(ms ÷ 1000) s"
    return "$(ms) ms"
end

const SUP = Dict(zip("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
_pow10(k) = "10" * map(c -> SUP[c], string(k))

_origin(r::GT3XRecording) = floor(first(r.acc.time), Day)
_ndays(r::GT3XRecording) = ceil(Int, days_since(last(r.acc.time), _origin(r)))

function _day_ticks!(ax, origin, ndays)
    ax.xtickformat = Makie.automatic           # marcas con etiqueta propia (ver _coma)
    ax.xticks = (collect(0.5:1:ndays - 0.5), day_labels(origin, ndays))
    ndays > 1 && vlines!(ax, 1:ndays-1; color = (:black, 0.25), linewidth = 0.6)
    xlims!(ax, 0, ndays)
end

function _plot_overview(r::GT3XRecording; bin::Period = Minute(1), ancho::Symbol = :completo)
    origin, nd = _origin(r), _ndays(r)
    tox(t) = days_since(t, origin)
    nrm = vector_norm(r.acc)
    s = binned_summary(r.acc.time, nrm, bin; origin = origin, stop = last(r.acc.time))
    half = Millisecond(div(Dates.value(convert(Millisecond, bin)), 2))
    ok = s.n .> 0
    x = tox.(s.start[ok] .+ half)

    fig = Figure(size = figure_size(ancho, 430))
    _gap_legend!(fig[1, 1], r.gaps)
    ax1 = Axis(fig[2, 1]; ylabel = "‖a‖ (g)\nmín–máx y media")
    _shade_gaps!(ax1, r.gaps, tox)
    linesegments!(ax1, vec(permutedims(hcat(x, x))), vec(permutedims(hcat(s.vmin[ok], s.vmax[ok])));
                  color = (AZUL, 0.35), linewidth = 0.5)
    scatter!(ax1, x, s.vmean[ok]; color = AZUL, markersize = 0.9)
    ylims!(ax1, 0, max(2.0, _pct(s.vmax, 0.999)))

    h = binned_summary(r.acc.time, ones(Float32, nrow(r.acc)), Hour(1); origin = origin,
                       stop = last(r.acc.time))
    ax2 = Axis(fig[3, 1]; ylabel = "muestreo (%)")
    barplot!(ax2, tox.(h.start) .+ 1 / 48, 100 .* h.n ./ (r.info.sample_rate * 3600);
             width = 1 / 24, gap = 0, color = TURQUESA)
    ylims!(ax2, 0, 105)

    ax3 = Axis(fig[4, 1]; ylabel = "log₁₀(1 + luz bruta)")
    if nrow(r.lux) > 0
        stairs!(ax3, tox.(r.lux.time), log10.(1 .+ Float64.(r.lux.lux_raw)); step = :post,
                color = AMBAR, linewidth = 0.7)
    end

    ax4 = Axis(fig[5, 1]; ylabel = "contacto", yticks = ([0, 1], ["no", "sí"]), ytickformat = Makie.automatic)
    if nrow(r.capsense) > 0
        stairs!(ax4, tox.(r.capsense.time), Float64.(r.capsense.worn); step = :post, color = MORADO)
    end
    ylims!(ax4, -0.2, 1.2)

    ax5 = Axis(fig[6, 1]; ylabel = "batería (V)")
    nrow(r.battery) > 0 && lines!(ax5, tox.(r.battery.time), r.battery.mv ./ 1000; color = AZUL)

    axs = (ax1, ax2, ax3, ax4, ax5)
    linkxaxes!(axs...)
    for ax in axs
        _day_ticks!(ax, origin, nd)
        ax !== ax5 && hidexdecorations!(ax; grid = false, minorgrid = false)
    end
    rowsize!(fig.layout, 2, Relative(0.34))
    rowgap!(fig.layout, 4)
    return fig
end

function _plot_raw_window(r::GT3XRecording, t0::DateTime, dur::Period; title::AbstractString = "",
                          ancho::Symbol = :texto)
    t1 = t0 + dur
    sel = (r.acc.time .>= t0) .& (r.acc.time .< t1)
    secs(t) = Dates.value(t - t0) / 1000
    ts = secs.(r.acc.time[sel])
    a = r.acc[sel, :]
    tt, x, y, z, n = break_at_gaps(ts, a.x, a.y, a.z, vector_norm(a); maxstep = 1.5 / r.info.sample_rate)
    T = Dates.value(convert(Millisecond, dur)) / 1000

    fig = Figure(size = figure_size(ancho, 240))
    gw = r.gaps[(r.gaps.stop .> t0) .& (r.gaps.start .< t1), :]
    # Leyenda superior: ejes y, si los hay, huecos.
    hz = [(lab, col) for (k, lab, col) in HUECOS if k in gw.kind]
    Legend(fig[1, 1], vcat([LineElement(color = c) for c in EJES], [PolyElement(color = c) for (_, c) in hz]),
           vcat(["x", "y", "z"], first.(hz)); orientation = :horizontal, tellheight = true)
    o = 1
    ax1 = Axis(fig[1 + o, 1]; ylabel = "aceleración (g)", title = title)
    ax2 = Axis(fig[2 + o, 1]; ylabel = "‖a‖ (g)", xlabel = "tiempo (s)")
    for ax in (ax1, ax2)
        _shade_gaps!(ax, r.gaps, secs; lo = 0, hi = T)
    end
    for (v, c, lab) in zip((x, y, z), EJES, ("x", "y", "z"))
        lines!(ax1, tt, v; color = c, linewidth = 0.7, label = lab)
    end
    hlines!(ax2, [1.0]; color = (:black, 0.5), linestyle = :dash, linewidth = 0.8)
    lines!(ax2, tt, n; color = AZUL, linewidth = 0.8)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, 4)
    xlims!(ax2, 0, T)
    hidexdecorations!(ax1; grid = false)
    return fig
end

const TRAMOS = [("no uso probable", RGBf(0.25, 0.25, 0.25)),
                ("inmóvil con contacto", MORADO),
                ("sin CAPSENSE", RGBf(0.6, 0.6, 0.6))]

function _plot_days24h(r::GT3XRecording; bin::Period = Minute(1), double::Bool = false,
                       runs::Union{Nothing,AbstractDataFrame} = nothing, ancho::Symbol = :completo)
    origin, nd = _origin(r), _ndays(r)
    dev = abs.(vector_norm(r.acc) .- 1f0)
    s = binned_summary(r.acc.time, dev, bin; origin = origin,
                       stop = origin + Day(nd) - Millisecond(1))
    w = Dates.value(convert(Millisecond, bin)) / 3_600_000           # anchura en horas
    ymax = max(0.1, _pct(s.vmean, 0.995))
    span = double ? 48 : 24
    labels = day_labels(origin, nd)

    fig = Figure(size = figure_size(ancho, 43 * nd + 66))
    if runs === nothing
        _gap_legend!(fig[1, 1:2], r.gaps)
    else                                   # leyenda conjunta: huecos y tramos inmóviles
        hz = [(lab, col) for (k, lab, col) in HUECOS if k in r.gaps.kind]
        tz = [(lab, col) for (lab, col) in TRAMOS if lab in runs.etiqueta]
        Legend(fig[1, 1:2], [PolyElement(color = c) for (_, c) in vcat(hz, tz)],
               vcat(first.(hz), "tramo inmóvil: " .* first.(tz)); orientation = :horizontal,
               tellheight = true)
    end
    axs = Axis[]
    for d in 0:nd-1
        d0 = origin + Day(d)
        hrs(t) = Dates.value(t - d0) / 3_600_000
        ax = Axis(fig[d + 2, 2]; yticks = [0, floor(ymax; digits = 2)], rightspinevisible = true, topspinevisible = true)
        _shade_gaps!(ax, r.gaps, hrs; lo = 0, hi = span)
        sel = (s.start .>= d0) .& (s.start .< d0 + Hour(span)) .& (s.n .> 0)
        any(sel) && barplot!(ax, hrs.(s.start[sel]) .+ w / 2, s.vmean[sel];
                             width = w, gap = 0, color = AZUL, strokewidth = 0)
        if runs !== nothing                # franja inferior con los tramos inmóviles (H4)
            for (lab, col) in TRAMOS
                rr = runs[runs.etiqueta .== lab, :]
                isempty(rr) && continue
                a = clamp.(hrs.(rr.start), 0, span); b = clamp.(hrs.(rr.stop), 0, span)
                for (x0, x1) in zip(a, b)
                    x1 > x0 && poly!(ax, Rect2f(x0, -0.16ymax, x1 - x0, 0.12ymax); color = col)
                end
            end
        end
        # Luz bruta en eje secundario (escala log), solo como contexto.
        axl = Axis(fig[d + 2, 2]; yaxisposition = :right, ylabelcolor = AMBAR, yticklabelcolor = AMBAR)
        hidexdecorations!(axl); hidespines!(axl); hideydecorations!(axl)
        lsel = (r.lux.time .>= d0) .& (r.lux.time .< d0 + Hour(span))
        any(lsel) && stairs!(axl, hrs.(r.lux.time[lsel]), log10.(1 .+ Float64.(r.lux.lux_raw[lsel]));
                             step = :post, color = (AMBAR, 0.6), linewidth = 0.5)
        ylims!(axl, 0, 4)
        linkxaxes!(ax, axl)
        ylims!(ax, runs === nothing ? 0 : -0.18ymax, ymax)
        xlims!(ax, 0, span)
        ax.xticks = 0:3:span
        double && vlines!(ax, [24]; color = (:black, 0.4), linewidth = 0.8)
        Label(fig[d + 2, 1], labels[d + 1]; halign = :right, tellheight = false)
        d < nd - 1 && hidexdecorations!(ax; grid = false)
        push!(axs, ax)
    end
    axs[end].xlabel = double ? "hora (doble gráfica: día k y k+1)" : "hora del día"
    _nota!(fig[nd + 2, 2], "Barras: media por $(_dur(bin)) de |‖a‖ − 1 g| (g). " *
           "Línea ámbar: log₁₀(1 + luz bruta), escala 0–4.")
    colsize!(fig.layout, 1, Fixed(58))
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_histograms(r::GT3XRecording; ancho::Symbol = :texto)
    nrm = vector_norm(r.acc)
    fig = Figure(size = figure_size(ancho, 180))
    edges = 0:0.02:max(3.0, _pct(nrm, 0.9999))
    cnt = zeros(Int, length(edges) - 1)
    for v in nrm
        k = searchsortedlast(edges, v)
        1 <= k <= length(cnt) && (cnt[k] += 1)
    end
    ax1 = Axis(fig[1, 1]; xlabel = "‖a‖ (g)", ylabel = "muestras registradas", yscale = log10)
    nz = cnt .> 0
    kmax = ceil(Int, log10(maximum(cnt)))
    ax1.ytickformat = Makie.automatic
    ax1.yticks = (10.0 .^ (0:kmax), _pow10.(0:kmax))
    barplot!(ax1, (edges[1:end-1] .+ 0.01)[nz], cnt[nz]; width = 0.02, gap = 0, color = AZUL,
             fillto = 0.5, strokewidth = 0)
    vlines!(ax1, [1.0]; color = (:black, 0.5), linestyle = :dash)

    ax2 = Axis(fig[1, 2]; xlabel = "log₁₀(1 + luz bruta)", ylabel = "registros de luz (al cambiar)")
    if nrow(r.lux) > 0
        hist!(ax2, log10.(1 .+ Float64.(r.lux.lux_raw)); bins = 40, color = AMBAR)
    end
    return fig
end

function _plot_metrics_window(r::GT3XRecording, t0::DateTime, dur::Period;
                                 epochs = (Second(1), Second(5), Second(60)), title::AbstractString = "",
                                 offset = (0.0, 0.0, 0.0), scale = (1.0, 1.0, 1.0), ancho::Symbol = :texto)
    t1 = t0 + dur
    secs(t) = Dates.value(t - t0) / 1000
    T = Dates.value(convert(Millisecond, dur)) / 1000
    sel = (r.acc.time .>= t0) .& (r.acc.time .< t1)
    a = r.acc[sel, :]
    cx = (a.x .+ offset[1]) .* scale[1]; cy = (a.y .+ offset[2]) .* scale[2]; cz = (a.z .+ offset[3]) .* scale[3]
    tt, n = break_at_gaps(secs.(a.time), vector_norm(cx, cy, cz); maxstep = 1.5 / r.info.sample_rate)
    fig = Figure(size = figure_size(ancho, 70 + 62 * (length(epochs) + 1)))
    gw = r.gaps[(r.gaps.stop .> t0) .& (r.gaps.start .< t1), :]
    o = isempty(gw) ? 0 : 1
    o == 1 && _gap_legend!(fig[1, 1], gw)
    ax0 = Axis(fig[1 + o, 1]; ylabel = "‖a‖ (g)", title = title)
    _shade_gaps!(ax0, r.gaps, secs; lo = 0, hi = T)
    band!(ax0, tt, fill(1.0, length(n)), max.(n, 1.0); color = (AMBAR, 0.45))   # área = ENMO
    lines!(ax0, tt, n; color = AZUL, linewidth = 0.6)
    hlines!(ax0, [1.0]; color = (:black, 0.5), linestyle = :dash, linewidth = 0.6)
    axs = Axis[ax0]
    # Todos los epochs en la misma escala, para ver cómo se diluyen los picos.
    eps = [epoch_metrics(r.acc, r.gaps, r.info.sample_rate, e; origin = t0 - Hour(1), stop = t1,
                         offset = offset, scale = scale) for e in epochs]
    eps = [e[(e.start .>= t0 - Millisecond(1)) .& (e.start .< t1), :] for e in eps]
    ymax = maximum(e -> maximum(filter(isfinite, e.enmo_idle0); init = 0.0), eps) * 1.05
    for (j, (p, e)) in enumerate(zip(epochs, eps))
        ms = Dates.value(convert(Millisecond, p)) / 1000
        ax = Axis(fig[1 + o + j, 1]; ylabel = "ENMO\n$(_dur(p)) (mg)")
        _shade_gaps!(ax, r.gaps, secs; lo = 0, hi = T)
        ok = isfinite.(e.enmo_idle0)
        any(ok) && barplot!(ax, secs.(e.start[ok]) .+ ms / 2, e.enmo_idle0[ok]; width = ms, gap = 0.05,
                            color = TURQUESA, strokewidth = 0)
        ylims!(ax, 0, max(ymax, 1))
        push!(axs, ax)
    end
    linkxaxes!(axs...)
    for ax in axs
        xlims!(ax, 0, T)
        ax !== axs[end] && hidexdecorations!(ax; grid = false)
    end
    axs[end].xlabel = "tiempo (s)"
    rowgap!(fig.layout, 3)
    return fig
end

# Tramo [a, a + L) en horas del reloj, partido en dos si cruza la medianoche.
_wrap_spans(a, L) = a + L <= 24 ? [(a, a + L)] : [(a, 24.0), (0.0, a + L - 24)]

function _plot_profile24h(perfiles, rar; smooth::Integer = 15, title::AbstractString = "", ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 200))
    ax = Axis(fig[2, 1]; xlabel = "hora del día", ylabel = "ENMO (mg)", title = title,
              xticks = 0:3:24)
    # Ventanas L5 (morado) y M10 (turquesa) del ancla elegida.
    for (a, b) in _wrap_spans(rar.L5_inicio_h, 5)
        vspan!(ax, [a], [b]; color = (MORADO, 0.15))
    end
    for (a, b) in _wrap_spans(rar.M10_inicio_h, 10)
        vspan!(ax, [a], [b]; color = (TURQUESA, 0.12))
    end
    vlines!(ax, [rar.L5_medio_h]; color = MORADO, linewidth = 1, linestyle = :dash)
    colores = (AZUL, AMBAR, RGBf(0.5, 0.5, 0.5))
    for (j, (lab, p)) in enumerate(perfiles)
        v = p.media; n = length(v); h = (0:n-1) ./ 60
        sm = [begin w = filter(isfinite, v[mod.((i - smooth ÷ 2):(i + smooth ÷ 2), n) .+ 1])
                  length(w) >= smooth ÷ 2 ? sum(w) / length(w) : NaN end for i in 0:n-1]
        lines!(ax, h, sm; color = colores[j], linewidth = 1, label = lab)
    end
    xlims!(ax, 0, 24); ylims!(ax, 0, nothing)
    Legend(fig[1, 1], vcat([LineElement(color = colores[j]) for j in eachindex(perfiles)],
                           [PolyElement(color = (MORADO, 0.3)), PolyElement(color = (TURQUESA, 0.25))]),
           vcat([lab for (lab, _) in perfiles], ["L5 (línea: punto medio)", "M10"]); orientation = :horizontal, tellheight = true)
    _nota!(fig[3, 1], "Media móvil de $(smooth) min.")
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_timing_circular(fase, dosis, dia, θ, R; title::AbstractString = "",
                               ancla::AbstractString = "punto medio de L5", ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 300))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = title)
    hidedecorations!(ax); hidespines!(ax)
    pos(r, φ) = (r * sin(φ), r * cos(φ))               # 0 arriba, sentido horario
    φs = range(0, 2π; length = 361)
    lines!(ax, [pos(1, φ)[1] for φ in φs], [pos(1, φ)[2] for φ in φs]; color = AZUL, linewidth = 0.8)
    for h in 0:3:21
        φ = 2π * h / 24
        a, b = pos(1, φ), pos(0.94, φ)
        lines!(ax, [a[1], b[1]], [a[2], b[2]]; color = AZUL)
        t = pos(1.12, φ)
        text!(ax, t[1], t[2]; text = h == 0 ? "0 h" : "+$(h) h", align = (:center, :center), fontsize = 7.5)
    end
    # Episodios: un anillo por día, tamaño ∝ √dosis.
    nd = isempty(dia) ? 1 : maximum(dia)
    r = 0.45 .+ 0.45 .* (dia .- 1) ./ max(nd - 1, 1)
    ms = 2 .+ 6 .* sqrt.(dosis ./ maximum(dosis; init = 1.0))
    scatter!(ax, [pos(ri, φ)[1] for (ri, φ) in zip(r, fase)], [pos(ri, φ)[2] for (ri, φ) in zip(r, fase)];
             color = (AMBAR, 0.55), markersize = ms, strokewidth = 0)
    e = pos(R, θ)
    arrows2d!(ax, [0.0], [0.0], [e[1]], [e[2]]; color = MORADO, shaftwidth = 2, tipwidth = 8)
    text!(ax, 0, -1.32; text = "θ = $(_coma([round(θ * 24 / (2π); digits = 1)])[1]) h tras el $(ancla) · R = $(_coma([round(R; digits = 2)])[1])",
          align = (:center, :center), fontsize = 8)
    xlims!(ax, -1.3, 1.3); ylims!(ax, -1.42, 1.3)
    return fig
end

const METODOS = Dict("postura" => (TURQUESA, "reposo nocturno estimado (postura)"),
                     "retirada" => (MORADO, "reposo supuesto (retirada, SPT*)"))

function _plot_sleep_nights(ep::AbstractDataFrame, noches::AbstractDataFrame; axis::Integer,
                            threshold::Real = 60, title::AbstractString = "", ancho::Symbol = :completo)
    col = (:anglex, :angley, :anglez)[axis]
    nn = nrow(noches)
    fig = Figure(size = figure_size(ancho, 52 * nn + 70))
    title == "" || Label(fig[0, 1:2], title; font = :bold)
    Legend(fig[1, 1:2], [LineElement(color = AZUL), PolyElement(color = (TURQUESA, 0.25)),
                         PolyElement(color = (MORADO, 0.25)), PolyElement(color = TURQUESA),
                         PolyElement(color = (RGBf(0.25, 0.25, 0.25), 0.35)), PolyElement(color = (AMBAR, 0.18))],
           ["|ángulo longitudinal|", "reposo nocturno estimado (postura)", "reposo supuesto (retirada)",
            "inactividad sostenida (SIB)", "no uso probable", "ángulo heredado (idle sleep)"];
           orientation = :horizontal, tellheight = true)
    nidx = Dict(d => k for (k, d) in enumerate(sort(unique(noches.noche))))
    axs = Axis[]
    for (r, n) in enumerate(eachrow(noches))
        d0 = DateTime(n.noche)
        hrs(t) = Dates.value(t - d0) / 3_600_000
        k = findall(t -> d0 + Hour(12) <= t < d0 + Hour(36), ep.start)
        ax = Axis(fig[r + 1, 2]; yticks = [0, threshold, 90], rightspinevisible = true, topspinevisible = true)
        # Medias por minuto (12 epochs de 5 s) para aligerar el dibujo.
        x = hrs.(ep.start[k]); y = abs.(ep[k, col])
        m = 12; nb = length(k) ÷ m
        xm = [x[(i-1)*m+1] for i in 1:nb]
        ym = [(v = filter(isfinite, y[(i-1)*m+1:i*m]); isempty(v) ? NaN : sum(v) / length(v)) for i in 1:nb]
        for (mask, c) in ((ep.nonwear[k], (RGBf(0.25, 0.25, 0.25), 0.35)), (ep.imputado[k], (AMBAR, 0.18)))
            for (i, j) in _bool_runs(mask)
                vspan!(ax, x[i], x[j] + 5 / 3600; color = c)
            end
        end
        if n.metodo in keys(METODOS)
            c = METODOS[n.metodo][1]
            vspan!(ax, n.inicio_h, n.fin_h; color = (c, 0.25))
            if n.metodo == "postura"
                ins = findall(t -> n.inicio <= t < n.fin, ep.start[k])
                sib = Brincadiano.Sleep.sib_vanhees2015(ep.anglez[k])
                for (i, j) in _bool_runs(sib[ins])
                    poly!(ax, Rect2f(x[ins[i]], -14, x[ins[j]] - x[ins[i]] + 5 / 3600, 10); color = TURQUESA)
                end
            end
        end
        hlines!(ax, [threshold]; color = (:black, 0.5), linestyle = :dash, linewidth = 0.5)
        lines!(ax, xm, ym; color = AZUL, linewidth = 0.6)
        ylims!(ax, -16, 95); xlims!(ax, 12, 36)
        ax.xtickformat = Makie.automatic
        ax.xticks = (12:3:36, [lpad(h % 24, 2, '0') * ":00" for h in 12:3:36])
        dn = Viz.DIAS_SEMANA[dayofweek(n.noche)]; dn2 = Viz.DIAS_SEMANA[dayofweek(n.noche + Day(1))]
        Label(fig[r + 1, 1], "noche $(nidx[n.noche]) ($dn→$dn2)"; halign = :right, tellheight = false)
        r < nn && hidexdecorations!(ax; grid = false)
        push!(axs, ax)
    end
    isempty(axs) || (axs[end].xlabel = "hora (de 12:00 a 12:00)")
    _nota!(fig[nn + 2, 2], "Ángulo del eje longitudinal con la horizontal (grados, media por minuto): " *
           "de pie ≈ 90°, tumbado < $(threshold)° (línea discontinua).")
    colsize!(fig.layout, 1, Fixed(72))
    rowgap!(fig.layout, 3)
    return fig
end

const FONDO_LUZ = [("sueño", (TURQUESA, 0.18)), ("retirado", (RGBf(0.25, 0.25, 0.25), 0.30)),
                   ("sin dato", (AMBAR, 0.22))]
const BARRAS_LUZ = [("válido", AZUL), ("oscuro y activo", RGBf(0.80, 0.15, 0.20)),
                    ("oscuro e inactivo", MORADO)]

function _plot_light_days(d::AbstractDataFrame; thr = (100, 1000), title::AbstractString = "", ancho::Symbol = :completo)
    origin = floor(first(d.start), Day)
    nd = ceil(Int, days_since(last(d.start), origin))
    labels = day_labels(origin, nd)
    fig = Figure(size = figure_size(ancho, 43 * nd + 66))
    title == "" || Label(fig[0, 1:2], title; font = :bold)
    Legend(fig[1, 1:2], vcat([PolyElement(color = c) for (_, c) in BARRAS_LUZ], [PolyElement(color = c) for (_, c) in FONDO_LUZ]),
           vcat(["luz: " * l for (l, _) in BARRAS_LUZ], ["fondo: " * (l == "sueño" ? "reposo nocturno estimado" : l) for (l, _) in FONDO_LUZ]);
           orientation = :horizontal, tellheight = true)
    ymax = log10(1 + 3000)
    axs = Axis[]
    for k in 0:nd-1
        d0 = origin + Day(k)
        sel = findall(t -> d0 <= t < d0 + Day(1), d.start)
        ax = Axis(fig[k + 2, 2]; yticks = ([0, log10(1 + thr[1]), log10(1 + thr[2])], ["0", string(thr[1]), string(thr[2])]),
                  ytickformat = Makie.automatic,
                  rightspinevisible = true, topspinevisible = true)
        x = [Dates.value(t - d0) / 3_600_000 for t in d.start[sel]]
        for (lab, col) in FONDO_LUZ
            for (i, j) in _bool_runs(d.estado[sel] .== lab)
                vspan!(ax, x[i], x[j] + 1 / 60; color = col)
            end
        end
        for (lab, col) in BARRAS_LUZ
            m = (d.estado[sel] .== lab)
            any(m) || continue
            y = log10.(1 .+ d.lux_media[sel][m])
            lab == "válido" || (y .= max.(y, 0.12))           # que lo oscuro se vea
            barplot!(ax, x[m] .+ 1 / 120, y; width = 1 / 60, gap = 0, color = col, strokewidth = 0)
        end
        hlines!(ax, [log10(1 + thr[1]), log10(1 + thr[2])]; color = (:black, 0.45), linestyle = :dash, linewidth = 0.5)
        ylims!(ax, 0, ymax); xlims!(ax, 0, 24)
        ax.xticks = 0:3:24
        Label(fig[k + 2, 1], labels[k + 1]; halign = :right, tellheight = false)
        k < nd - 1 && hidexdecorations!(ax; grid = false)
        push!(axs, ax)
    end
    axs[end].xlabel = "hora del día"
    _nota!(fig[nd + 2, 2], "Barras: log₁₀(1 + luz bruta) por minuto (escala bruta: la conversión a lux no está resuelta). " *
           "Líneas: $(thr[1]) y $(thr[2]) brutos.")
    colsize!(fig.layout, 1, Fixed(58))
    rowgap!(fig.layout, 3)
    return fig
end


# ─── Figuras de síntesis (H13) ────────────────────────────────────────────────────────

# = BREstNoValido de report/styles/paleta.tex (rojo teja), para que el rojo de las
# figuras sea el mismo que el de las etiquetas «no válido» del informe.
const ROJO = RGBf(170 / 255, 52 / 255, 40 / 255)
const COLORES_ID = (AZUL, AMBAR, MORADO, TURQUESA)
const MAPA_LUZ = cgrad([RGBf(0.12, 0.12, 0.12), AMBAR, RGBf(1.0, 0.92, 0.35)])
const SIN_DATO = RGBf(0.86, 0.86, 0.86)
_hhmm_ticks(hs) = (collect(hs), [lpad(h % 24, 2, '0') * ":00" for h in hs])

# Media móvil circular de `k` minutos (al menos la mitad de valores finitos).
function _suaviza(v, k)
    n = length(v)
    return [begin w = filter(isfinite, v[mod.((i - k ÷ 2):(i + k ÷ 2), n) .+ 1])
                length(w) >= k ÷ 2 ? sum(w) / length(w) : NaN end for i in 0:n-1]
end

function _plot_actogram(e60, luz, noches, episodios, capsense; ancho::Symbol = :completo)
    origin = floor(first(e60.start), Day)
    nd = ceil(Int, days_since(last(e60.start), origin))
    labels = day_labels(origin, nd)
    v = e60.enmo_idle0
    ymax = max(50.0, _pct(v[.!e60.nonwear .& isfinite.(v)], 0.995))
    fig = Figure(size = figure_size(ancho, 50 * nd + 85))
    Legend(fig[1, 1:2], [PolyElement(color = AZUL), MarkerElement(marker = :dtriangle, color = ROJO, markersize = 6),
                         PolyElement(color = (TURQUESA, 0.25)), PolyElement(color = (GRIS, 0.3)),
                         PolyElement(color = MAPA_LUZ[1]), PolyElement(color = get(MAPA_LUZ, 0.5)),
                         PolyElement(color = MAPA_LUZ[end]), PolyElement(color = SIN_DATO), PolyElement(color = MORADO)],
           ["ENMO por minuto", "episodio vigoroso", "reposo nocturno estimado", "no uso probable",
            "luz 0", "luz ≈ 30", "luz ≥ 1000 (bruta)", "luz sin dato", "contacto con la piel"];
           orientation = :horizontal, tellheight = true, nbanks = 2)
    tox(t, d0) = Dates.value(t - d0) / 3_600_000
    axs = Axis[]
    for k in 0:nd-1
        d0 = origin + Day(k); d1 = d0 + Day(1)
        ax = Axis(fig[k + 2, 2]; yticks = [0, floor(Int, ymax)], rightspinevisible = true, topspinevisible = true)
        sel = findall(t -> d0 <= t < d1, e60.start)
        x = tox.(e60.start[sel], d0)
        for n in eachrow(noches)
            n.metodo == "sin ventana" && continue
            a = max(tox(n.inicio, d0), 0.0); b = min(tox(n.fin, d0), 24.0)
            b > a && vspan!(ax, a, b; color = (TURQUESA, 0.2))
        end
        for (i, j) in _bool_runs(e60.nonwear[sel])
            vspan!(ax, x[i], x[j] + 1 / 60; color = (GRIS, 0.25))
        end
        vs = v[sel]; ok = isfinite.(vs)
        any(ok) && barplot!(ax, x[ok] .+ 1 / 120, min.(vs[ok], ymax); width = 1 / 60, gap = 0, color = AZUL, strokewidth = 0)
        ls = findall(t -> d0 <= t < d1, luz.start)
        if !isempty(ls)
            col = [luz.lux_obs[i] >= 0.5 && isfinite(luz.lux_media[i]) ?
                   get(MAPA_LUZ, clamp(log10(1 + luz.lux_media[i]) / 3, 0, 1)) : SIN_DATO for i in ls]
            barplot!(ax, tox.(luz.start[ls], d0) .+ 1 / 120, fill(-0.06ymax, length(ls)); fillto = -0.2ymax,
                     width = 1 / 60, gap = 0, color = col, strokewidth = 0)
        end
        cs = findall(t -> d0 <= t < d1, capsense.time)
        if !isempty(cs)
            m = Bool.(capsense.worn[cs])
            any(m) && barplot!(ax, tox.(capsense.time[cs][m], d0) .+ 1 / 120, fill(-0.24ymax, count(m));
                               fillto = -0.32ymax, width = 1 / 60, gap = 0, color = MORADO, strokewidth = 0)
        end
        es = findall(t -> d0 <= t < d1, episodios.medio)
        isempty(es) || scatter!(ax, tox.(episodios.medio[es], d0), fill(0.93ymax, length(es));
                                marker = :dtriangle, color = ROJO, markersize = 4)
        ylims!(ax, -0.34ymax, ymax); xlims!(ax, 0, 24)
        ax.xticks = 0:3:24
        Label(fig[k + 2, 1], labels[k + 1]; halign = :right, tellheight = false)
        k < nd - 1 && hidexdecorations!(ax; grid = false)
        push!(axs, ax)
    end
    axs[end].xlabel = "hora del día"
    _nota!(fig[nd + 2, 2], "Barras: ENMO por minuto (mg, idle sleep = 0). Franjas inferiores: luz (log, por minuto) y contacto con la piel (CAPSENSE).")
    colsize!(fig.layout, 1, Fixed(58))
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_profiles_compare(perfiles; smooth::Integer = 15, ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 200))
    ax = Axis(fig[2, 1]; xlabel = "hora del día", ylabel = "ENMO (mg)", xticks = 0:3:24)
    for (j, (lab, p, l5)) in enumerate(perfiles)
        sm = _suaviza(p.media, smooth)
        lines!(ax, (0:length(sm)-1) ./ 60, sm; color = COLORES_ID[j], linewidth = 1)
        vlines!(ax, [l5]; color = COLORES_ID[j], linestyle = :dash, linewidth = 0.8)
    end
    xlims!(ax, 0, 24); ylims!(ax, 0, nothing)
    Legend(fig[1, 1], [LineElement(color = COLORES_ID[j]) for j in eachindex(perfiles)], [p[1] for p in perfiles];
           orientation = :horizontal, tellheight = true)
    _nota!(fig[3, 1], "Media móvil de $(smooth) min. Discontinuas: punto medio de L5 de cada registro (el ancla del timing).")
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_enmo_ecdf(series; umbrales = (63.3, 142.6, 464.6), ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 200))
    xt = [0.1, 1, 10, 100, 1000]
    ax = Axis(fig[2, 1]; xscale = log10, xlabel = "ENMO (mg, escala logarítmica)", ylabel = "fracción acumulada",
              xticks = (xt, ["0,1", "1", "10", "100", "1000"]), xtickformat = Makie.automatic)
    for (j, (lab, v)) in enumerate(series)
        w = sort!(max.(filter(isfinite, collect(Float64, v)), 0.1))
        isempty(w) && continue
        q = range(0, 1; length = 2001)
        xs = [w[clamp(ceil(Int, p * length(w)), 1, length(w))] for p in q]
        stairs!(ax, xs, collect(q); step = :post, color = COLORES_ID[j], linewidth = 1)
    end
    for (u, nom) in zip(umbrales, ("ligera", "moderada", "vigorosa"))
        vlines!(ax, [u]; color = (:black, 0.4), linestyle = :dash, linewidth = 0.6)
        text!(ax, u, 0.04; text = " " * nom, align = (:left, :bottom), fontsize = 7, color = (:black, 0.6))
    end
    xlims!(ax, 0.1, 3000); ylims!(ax, 0, 1)
    Legend(fig[1, 1], [LineElement(color = COLORES_ID[j]) for j in eachindex(series)], [s[1] for s in series];
           orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_timing_compare(datos; ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 300))
    ax = Axis(fig[2, 1]; aspect = DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    pos(r, φ) = (r * sin(φ), r * cos(φ))
    φs = range(0, 2π; length = 361)
    lines!(ax, [pos(1, φ)[1] for φ in φs], [pos(1, φ)[2] for φ in φs]; color = (:black, 0.5), linewidth = 0.8)
    for h in 0:3:21
        φ = 2π * h / 24
        a, b = pos(1, φ), pos(0.95, φ)
        lines!(ax, [a[1], b[1]], [a[2], b[2]]; color = (:black, 0.5))
        t = pos(1.12, φ)
        text!(ax, t[1], t[2]; text = h == 0 ? "0 h" : "+$(h) h", align = (:center, :center), fontsize = 7.5)
    end
    nd = length(datos)
    for (j, (lab, fase, dosis, θ, R)) in enumerate(datos)
        r0 = 0.3 + 0.6 * (j - 1) / max(nd, 1); dr = 0.5 / max(nd, 1)
        r = [r0 + dr * ((i * 37) % 11) / 10 for i in eachindex(fase)]        # dispersión fija, solo visual
        ms = 1.5 .+ 4 .* sqrt.(dosis ./ maximum(dosis; init = 1.0))
        scatter!(ax, [pos(ri, φ)[1] for (ri, φ) in zip(r, fase)], [pos(ri, φ)[2] for (ri, φ) in zip(r, fase)];
                 color = (COLORES_ID[j], 0.45), markersize = ms, strokewidth = 0)
        e = pos(R, θ)
        arrows2d!(ax, [0.0], [0.0], [e[1]], [e[2]]; color = COLORES_ID[j], shaftwidth = 2, tipwidth = 7)
    end
    xlims!(ax, -1.25, 1.25); ylims!(ax, -1.25, 1.25)
    Legend(fig[1, 1], [MarkerElement(marker = :circle, color = COLORES_ID[j], markersize = 6) for j in 1:nd],
           ["$(d[1]): θ = $(_coma([round(d[4] * 24 / (2π); digits = 1)])[1]) h, R = $(_coma([round(d[5]; digits = 2)])[1])" for d in datos];
           orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 2)
    return fig
end

function _plot_sleep_raster(datos; ggir::Bool = false, ancho::Symbol = :texto)
    nrows = sum(nrow(n) for (_, n) in datos)
    fig = Figure(size = figure_size(ancho, 60 + (ggir ? 18 : 13) * nrows + 8 * length(datos)))
    ax = Axis(fig[2, 1]; xlabel = "hora", xticks = _hhmm_ticks(18:3:36), xtickformat = Makie.automatic,
              ytickformat = Makie.automatic, yreversed = true, ygridvisible = false)
    y = 0.0; ticks = Float64[]; labs = String[]
    for (lab, n) in datos
        for (k, r) in enumerate(eachrow(n))
            y += 1
            push!(ticks, y); push!(labs, "$(lab) · $(k) ($(Viz.DIAS_SEMANA[r.dia_semana]))")
            isfinite(r.inicio_h) || continue
            c = r.libre ? AMBAR : AZUL
            α = r.metodo == "retirada" ? 0.35 : 0.9
            h = ggir ? 0.34 : 0.6
            y0 = ggir ? y - 0.38 : y - h / 2
            poly!(ax, Rect2f(r.inicio_h, y0, r.fin_h - r.inicio_h, h); color = (c, α))
            if ggir && hasproperty(r, :inicio_ggir) && isfinite(r.inicio_ggir)
                poly!(ax, Rect2f(r.inicio_ggir, y + 0.04, r.fin_ggir - r.inicio_ggir, 0.34); color = (GRIS, 0.55))
            end
            if hasproperty(r, :L5_noche_h) && isfinite(r.L5_noche_h)
                scatter!(ax, [r.L5_noche_h], [y0 + h / 2]; marker = :diamond, color = MORADO, markersize = 5)
            end
        end
        y += 0.8
    end
    ax.yticks = (ticks, labs)
    xlims!(ax, 18, 37); ylims!(ax, y + 0.2, 0.2)
    els = Any[PolyElement(color = AZUL), PolyElement(color = AMBAR), PolyElement(color = (AZUL, 0.35))]
    nms = ["noche de diario", "noche libre (vie, sáb)", "SPT* (retirada)"]
    any(hasproperty(n, :L5_noche_h) for (_, n) in datos) &&
        (push!(els, MarkerElement(marker = :diamond, color = MORADO, markersize = 6)); push!(nms, "punto medio de L5 de la noche"))
    ggir && (push!(els, PolyElement(color = (GRIS, 0.55))); push!(nms, "GGIR (tiempo en cama, eje Y)"))
    Legend(fig[1, 1], els, nms; orientation = :horizontal, tellheight = true, nbanks = ggir ? 2 : 1)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_vigorous_summary(r::GT3XRecording, ventana, datos; offset = (0.0, 0.0, 0.0), scale = (1.0, 1.0, 1.0),
                                umbral_mg::Real = 464.6, ancho::Symbol = :completo)
    t0, t1, e0, e1 = ventana
    fig = Figure(size = figure_size(ancho, 300))
    secs(t) = Dates.value(t - t0) / 1000
    sel = (r.acc.time .>= t0) .& (r.acc.time .< t1)
    a = r.acc[sel, :]
    nr = vector_norm((a.x .+ offset[1]) .* scale[1], (a.y .+ offset[2]) .* scale[2], (a.z .+ offset[3]) .* scale[3])
    tt, nn = break_at_gaps(secs.(a.time), nr; maxstep = 1.5 / r.info.sample_rate)
    ax1 = Axis(fig[2, 1:2]; xlabel = "tiempo (s)", ylabel = "‖a‖ (g)")
    vspan!(ax1, secs(e0), secs(e1); color = (ROJO, 0.12))
    hlines!(ax1, [1.0]; color = (:black, 0.4), linestyle = :dash, linewidth = 0.6)
    hlines!(ax1, [1 + umbral_mg / 1000]; color = ROJO, linestyle = :dot, linewidth = 0.8)
    lines!(ax1, tt, nn; color = AZUL, linewidth = 0.6)
    xlims!(ax1, 0, secs(t1))
    # Duraciones (s) y dosis por hora del día.
    edges = [3, 4, 5, 6, 8, 10, 15, 20, 30, 1_000]
    labsd = ["3", "4", "5", "6–7", "8–9", "10–14", "15–19", "20–29", "≥ 30"]
    ax2 = Axis(fig[3, 1]; xlabel = "duración del episodio (s)", ylabel = "% de episodios",
               xticks = (1:length(labsd), labsd), xtickformat = Makie.automatic)
    ax3 = Axis(fig[3, 2]; xlabel = "hora del día", ylabel = "dosis por día (g·s)", xticks = 0:3:24)
    nd = length(datos)
    for (j, (lab, ep, dias)) in enumerate(datos)
        c = [count(d -> edges[k] <= d < edges[k+1], ep.duracion_s) for k in 1:length(labsd)]
        barplot!(ax2, (1:length(labsd)) .+ (j - (nd + 1) / 2) * 0.8 / nd, 100 .* c ./ max(nrow(ep), 1);
                 width = 0.8 / nd, color = COLORES_ID[j], strokewidth = 0)
        h = [sum(ep.dosis[hour.(ep.medio) .== k]; init = 0.0) / 1000 / max(dias, 1) for k in 0:23]
        barplot!(ax3, (0:23) .+ 0.5 .+ (j - (nd + 1) / 2) * 0.8 / nd, h; width = 0.8 / nd, color = COLORES_ID[j], strokewidth = 0)
    end
    xlims!(ax3, 0, 24); ylims!(ax2, 0, nothing); ylims!(ax3, 0, nothing)
    Legend(fig[1, 1:2], vcat([PolyElement(color = COLORES_ID[j]) for j in 1:nd],
                             [PolyElement(color = (ROJO, 0.2)), LineElement(color = ROJO, linestyle = :dot)]),
           vcat([d[1] for d in datos], ["episodio del ejemplo", "‖a‖ = 1 g + umbral vigoroso"]);
           orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 4)
    return fig
end

function _plot_circular_plane(A::Real, θ0_h::Real, puntos; ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 290))
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    pos(R, h) = (R * sin(2π * h / 24), R * cos(2π * h / 24))
    g = range(-1, 1; length = 601)
    φ0 = 2π * θ0_h / 24
    z = [hypot(x, y) <= 1 ? exp(A * hypot(x, y) * cos(atan(x, y) - φ0)) : NaN for x in g, y in g]
    hm = heatmap!(ax, g, g, z; colormap = cgrad([TURQUESA, RGBf(1, 1, 1), AMBAR]), colorrange = (exp(-A), exp(A)), nan_color = :transparent)
    φs = range(0, 2π; length = 361)
    for (R, ls) in ((1.0, :solid), (0.5, :dash))
        lines!(ax, [R * sin(φ) for φ in φs], [R * cos(φ) for φ in φs]; color = (:black, 0.5), linewidth = 0.6, linestyle = ls)
    end
    for h in 0:3:21
        t = pos(1.12, h)
        text!(ax, t[1], t[2]; text = h == 0 ? "0 h" : "+$(h) h", align = (:center, :center), fontsize = 7.5)
    end
    e = pos(1.0, θ0_h)
    lines!(ax, [0, e[1]], [0, e[2]]; color = AMBAR, linestyle = :dash, linewidth = 0.8)
    text!(ax, e[1] * 1.02, e[2] * 1.02; text = "θ₀", align = (:left, :bottom), fontsize = 8, color = AMBAR)
    for (j, (lab, R, θh)) in enumerate(puntos)
        p = pos(R, θh)
        scatter!(ax, [p[1]], [p[2]]; color = COLORES_ID[j], markersize = 7, strokecolor = :white, strokewidth = 0.8)
        text!(ax, p[1], p[2]; text = "  " * lab, align = (:left, :center), fontsize = 7.5, color = COLORES_ID[j])
    end
    xlims!(ax, -1.25, 1.25); ylims!(ax, -1.25, 1.25)
    Colorbar(fig[1, 2], hm; label = "HOMA-IR relativo predicho\n(cohorte sintética)", ticklabelsize = 7.5, labelsize = 8,
             width = 6, height = Relative(0.6))
    return fig
end

function _plot_bland_altman(datos; unidad::AbstractString = "mg", ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 200))
    for (j, (lab, a, b)) in enumerate(datos)
        m = (a .+ b) ./ 2; d = a .- b
        ax = Axis(fig[1, j]; xlabel = "media de los dos ($unidad)", ylabel = j == 1 ? "diferencia, nuestro − GGIR ($unidad)" : "", title = lab)
        xmax = _pct(m, 0.995); dq = 1.15 * _pct(abs.(d), 0.999)
        keep = (m .<= xmax) .& (abs.(d) .<= dq)
        nx, ny = 70, 50
        cnt = zeros(nx, ny)
        for (mi, di) in zip(m[keep], d[keep])
            ix = clamp(floor(Int, mi / xmax * nx) + 1, 1, nx); iy = clamp(floor(Int, (di + dq) / (2dq) * ny) + 1, 1, ny)
            cnt[ix, iy] += 1
        end
        zz = [c > 0 ? log10(c) : NaN for c in cnt]
        heatmap!(ax, range(0, xmax; length = nx), range(-dq, dq; length = ny), zz;
                 colormap = cgrad([RGBf(0.80, 0.86, 0.92), AZUL]), nan_color = :transparent)
        μ = sum(d) / length(d); sd = sqrt(sum((d .- μ) .^ 2) / (length(d) - 1))
        hlines!(ax, [μ]; color = AMBAR, linewidth = 0.9)
        hlines!(ax, [μ - 1.96sd, μ + 1.96sd]; color = AMBAR, linestyle = :dash, linewidth = 0.8)
        ylims!(ax, -dq, dq); xlims!(ax, 0, xmax)
    end
    _nota!(fig[2, :], "Color: número de epochs (escala log). Línea ámbar: diferencia media; discontinuas: ±1,96 DE.")
    return fig
end

function _plot_agreement_timeline(datos; paso::Period = Minute(15), ancho::Symbol = :completo,
                                  nombres = ("Brincadiano", "GGIR"))
    n = length(datos)
    fig = Figure(size = figure_size(ancho, 50 + 48n))
    Legend(fig[1, 1:2], [PolyElement(color = GRIS), PolyElement(color = (ROJO, 0.3))], ["no uso", "discrepancia"];
           orientation = :horizontal, tellheight = true)
    w = Dates.value(convert(Millisecond, paso)) / 86_400_000
    for (j, (lab, t, a, b)) in enumerate(datos)
        origin = floor(first(t), Day)
        x = days_since(t, origin)
        nd = ceil(Int, maximum(x) + w)
        ax = Axis(fig[j + 1, 2]; yticks = ([1, 2], [nombres[2], nombres[1]]), ytickformat = Makie.automatic,
                  ygridvisible = false, rightspinevisible = true, topspinevisible = true)
        for (i, k) in _bool_runs(a .!= b)
            vspan!(ax, x[i], x[k] + w; color = (ROJO, 0.3))
        end
        for (yy, v) in ((2, a), (1, b))
            for (i, k) in _bool_runs(v)
                poly!(ax, Rect2f(x[i], yy - 0.32, x[k] + w - x[i], 0.64); color = GRIS)
            end
        end
        ylims!(ax, 0.4, 2.6)
        _day_ticks!(ax, origin, nd)
        Label(fig[j + 1, 1], lab; halign = :right, tellheight = false)
    end
    colsize!(fig.layout, 1, Fixed(58))
    rowgap!(fig.layout, 4)
    return fig
end

function _plot_timing_sensitivity(filas; ancho::Symbol = :texto)
    vars = unique(filas.variante); ids = unique(filas.id); nv = length(vars); ni = length(ids)
    fig = Figure(size = figure_size(ancho, 50 + 15 * nv))
    yt = (collect(1.0:nv), vars)
    axθ = Axis(fig[2, 1]; xlabel = "θ (h tras el ancla)", yticks = yt, yreversed = true,
               xticks = (0:3:24, ["$(h)" for h in 0:3:24]), xtickformat = Makie.automatic,
               ytickformat = Makie.automatic, ygridvisible = false)
    axR = Axis(fig[2, 2]; xlabel = "R", yticks = yt, ytickformat = Makie.automatic, yreversed = true, ygridvisible = false,
               yticklabelsvisible = false, yticksvisible = false)
    off(j) = ni == 1 ? 0.0 : 0.36 * ((j - 1) / (ni - 1) - 0.5)
    if hasproperty(filas, :grupo)                    # separadores entre grupos de variantes
        g = [filas.grupo[findfirst(==(v), filas.variante)] for v in vars]
        sep = [k + 0.5 for k in 1:nv-1 if g[k] != g[k+1]]
        for ax in (axθ, axR); hlines!(ax, sep; color = (:black, 0.25), linewidth = 0.5); end
    end
    for (j, id) in enumerate(ids)
        c = COLORES_ID[j]
        f = filas[filas.id .== id, :]
        isempty(f) || vlines!(axθ, [f.θ_h[1]]; color = (c, 0.5), linestyle = :dash, linewidth = 0.8)
        for r in eachrow(f)
            y = findfirst(==(r.variante), vars) + off(j)
            isfinite(r.θ_h) || continue
            a, b = r.θ_ic_desde_h, r.θ_ic_hasta_h
            if isfinite(a) && isfinite(b)
                segs = a <= b ? [(a, b)] : [(a, 24.0), (0.0, b)]          # el arco puede cruzar las 24 h
                for (u, v) in segs
                    lines!(axθ, [u, v], [y, y]; color = c, linewidth = 1.2)
                end
            end
            scatter!(axθ, [r.θ_h], [y]; color = c, markersize = 5)
            isfinite(r.R_ic_inf) && lines!(axR, [r.R_ic_inf, r.R_ic_sup], [y, y]; color = c, linewidth = 1.2)
            scatter!(axR, [r.R], [y]; color = c, markersize = 5)
        end
    end
    xlims!(axθ, 0, 24); xlims!(axR, 0, 1); ylims!(axθ, nv + 0.6, 0.4); ylims!(axR, nv + 0.6, 0.4)
    Legend(fig[1, 1:2], [MarkerElement(marker = :circle, color = COLORES_ID[j], markersize = 6) for j in 1:ni],
           ids; orientation = :horizontal, tellheight = true)
    colsize!(fig.layout, 2, Relative(0.22))
    rowgap!(fig.layout, 2)
    return fig
end

function _plot_rhythm_spectrum(datos; fmax::Real = 8, ancho::Symbol = :completo)
    fig = Figure(size = figure_size(ancho, 190))
    ax1 = Axis(fig[2, 1]; xlabel = "frecuencia (ciclos/día)", ylabel = "varianza explicada",
               xticks = 0:1:fmax)
    ax2 = Axis(fig[2, 2]; xlabel = "frecuencia (ciclos/día)", ylabel = "ventana espectral",
               xticks = 0:1:fmax)
    ax3 = Axis(fig[2, 3]; xlabel = "periodo (h)", ylabel = "Qp / umbral",
               xticks = (0:6:36, string.(0:6:36)), xtickformat = Makie.automatic)
    for ax in (ax1, ax2)
        vlines!(ax, 1:floor(Int, fmax); color = (:black, 0.12), linewidth = 0.6)
    end
    hlines!(ax3, [1.0]; color = (:black, 0.4), linestyle = :dash, linewidth = 0.8)
    for (j, (lab, g, w, q)) in enumerate(datos)
        c = COLORES_ID[j]; k = g.frecuencia .<= fmax
        lines!(ax1, g.frecuencia[k], g.potencia[k]; color = c, linewidth = 1)
        lines!(ax2, g.frecuencia[k], w[k]; color = c, linewidth = 1)
        lines!(ax3, q.periodo_h, q.Qp_rel; color = c, linewidth = 1)
    end
    xlims!(ax1, 0, fmax); xlims!(ax2, 0, fmax); xlims!(ax3, 0, 36)
    ylims!(ax1, 0, nothing); ylims!(ax2, 0, nothing); ylims!(ax3, 0, nothing)
    Legend(fig[1, 1:3], [LineElement(color = COLORES_ID[j], linewidth = 1.5) for j in eachindex(datos)],
           [d[1] for d in datos]; orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 2)
    return fig
end

function _plot_harmonic_profile(datos; K = (1, 4), ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 170))
    hs = range(0, 24; length = 241)
    for (j, (lab, h)) in enumerate(datos)
        c = COLORES_ID[j]
        ax = Axis(fig[1, j]; xlabel = "hora del día", ylabel = j == 1 ? "ENMO (mg)" : "",
                  xticks = _hhmm_ticks(0:6:24), xtickformat = Makie.automatic)
        scatter!(ax, collect(0:23) .+ 0.5, h.perfil; color = (c, 0.8), markersize = 4, label = lab)
        for (i, kk) in enumerate(K)
            a = h.armonicos[h.armonicos.k .<= kk, :]
            rec = [h.media + sum(r.amplitud * cos(2π * r.k * (x - r.acrofase_h) / 24) for r in eachrow(a)) for x in hs]
            pct = 100 * sum(a.IS_k) / h.IS
            nombre = kk == 1 ? "1 armónico" : "$(kk) armónicos"
            lines!(ax, collect(hs), rec; color = c, linewidth = i == 1 ? 0.8 : 1.4, linestyle = i == 1 ? :dash : :solid,
                   label = "$(nombre): $(_coma([round(pct; digits = 0)])[1]) % de IS")
        end
        xlims!(ax, 0, 24); ylims!(ax, 0, nothing)
        axislegend(ax; position = :lt, framevisible = false, labelsize = 7, rowgap = 0, padding = (2, 2, 2, 2))
    end
    return fig
end

const COLORES_CLASE = Dict("sedentaria" => RGBf(0.55, 0.55, 0.55), "ligera" => TURQUESA,
                           "moderada" => AMBAR, "vigorosa" => RGBf(0.80, 0.15, 0.20))

function _plot_psd_classes(datos, fr; suelo = nothing, cuantizacion = nothing, ancho::Symbol = :completo)
    fig = Figure(size = figure_size(ancho, 200))
    k = fr .>= 0.25
    els = Any[]; nombres = String[]
    for (j, (lab, clases)) in enumerate(datos)
        ax = Axis(fig[2, j]; xlabel = "frecuencia (Hz)", ylabel = j == 1 ? "PSD (mg²/Hz)" : "",
                  yscale = log10, xticks = 0:3:15, ytickformat = Makie.automatic,
                  yticks = (10.0 .^ (-1:6), ["10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "10⁵", "10⁶"]))
        for (nombre, p) in clases
            l = lines!(ax, fr[k], max.(p[k], 1e-3); color = COLORES_CLASE[nombre], linewidth = 1.2)
            j == 1 && (push!(els, l); push!(nombres, nombre))
        end
        if suelo !== nothing && suelo[j] !== nothing
            l = lines!(ax, fr[k], suelo[j][k]; color = :black, linestyle = :dash, linewidth = 0.8)
            j == 1 && (push!(els, l); push!(nombres, "1 % de ventanas más quietas"))
        end
        if cuantizacion !== nothing
            l = hlines!(ax, [cuantizacion]; color = :black, linestyle = :dot, linewidth = 0.8)
            j == 1 && (push!(els, l); push!(nombres, "ruido de cuantización"))
        end
        text!(ax, 0.98, 0.97; text = lab, space = :relative, align = (:right, :top), fontsize = 7.5)
        xlims!(ax, 0, 15); ylims!(ax, 0.04, 2e6)
    end
    Legend(fig[1, 1:length(datos)], els, nombres; orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 2)
    return fig
end

function _plot_cadence(datos, fr; ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 170))
    els = Any[]; nombres = String[]
    for (j, (lab, clases)) in enumerate(datos)
        ax = Axis(fig[2, j]; xlabel = "frecuencia dominante (Hz)", ylabel = j == 1 ? "fracción de ventanas" : "",
                  xticks = 0:1:5)
        for (nombre, f) in clases
            isempty(f) && continue
            c = countmap_f(f, fr)
            l = lines!(ax, c[1], c[2]; color = COLORES_CLASE[nombre], linewidth = 1.2)
            j == 1 && (push!(els, l); push!(nombres, nombre))
        end
        text!(ax, 0.02, 0.97; text = lab, space = :relative, align = (:left, :top), fontsize = 7.5)
        xlims!(ax, 0.5, 5); ylims!(ax, 0, nothing)
    end
    Legend(fig[1, 1:length(datos)], els, nombres; orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 2)
    return fig
end
# Fracción de valores en cada bin de frecuencia (las f_dom caen exactamente en `fr`).
function countmap_f(f, fr)
    b = [f_ for f_ in fr if 0.5 <= f_ <= 5]
    n = [count(x -> abs(x - c) < 1e-9, f) for c in b]
    return b, n ./ length(f)
end

function _plot_episode_spectral_qc(datos; umbral::Real = 0.5, ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 170))
    for (j, (lab, ep, ais)) in enumerate(datos)
        ax = Axis(fig[2, j]; xlabel = "frecuencia dominante (Hz)", ylabel = j == 1 ? "potencia en 10–15 Hz" : "",
                  xticks = 0:1:5)
        scatter!(ax, ais.f_dom, ais.frac_3; color = (GRIS, 0.25), markersize = 3, strokewidth = 0)
        scatter!(ax, ep.f_dom, ep.frac_3; color = (COLORES_ID[j], 0.7), markersize = 4, strokewidth = 0)
        hlines!(ax, [umbral]; color = ROJO, linestyle = :dash, linewidth = 0.8)
        text!(ax, 0.02, 0.97; text = lab, space = :relative, align = (:left, :top), fontsize = 7.5)
        xlims!(ax, 0.4, 5.1); ylims!(ax, 0, 1)
    end
    Legend(fig[1, 1:length(datos)],
           [MarkerElement(marker = :circle, color = (AZUL, 0.7), markersize = 6),
            MarkerElement(marker = :circle, color = (GRIS, 0.4), markersize = 6),
            LineElement(color = ROJO, linestyle = :dash)],
           ["episodio vigoroso", "segundo vigoroso aislado", "umbral de vibración"];
           orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 2)
    return fig
end

# ─── Figuras complementarias (H26) ────────────────────────────────────────────

const PARTES_DIA = [("uso", AZUL), ("inmóvil con contacto (uso)", MORADO),
                    ("no uso probable", RGBf(0.25, 0.25, 0.25)), ("sin datos", RGBf(0.85, 0.35, 0.3)),
                    ("fuera del registro", RGBf(0.88, 0.88, 0.88))]

function _plot_qc_days(datos; minimo_h::Real = 10, ancho::Symbol = :completo)
    nmax = maximum(nrow(d) for (_, d) in datos)
    fig = Figure(size = figure_size(ancho, 62 + 15 * nmax))
    for (j, (lab, d)) in enumerate(datos)
        n = nrow(d)
        labs = ["día $(k) ($(Viz.DIAS_SEMANA[dayofweek(d.dia[k])]))" for k in 1:n]
        ax = Axis(fig[2, j]; title = lab, xlabel = "horas del día natural", xticks = 0:4:24,
                  yticks = (collect(1.0:n), labs), ytickformat = Makie.automatic, yreversed = true,
                  ygridvisible = false)
        for (k, r) in enumerate(eachrow(d))
            partes = (r.uso_h - r.inmovil_contacto_h, r.inmovil_contacto_h, r.no_uso_h, r.sin_datos_h,
                      24 - r.cubiertas_h)
            x = 0.0
            for (v, (_, c)) in zip(partes, PARTES_DIA)
                v > 1e-6 && poly!(ax, Rect2f(x, k - 0.36, v, 0.72); color = c)
                x += max(v, 0)
            end
            text!(ax, 24.6, k; text = r.valido ? "válido" : "no válido", align = (:left, :center),
                  fontsize = 7, color = r.valido ? TURQUESA : (:black, 0.5))
        end
        vlines!(ax, [minimo_h]; color = (:black, 0.7), linestyle = :dash, linewidth = 0.8)
        xlims!(ax, 0, 29.5); ylims!(ax, nmax + 0.6, 0.4)
        j > 1 && hideydecorations!(ax; ticklabels = false)
    end
    Legend(fig[1, 1:length(datos)], vcat([PolyElement(color = c) for (_, c) in PARTES_DIA],
                                         [LineElement(color = (:black, 0.7), linestyle = :dash)]),
           vcat(first.(PARTES_DIA), ["mínimo de uso ($(_coma([minimo_h])[1]) h)"]);
           orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_capsense_runs(datos; ratio::Real = 0.25, ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 215))
    xt = [0.25, 0.5, 1, 2, 4, 8, 16]
    for (j, (lab, t, cref)) in enumerate(datos)
        ax = Axis(fig[2, j]; title = lab, xlabel = "duración del tramo (h)",
                  ylabel = j == 1 ? "contacto con la piel (%)" : "", xscale = log10,
                  xticks = (xt, ["0,25", "0,5", "1", "2", "4", "8", "16"]), xtickformat = Makie.automatic,
                  yticks = 0:20:100)
        hlines!(ax, [100cref]; color = (TURQUESA, 0.8), linestyle = :dot, linewidth = 0.9)
        hlines!(ax, [100ratio * cref]; color = (:black, 0.7), linestyle = :dash, linewidth = 0.8)
        for (et, c) in TRAMOS
            for (noct, mk) in ((true, :circle), (false, :utriangle))
                sel = (t.etiqueta .== et) .& ((t.nocturno .>= 0.5) .== noct)
                any(sel) && scatter!(ax, t.horas[sel], 100 .* t.contacto[sel]; color = (c, 0.85),
                                     marker = mk, markersize = 6, strokewidth = 0)
            end
        end
        xlims!(ax, 0.2, 20); ylims!(ax, -3, 103)
    end
    Legend(fig[1, 1:length(datos)],
           vcat(Any[MarkerElement(marker = :circle, color = c, markersize = 6) for (_, c) in TRAMOS[1:2]],
                Any[MarkerElement(marker = :utriangle, color = (:black, 0.5), markersize = 6),
                    LineElement(color = (:black, 0.7), linestyle = :dash), LineElement(color = TURQUESA, linestyle = :dot)]),
           vcat(first.(TRAMOS[1:2]), ["diurno (triángulo)",
                "regla: $(_coma([ratio])[1]) × referencia", "contacto de referencia (en uso)"]);
           orientation = :horizontal, tellheight = true, nbanks = 2)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_anchor_nights(datos; ancho::Symbol = :texto)
    nrows = sum(nrow(n) for (_, n) in datos)
    fig = Figure(size = figure_size(ancho, 58 + 12 * nrows + 8 * length(datos)))
    ax1 = Axis(fig[2, 1]; xlabel = "punto medio del periodo − punto medio de L5 (h)",
               yreversed = true, ygridvisible = false, ytickformat = Makie.automatic)
    ax2 = Axis(fig[2, 2]; xlabel = "% del periodo en idle sleep", yreversed = true,
               ygridvisible = false, yticklabelsvisible = false, yticksvisible = false)
    y = 0.0; ticks = Float64[]; labs = String[]
    for (j, (lab, n)) in enumerate(datos)
        c = COLORES_ID[j]
        for (k, r) in enumerate(eachrow(n))
            y += 1
            push!(ticks, y); push!(labs, "$(lab) · $(k) ($(Viz.DIAS_SEMANA[r.dia_semana]))")
            isfinite(r.spt_menos_L5_h) && scatter!(ax1, [r.spt_menos_L5_h], [y];
                color = r.metodo == "retirada" ? (c, 0.4) : c, markersize = 6,
                marker = r.metodo == "retirada" ? :rect : :circle, strokewidth = 0)
            isfinite(r.pct_spt_idle) && poly!(ax2, Rect2f(0, y - 0.32, r.pct_spt_idle, 0.64); color = (c, 0.8))
        end
        y += 0.8
    end
    ax1.yticks = (ticks, labs)
    vlines!(ax1, [0.0]; color = (:black, 0.5), linewidth = 0.8)
    ylims!(ax1, y + 0.2, 0.2); ylims!(ax2, y + 0.2, 0.2); xlims!(ax2, 0, 100)
    Legend(fig[1, 1:2], [MarkerElement(marker = :circle, color = (:black, 0.7), markersize = 6),
                         MarkerElement(marker = :rect, color = (:black, 0.3), markersize = 6)],
           ["periodo estimado por la postura", "SPT* (retirada del dispositivo)"]; orientation = :horizontal, tellheight = true)
    colsize!(fig.layout, 2, Relative(0.28))
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_timing_days(datos; ancla::AbstractString = "punto medio de L5", ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 240))
    pos(r, φ) = (r * sin(φ), r * cos(φ))               # 0 arriba, sentido horario (como la diana)
    φs = range(0, 2π; length = 361)
    for (j, (lab, d, θh, R, ic)) in enumerate(datos)
        ax = Axis(fig[2, j]; aspect = DataAspect(), title = lab)
        hidedecorations!(ax); hidespines!(ax)
        lines!(ax, [pos(1, φ)[1] for φ in φs], [pos(1, φ)[2] for φ in φs]; color = (AZUL, 0.6), linewidth = 0.7)
        for h in 0:3:21
            φ = 2π * h / 24; a, b = pos(1, φ), pos(0.94, φ); t = pos(1.14, φ)
            lines!(ax, [a[1], b[1]], [a[2], b[2]]; color = (AZUL, 0.6))
            text!(ax, t[1], t[2]; text = h == 0 ? "0 h" : "+$(h) h", align = (:center, :center), fontsize = 7)
        end
        # Intervalo bootstrap de θ de la semana, como arco exterior.
        a, b = ic
        if isfinite(a) && isfinite(b)
            span = mod(b - a, 24)
            arco = range(2π * a / 24, 2π * (a + span) / 24; length = 100)
            lines!(ax, [pos(1.04, φ)[1] for φ in arco], [pos(1.04, φ)[2] for φ in arco];
                   color = (MORADO, 0.55), linewidth = 3)
        end
        for r in eachrow(d)
            isfinite(r.θ_h) || continue
            e = pos(r.R, 2π * r.θ_h / 24)
            c = dayofweek(r.dia) >= 6 ? AMBAR : AZUL
            arrows2d!(ax, [0.0], [0.0], [e[1]], [e[2]]; color = (c, 0.75), shaftwidth = 0.8, tipwidth = 4.5)
        end
        e = pos(R, 2π * θh / 24)
        arrows2d!(ax, [0.0], [0.0], [e[1]], [e[2]]; color = MORADO, shaftwidth = 2, tipwidth = 8)
        xlims!(ax, -1.3, 1.3); ylims!(ax, -1.3, 1.3)
    end
    Legend(fig[1, 1:length(datos)], [LineElement(color = AZUL), LineElement(color = AMBAR),
                                     LineElement(color = MORADO, linewidth = 2), LineElement(color = (MORADO, 0.55), linewidth = 3)],
           ["un día de diario", "un día de fin de semana", "la semana (corregida)", "IC 95 % de θ (bootstrap de días)"];
           orientation = :horizontal, tellheight = true, nbanks = 2)
    _nota!(fig[3, 1:length(datos)], "Ángulo: horas tras el $(ancla). Longitud: R (el círculo es R = 1).")
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_clock_coverage(datos; ancho::Symbol = :completo)
    fig = Figure(size = figure_size(ancho, 230))
    for (j, (lab, cob, ep, nd, ancla_h)) in enumerate(datos)
        ax1 = Axis(fig[2, j]; title = lab, ylabel = j == 1 ? "días observados" : "", xticks = 0:3:24)
        ax2 = Axis(fig[3, j]; xlabel = "hora del reloj", ylabel = j == 1 ? "dosis vigorosa (g·s/día)" : "",
                   xticks = 0:3:24)
        m = [sum(cob[(60i + 1):(60i + 60)]) / 60 for i in 0:1439]
        stairs!(ax1, collect(0:1439) ./ 60, m; step = :post, color = AZUL, linewidth = 0.9)
        hlines!(ax1, [nd]; color = (:black, 0.4), linestyle = :dot, linewidth = 0.8)
        h = hour.(ep.medio)
        sin_c = [sum(ep.dosis[h .== k]; init = 0.0) / 1000 / max(nd, 1) for k in 0:23]
        con_c = [sum(ep.peso[h .== k]; init = 0.0) / 1000 for k in 0:23]
        barplot!(ax2, (0:23) .+ 0.3, sin_c; width = 0.4, color = (GRIS, 0.45), strokewidth = 0)
        barplot!(ax2, (0:23) .+ 0.7, con_c; width = 0.4, color = AMBAR, strokewidth = 0)
        for ax in (ax1, ax2)
            vlines!(ax, [ancla_h]; color = MORADO, linestyle = :dash, linewidth = 0.8)
            xlims!(ax, 0, 24)
        end
        ylims!(ax1, 0, nd + 0.6)
        linkxaxes!(ax1, ax2); hidexdecorations!(ax1; grid = false)
    end
    Legend(fig[1, 1:length(datos)], [LineElement(color = AZUL), LineElement(color = (:black, 0.4), linestyle = :dot),
                                     PolyElement(color = (GRIS, 0.45)), PolyElement(color = AMBAR),
                                     LineElement(color = MORADO, linestyle = :dash)],
           ["días en que se observó esa hora", "días válidos", "dosis por día válido (sin corregir)",
            "dosis por día observado a esa hora (corregida)", "punto medio de L5 (ancla)"];
           orientation = :horizontal, tellheight = true, nbanks = 2)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_threshold_sweep(datos; referencia::Real = 464.6, marcas = (), ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 300))
    ax1 = Axis(fig[2, 1]; ylabel = "Δθ (h)")
    ax2 = Axis(fig[3, 1]; ylabel = "R", yticks = 0:0.2:1)
    ax3 = Axis(fig[4, 1]; ylabel = "episodios", xlabel = "umbral de ENMO de 1 s (mg)", yscale = log10)
    for (j, (lab, b)) in enumerate(datos)
        c = COLORES_ID[j]
        k = argmin(abs.(b.umbral_mg .- referencia))
        dθ = [mod(θ - b.θ_h[k] + 12, 24) - 12 for θ in b.θ_h]
        lines!(ax1, b.umbral_mg, dθ; color = c, linewidth = 1, label = lab)
        lines!(ax2, b.umbral_mg, b.R; color = c, linewidth = 1)
        ok = b.episodios .> 0
        lines!(ax3, b.umbral_mg[ok], Float64.(b.episodios[ok]); color = c, linewidth = 1)
    end
    for ax in (ax1, ax2, ax3)
        vlines!(ax, [referencia]; color = (:black, 0.6), linewidth = 0.8)
        isempty(marcas) || vlines!(ax, collect(marcas); color = (:black, 0.45), linestyle = :dot, linewidth = 0.8)
    end
    hlines!(ax1, [0.0]; color = (:black, 0.3), linewidth = 0.6)
    ylims!(ax2, 0, 1.08)
    linkxaxes!(ax1, ax2, ax3); hidexdecorations!(ax1; grid = false); hidexdecorations!(ax2; grid = false)
    Legend(fig[1, 1], ax1; orientation = :horizontal, tellheight = true)
    _nota!(fig[5, 1], "Línea continua: umbral de referencia; de puntos: las variantes del análisis de sensibilidad.")
    rowgap!(fig.layout, 9)
    return fig
end

function _plot_nights_subsample(datos; variables = (("punto_medio_h", "punto medio del reposo"),
                                                    ("duracion_h", "duración del periodo")),
                                ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 215))
    axs = Axis[]
    for (j, (_, lab)) in enumerate(variables)
        ax = Axis(fig[2, j]; xlabel = "noches usadas", ylabel = j == 1 ? "error absoluto (h)" : "", title = lab,
                  titlesize = 8.5, titlefont = :regular, ygridvisible = true)
        push!(axs, ax)
    end
    for (i, (id, r)) in enumerate(datos), (j, (v, _)) in enumerate(variables)
        d = r[r.variable .== v, :]
        isempty(d) && continue
        c = COLORES_ID[i]
        lines!(axs[j], d.k, d.mediana_abs; color = c, linewidth = 1.1, label = id)
        lines!(axs[j], d.k, d.max_abs; color = c, linewidth = 0.9, linestyle = :dash)
        scatter!(axs[j], d.k, d.mediana_abs; color = c, markersize = 4)
    end
    for ax in axs
        hlines!(ax, [0.5]; color = (:black, 0.3), linewidth = 0.6)
        ylims!(ax, 0, nothing)
    end
    linkyaxes!(axs...)
    n = maximum(maximum(r.k) for (_, r) in datos)
    foreach(ax -> ax.xticks = 1:n, axs)
    Legend(fig[1, 1:length(variables)], axs[1]; orientation = :horizontal, tellheight = true)
    _nota!(fig[3, 1:length(variables)], "Continua: mediana de todos los subconjuntos de k noches; discontinua: el peor. Gris: 0,5 h.")
    rowgap!(fig.layout, 6)
    return fig
end

# ─── Figuras de luz × actividad × sueño (H27) ─────────────────────────────────

# Colores de los seis estados de Light.ESTADOS_LUZ (los mismos que plot_light_days).
const COLOR_ESTADO_LUZ = Dict("válido" => AZUL, "oscuro y activo" => RGBf(0.80, 0.15, 0.20),
                              "oscuro e inactivo" => MORADO, "sin dato" => RGBf(0.93, 0.70, 0.48),
                              "sueño" => RGBf(0.70, 0.87, 0.85), "retirado" => RGBf(0.55, 0.55, 0.55))
const ORDEN_ESTADOS_LUZ = ("válido", "oscuro y activo", "oscuro e inactivo", "sin dato", "sueño", "retirado")
_despierto_puesto(d) = (d.estado .!= "sueño") .& (d.estado .!= "retirado")

# Filas «ID · día k (dow)» de varios registros con un hueco entre registros.
function _filas_dias(datos)
    y = 0.0; filas = Tuple{Float64,Int,Any}[]; ticks = Float64[]; labs = String[]
    for (j, (lab, d)) in enumerate(datos)
        for (k, r) in enumerate(eachrow(d))
            y += 1
            push!(filas, (y, j, r)); push!(ticks, y)
            push!(labs, "$(lab) · día $(k) ($(Viz.DIAS_SEMANA[r.dia_semana]))")
        end
        y += 0.8
    end
    return filas, ticks, labs, y
end

function _plot_light_anchored(datos; ancho::Symbol = :completo)
    filas, ticks, labs, ymax = _filas_dias(datos)
    fig = Figure(size = figure_size(ancho, 62 + 12.5 * length(filas) + 8 * length(datos)))
    ax1 = Axis(fig[2, 1]; xlabel = "hora del día", xticks = _hhmm_ticks(0:3:27), xtickformat = Makie.automatic,
               yticks = (ticks, labs), ytickformat = Makie.automatic, yreversed = true, ygridvisible = false)
    ax2 = Axis(fig[2, 2]; xlabel = "% de luz válida en la ventana", yreversed = true, ygridvisible = false,
               yticklabelsvisible = false, yticksvisible = false)
    for (y, j, r) in filas
        a, b = r.despertar_h, r.dormir_h
        if isfinite(a)
            poly!(ax1, Rect2f(a, y - 0.32, 2, 0.64); color = (TURQUESA, 0.25))
        end
        if isfinite(b)
            poly!(ax1, Rect2f(b - 3, y - 0.32, 3, 0.64); color = (MORADO, 0.2))
        end
        isfinite(a) && isfinite(b) && lines!(ax1, [a, b], [y, y]; color = (:black, 0.45), linewidth = 0.8)
        for (h, mk) in ((r.primera_alta_h, :rtriangle), (r.ultima_alta_h, :ltriangle))
            isfinite(h) && scatter!(ax1, [h], [y]; marker = mk, color = AMBAR, markersize = 6)
        end
        isfinite(r.hora_media_alta_h) && scatter!(ax1, [r.hora_media_alta_h], [y]; marker = :diamond,
                                                   color = AZUL, markersize = 5)
        isfinite(r.tras_despertar_pct_valido) &&
            poly!(ax2, Rect2f(0, y - 0.36, r.tras_despertar_pct_valido, 0.34); color = TURQUESA)
        isfinite(r.antes_dormir_pct_valido) &&
            poly!(ax2, Rect2f(0, y + 0.02, r.antes_dormir_pct_valido, 0.34); color = MORADO)
    end
    xlims!(ax1, 0, 27); xlims!(ax2, 0, 100); ylims!(ax1, ymax + 0.2, 0.2); ylims!(ax2, ymax + 0.2, 0.2)
    Legend(fig[1, 1:2], [LineElement(color = (:black, 0.45)), PolyElement(color = (TURQUESA, 0.5)),
                         PolyElement(color = (MORADO, 0.4)),
                         MarkerElement(marker = :rtriangle, color = AMBAR, markersize = 6),
                         MarkerElement(marker = :ltriangle, color = AMBAR, markersize = 6),
                         MarkerElement(marker = :diamond, color = AZUL, markersize = 5)],
           ["despierto (fin e inicio del reposo nocturno)", "2 h tras despertar", "3 h antes de dormir",
            "primera luz alta", "última luz alta", "hora media de luz alta"];
           orientation = :horizontal, tellheight = true, nbanks = 2)
    colsize!(fig.layout, 2, Relative(0.2))
    rowgap!(fig.layout, 3)
    return fig
end

const PARTES_AFMV = [("a oscuras", RGBf(0.80, 0.15, 0.20)), ("interior", AZUL),
                     ("compatible con exterior", AMBAR), ("sin luz válida", RGBf(0.82, 0.82, 0.82))]

function _plot_mvpa_light(datos; ancho::Symbol = :texto)
    filas, ticks, labs, ymax = _filas_dias(datos)
    fig = Figure(size = figure_size(ancho, 58 + 12 * length(filas) + 8 * length(datos)))
    ax = Axis(fig[2, 1]; xlabel = "minutos de AFMV del día", yticks = (ticks, labs), ytickformat = Makie.automatic,
              yreversed = true, ygridvisible = false)
    for (y, _, r) in filas
        x = 0.0
        for (v, (_, c)) in zip((r.afmv_oscuro, r.afmv_interior, r.afmv_exterior, r.afmv_sin_luz), PARTES_AFMV)
            v > 0 && poly!(ax, Rect2f(x, y - 0.34, v, 0.68); color = c)
            x += v
        end
    end
    ylims!(ax, ymax + 0.2, 0.2); xlims!(ax, 0, nothing)
    Legend(fig[1, 1], [PolyElement(color = c) for (_, c) in PARTES_AFMV], first.(PARTES_AFMV);
           orientation = :horizontal, tellheight = true)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_light_enmo_density(datos; dark::Real = 1, active_mg::Real = 63.3, ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 215))
    xe = range(0, 3.3; length = 45); ye = range(0, 3.6; length = 37)
    xt = [0, 10, 63.3, 300, 1000]; yt = [0, 1, 10, 100, 1000]
    mapa = cgrad([RGBf(0.93, 0.95, 0.97), AZUL])
    hm = nothing; cmax = 0.0
    mats = map(datos) do (lab, d)
        sel = _despierto_puesto(d) .& (d.lux_obs .>= 0.5) .& isfinite.(d.enmo_idle0) .& isfinite.(d.lux_media)
        x = log10.(1 .+ max.(d.enmo_idle0[sel], 0)); y = log10.(1 .+ max.(d.lux_media[sel], 0))
        m = zeros(length(xe) - 1, length(ye) - 1)
        for (u, v) in zip(x, y)
            i = clamp(searchsortedlast(xe, u), 1, length(xe) - 1); k = clamp(searchsortedlast(ye, v), 1, length(ye) - 1)
            m[i, k] += 1
        end
        lm = [c > 0 ? log10(c) : NaN for c in m]
        cmax = max(cmax, maximum(filter(isfinite, lm); init = 0.0))
        (lab, lm)
    end
    for (j, (lab, lm)) in enumerate(mats)
        ax = Axis(fig[1, j]; title = lab, xlabel = "ENMO de 60 s (mg)", ylabel = j == 1 ? "luz media (bruta)" : "",
                  xticks = (log10.(1 .+ xt), ["0", "10", "63", "300", "1000"]), xtickformat = Makie.automatic,
                  yticks = (log10.(1 .+ yt), ["0", "1", "10", "100", "1000"]), ytickformat = Makie.automatic)
        hm = heatmap!(ax, collect(xe), collect(ye), lm; colormap = mapa, colorrange = (0, cmax), nan_color = :transparent)
        vlines!(ax, [log10(1 + active_mg)]; color = (:black, 0.6), linestyle = :dash, linewidth = 0.8)
        hlines!(ax, [log10(1 + dark)]; color = (:black, 0.6), linestyle = :dash, linewidth = 0.8)
        text!(ax, 3.25, 0.13; text = "oscuro y activo", align = (:right, :bottom), fontsize = 7, color = RGBf(0.80, 0.15, 0.20))
        text!(ax, 0.05, 0.13; text = "oscuro e inactivo", align = (:left, :bottom), fontsize = 7, color = MORADO)
        xlims!(ax, 0, 3.3); ylims!(ax, 0, 3.6)
    end
    Colorbar(fig[1, length(datos) + 1], hm; label = "minutos (log₁₀)", width = 7,
             ticks = 0:1:ceil(Int, cmax), tickformat = _coma)
    return fig
end

function _plot_light_states_by_hour(datos; ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 200))
    for (j, (lab, d)) in enumerate(datos)
        ax = Axis(fig[2, j]; title = lab, xlabel = "hora del reloj", ylabel = j == 1 ? "% de los minutos" : "",
                  xticks = 0:3:24, yticks = 0:25:100)
        h = hour.(d.start)
        for k in 0:23
            sel = h .== k; n = count(sel)
            n == 0 && continue
            y = 0.0
            for e in ORDEN_ESTADOS_LUZ
                f = 100 * count(d.estado[sel] .== e) / n
                f > 0 && poly!(ax, Rect2f(k, y, 1, f); color = COLOR_ESTADO_LUZ[e])
                y += f
            end
        end
        xlims!(ax, 0, 24); ylims!(ax, 0, 100)
        j > 1 && hideydecorations!(ax; ticks = false, ticklabels = false, grid = false)
    end
    Legend(fig[1, 1:length(datos)], [PolyElement(color = COLOR_ESTADO_LUZ[e]) for e in ORDEN_ESTADOS_LUZ],
           [e == "sueño" ? "reposo nocturno" : e for e in ORDEN_ESTADOS_LUZ]; orientation = :horizontal,
           tellheight = true, nbanks = 2)
    rowgap!(fig.layout, 3)
    return fig
end

function _plot_light_profile(datos; min_validos::Integer = 10, ancho::Symbol = :texto)
    fig = Figure(size = figure_size(ancho, 200))
    yt = [0, 1, 10, 100, 1000]
    for (j, (lab, d)) in enumerate(datos)
        ax = Axis(fig[2, j]; title = lab, xlabel = "hora del reloj", ylabel = j == 1 ? "luz válida (bruta)" : "",
                  xticks = 0:3:24, yticks = (log10.(1 .+ yt), ["0", "1", "10", "100", "1000"]),
                  ytickformat = Makie.automatic)
        axr = Axis(fig[2, j]; yaxisposition = :right, ylabel = j == length(datos) ? "ENMO medio (mg)" : "",
                   ylabelcolor = AMBAR, yticklabelcolor = AMBAR)
        hidexdecorations!(axr); hidespines!(axr); axr.ygridvisible = false
        h = hour.(d.start); dp = _despierto_puesto(d)
        for k in 0:23
            v = d.lux_media[(h .== k) .& (d.estado .== "válido")]
            length(v) >= min_validos || continue
            q = sort(log10.(1 .+ v)); n = length(q); qq(p) = q[clamp(ceil(Int, p * n), 1, n)]
            rangebars!(ax, [k + 0.5], [qq(0.25)], [qq(0.75)]; color = (AZUL, 0.6), linewidth = 1.2)
            scatter!(ax, [k + 0.5], [qq(0.5)]; color = AZUL, markersize = 5)
        end
        en = [begin w = d.enmo_idle0[(h .== k) .& dp .& isfinite.(d.enmo_idle0)]; isempty(w) ? NaN : sum(w) / length(w) end
              for k in 0:23]
        stairs!(axr, collect(0:24), vcat(en, en[end]); step = :post, color = (AMBAR, 0.9), linewidth = 1)
        xlims!(ax, 0, 24); xlims!(axr, 0, 24); ylims!(ax, 0, 3.6); ylims!(axr, 0, nothing)
        linkxaxes!(ax, axr)
    end
    Legend(fig[1, 1:length(datos)], [MarkerElement(marker = :circle, color = AZUL, markersize = 5),
                                     LineElement(color = (AZUL, 0.6), linewidth = 1.2), LineElement(color = AMBAR)],
           ["mediana de la luz válida", "rango intercuartílico", "ENMO medio despierto (eje derecho)"];
           orientation = :horizontal, tellheight = true)
    _nota!(fig[3, 1:length(datos)], "Solo horas con al menos $(min_validos) minutos de luz válida.")
    rowgap!(fig.layout, 3)
    return fig
end

# Arco corto respecto a una referencia, para que 23 h y 1 h queden a 2 h y no a 22.
_arco(h, ref) = ref + (mod(h - ref + 12, 24) - 12)

function _plot_imputation_dumbbell(filas; ancho::Symbol = :completo)
    vars = ["huecos", "ceros", "última"]
    marca = Dict("huecos" => :circle, "ceros" => :rect, "última" => :diamond)
    color = Dict("huecos" => AZUL, "ceros" => GRIS, "última" => AMBAR)
    ids = unique(filas.id)
    fig = Figure(size = figure_size(ancho, 52 + 36length(ids)))
    paneles = ((:ancla_h, "punto medio del ancla (h del reloj)"), (:θ_h, "θ (h tras el ancla)"))
    for (k, (col, xlab)) in enumerate(paneles)
        ax = Axis(fig[1, k]; xlabel = xlab,
                  yticks = (collect(Float64, eachindex(ids)), String.(ids)),
                  ytickformat = Makie.automatic, yreversed = true, ygridvisible = false)
        for (j, id) in enumerate(ids)
            sub = filas[filas.id .== id, :]
            ref = sub[sub.variante .== "huecos", col]
            href = isempty(ref) ? 0.0 : ref[1]
            xs = Float64[]
            cols = []
            for v in vars
                r = sub[sub.variante .== v, :]
                isempty(r) && continue
                push!(xs, _arco(r[1, col], href))
                push!(cols, v)
            end
            length(xs) >= 2 && lines!(ax, xs, fill(j, length(xs)); color = (GRIS, 0.4), linewidth = 1.1)
            for (x, v) in zip(xs, cols)
                scatter!(ax, [x], [j]; marker = marca[v], color = color[v], markersize = 9,
                         strokecolor = :white, strokewidth = 0.6)
            end
        end
        ylims!(ax, 0.4, length(ids) + 0.6)
    end
    Legend(fig[2, :], [MarkerElement(marker = marca[v], color = color[v], markersize = 8) for v in vars],
           ["sin imputar", "ceros", "última muestra"]; orientation = :horizontal, tellheight = true)
    _nota!(fig[3, :], "Cada punto es una variante de la misma señal. La línea une el arco corto, no el salto a través de la medianoche.")
    return fig
end

function _plot_intensity_ggir(filas; ancho::Symbol = :texto)
    ids = unique(filas.id)
    bandas = unique(filas.banda)
    fig = Figure(size = figure_size(ancho, 36 + 28length(bandas)))
    for (j, id) in enumerate(ids)
        ax = Axis(fig[1, j]; xlabel = "minutos al día", title = id,
                  yticks = (collect(Float64, eachindex(bandas)), String.(bandas)),
                  ytickformat = Makie.automatic, yreversed = true, ygridvisible = false)
        for (i, b) in enumerate(bandas)
            r = filas[(filas.id .== id) .& (filas.banda .== b), :]
            isempty(r) && continue
            a, g = r.nuestro[1], r.ggir[1]
            lines!(ax, [a, g], [i, i]; color = (GRIS, 0.45), linewidth = 1.2)
            scatter!(ax, [a], [i]; color = AZUL, markersize = 8, strokecolor = :white, strokewidth = 0.5)
            scatter!(ax, [g], [i]; color = AMBAR, markersize = 8, marker = :rect, strokecolor = :white, strokewidth = 0.5)
        end
        ylims!(ax, 0.4, length(bandas) + 0.6)
    end
    Legend(fig[2, 1:length(ids)],
           [MarkerElement(marker = :circle, color = AZUL, markersize = 8),
            MarkerElement(marker = :rect, color = AMBAR, markersize = 8)],
           ["Brincadiano", "GGIR"]; orientation = :horizontal, tellheight = true)
    _nota!(fig[3, 1:length(ids)], "Media de los días comparados. La distancia es la diferencia de minutos.")
    return fig
end

# Funciones públicas: cada una dibuja con el tema del informe.
for f in (:plot_overview, :plot_raw_window, :plot_days24h, :plot_histograms, :plot_metrics_window,
          :plot_profile24h, :plot_timing_circular, :plot_sleep_nights, :plot_light_days,
          :plot_actogram, :plot_profiles_compare, :plot_enmo_ecdf, :plot_timing_compare, :plot_sleep_raster,
          :plot_vigorous_summary, :plot_circular_plane, :plot_bland_altman, :plot_agreement_timeline,
          :plot_timing_sensitivity, :plot_rhythm_spectrum, :plot_harmonic_profile,
          :plot_psd_classes, :plot_cadence, :plot_episode_spectral_qc,
          :plot_qc_days, :plot_capsense_runs, :plot_anchor_nights, :plot_timing_days,
          :plot_clock_coverage, :plot_threshold_sweep, :plot_light_anchored, :plot_mvpa_light,
          :plot_light_enmo_density, :plot_light_states_by_hour, :plot_light_profile,
          :plot_imputation_dumbbell, :plot_intensity_ggir, :plot_nights_subsample)
    impl = Symbol(:_, f)
    @eval Viz.$f(args...; kw...) = with_theme(() -> $impl(args...; kw...), Viz.tema_informe())
end

# Tramos (i, j) de `true` consecutivos.
function _bool_runs(v)
    r = Tuple{Int,Int}[]; i = 1; N = length(v)
    while i <= N
        if v[i]
            j = i
            while j < N && v[j + 1]; j += 1; end
            push!(r, (i, j)); i = j + 1
        else
            i += 1
        end
    end
    return r
end

end # module
