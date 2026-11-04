# 93 · Validación del periodo de reposo nocturno frente a GGIR en cadera.
#
# Qué hace
#   Lanza GGIR 3.3.9 dos veces (scripts/R/ggir_sueno.R, partes 1–4, sensor en cadera):
#   una con el eje longitudinal que elige GGIR y otra con el eje Y forzado (--eje=2).
#   Compara con los ángulos que guardó el script 08 (solo epochs de muestreo completo
#   y en uso) y, noche a noche, inicio, fin, duración, solapamiento (Jaccard) y sueño
#   dentro de la ventana. Las noches se reconstruyen aquí con sleep_nights; no se lee
#   sueno_<ID>_noches.csv, así que hace falta que angulos_5s.arrow sea el del 08.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/93_validar_sueno_ggir.jl --confirmo
#   julia --project=scripts scripts/93_validar_sueno_ggir.jl --confirmo --sin-ggir
#   Sin --confirmo solo describe lo que haría. --sin-ggir reutiliza
#   outputs/validation/ggir_sueno/ y ggir_sueno_ejeY/ si ya existen. El 96 llama al
#   mismo R con --parte5 y escribe en otro directorio: no comparte esta extracción.
#
# Depende de
#   01 y 08: data/interim/<ID>/angulos_5s.arrow (y el registro, para el inicio).
#   data/raw/*.gt3x, que recibe GGIR.
#   scripts/R/ggir_sueno.R y el renv de scripts/R/.
#
# Lo usan
#   10, que lee ggir_sueno_noches.csv (filas config = ejeY) para la figura de ventanas
#   frente a GGIR. Capítulo 14, vía `make sync-assets`.
#
# Entrada
#   data/raw/*.gt3x
#   data/interim/<ID>/angulos_5s.arrow
#
# Salida de la lista blanca de GGIR (no se versiona; el temporal se borra)
#   outputs/validation/ggir_sueno/ggir_<ID>_angulos5s.csv       eje elegido por GGIR
#   outputs/validation/ggir_sueno/ggir_part4_noches.csv
#   outputs/validation/ggir_sueno_ejeY/                         la misma pareja, con eje Y
# Resumen (no se versiona)
#   outputs/validation/tables/ggir_sueno_angulos.csv
#   outputs/validation/tables/ggir_sueno_noches.csv             una fila por noche y configuración
#   outputs/validation/tables/ggir_sueno_resumen.csv
#   outputs/validation/tables/ggir_sueno.tex
#
# La fecha de calendario de cada noche se convierte a «días desde el inicio» antes de
# escribirla. Consola: solo el ID y los estadísticos (correlación de ángulos, Jaccard).

using Brincadiano.GT3X, Brincadiano.Sleep
using Brincadiano.Rutas: Rutas, raw_files
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const INTERIM = joinpath(ROOT, "data", "interim")
const VAL = joinpath(ROOT, "outputs", "validation")
const OUT = joinpath(VAL, "tables")
const RENV = joinpath(ROOT, "scripts", "R")
const CONFIGS = (("auto", "ggir_sueno", String[]), ("ejeY", "ggir_sueno_ejeY", ["--eje=2"]))
const EJES = ("X", "Y", "Z")

fmt(x; d = 2) = (x isa Real && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
mean_(v) = isempty(v) ? NaN : sum(v) / length(v)
function pearson(a, b)
    ma, mb = mean_(a), mean_(b)
    return sum((a .- ma) .* (b .- mb)) / sqrt(sum((a .- ma) .^ 2) * sum((b .- mb) .^ 2))
end
ggtime(s) = DateTime(first(s, 19))
# Índice de Jaccard de dos intervalos [a0, a1] y [b0, b1] (horas).
function jaccard(a0, a1, b0, b1)
    (isnan(a0) || isnan(b0)) && return NaN
    inter = max(0.0, min(a1, b1) - max(a0, b0))
    uni = (a1 - a0) + (b1 - b0) - inter
    return uni > 0 ? inter / uni : NaN
end

function correr_ggir(args, dir)
    files = raw_files("gt3x")
    err = tempname()
    try
        run(pipeline(Cmd(`Rscript ggir_sueno.R $args $dir $files`; dir = RENV); stdout = devnull, stderr = err))
    catch
        foreach(l -> occursin(r"Error"i, l) && println("  ", l), eachline(err))
        rethrow()
    finally
        rm(err; force = true)
    end
end

function comparar(id, p4s)
    dir = joinpath(INTERIM, id)
    r = read_interim(dir)
    epf = load_angles(dir)
    eje = longitudinal_axis(epf; mask = .!epf.nonwear)
    nn = sleep_nights(epf; axis = eje, nonwear = epf.nonwear, rec_start = first(r.acc.time), rec_stop = last(r.acc.time))
    nn.dia_desde_inicio = [Dates.value(d - Date(first(r.acc.time))) for d in nn.noche]

    # 1. Ángulos (configuración automática; los ángulos no dependen del eje elegido).
    ga = CSV.read(joinpath(VAL, "ggir_sueno", "ggir_$(id)_angulos5s.csv"), DataFrame)
    ga.start = ggtime.(String.(ga.timestamp))
    j = innerjoin(select(epf, :start, :n, :nonwear, :anglex, :angley, :anglez),
                  select(ga, :start, :anglex => :gx, :angley => :gy, :anglez => :gz), on = :start)
    ok = (j.n .== maximum(j.n)) .& .!j.nonwear
    ang = (id = id, epochs = count(ok),
           r_x = pearson(j.anglex[ok], j.gx[ok]), r_y = pearson(j.angley[ok], j.gy[ok]), r_z = pearson(j.anglez[ok], j.gz[ok]),
           mae_x = mean_(abs.(j.anglex[ok] .- j.gx[ok])), mae_y = mean_(abs.(j.angley[ok] .- j.gy[ok])),
           mae_z = mean_(abs.(j.anglez[ok] .- j.gz[ok])))

    # 2. Noches, para cada configuración de GGIR.
    filas = NamedTuple[]
    for (nombre, p4) in p4s
        g = p4[p4.id .== id, :]
        for n in eachrow(nn)
            k = findfirst(==(n.dia_desde_inicio), g.dia_desde_inicio)
            k === nothing && continue
            # Con HorAngle, el tramo «tumbado» de GGIR (equivalente a nuestra ventana) es el
            # guider (tiempo en cama); sleeponset/wakeup son el primer y el último SIB dentro.
            gi, gf = Float64(g.guider_inbedStart[k]), Float64(g.guider_inbedEnd[k])
            gi < 12 && (gi += 24); gf < 12 && (gf += 24)     # horas desde la medianoche del día de la noche
            push!(filas, (id = id, config = nombre, eje_ggir = EJES[g.longitudinal_axis[k]], noche = n.dia_desde_inicio + 1,
                          metodo = n.metodo, inicio_nuestro = n.inicio_h, inicio_ggir = gi,
                          fin_nuestro = n.fin_h, fin_ggir = gf, dur_nuestra = n.duracion_h, dur_ggir = Float64(g.guider_inbedDuration[k]),
                          inicio_sueno_ggir = Float64(g.sleeponset[k]), fin_sueno_ggir = Float64(g.wakeup[k]),
                          sueno_nuestro = n.sueno_h, sueno_ggir = Float64(g.SleepDurationInSpt[k]),
                          jaccard = jaccard(n.inicio_h, n.fin_h, gi, gf),
                          invalida_ggir = Float64(g.fraction_night_invalid[k])))
        end
    end
    return ang, DataFrame(filas), EJES[eje]
end

function resumen(noches)
    combine(groupby(noches, [:id, :config, :eje_ggir]),
            nrow => :noches,
            [:inicio_nuestro, :inicio_ggir] => ((a, b) -> mean_(abs.(a .- b))) => :dif_inicio_h,
            [:fin_nuestro, :fin_ggir] => ((a, b) -> mean_(abs.(a .- b))) => :dif_fin_h,
            [:dur_nuestra, :dur_ggir] => ((a, b) -> mean_(a .- b)) => :dif_dur_h,
            :jaccard => (x -> mean_(filter(isfinite, x))) => :jaccard_medio,
            [:sueno_nuestro, :sueno_ggir] => ((a, b) -> (ok = isfinite.(a); mean_(a[ok] .- b[ok]))) => :dif_sueno_h,
            :invalida_ggir => mean_ => :invalida_ggir)
end

function tabla_tex(ang, res)
    io = IOBuffer()
    ids = ang.id
    println(io, "% Generado por scripts/93_validar_sueno_ggir.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X *{$(length(ids))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(replace.(ids, "_" => "\\_"), " & "), " \\\\")
    println(io, "\\midrule")
    fila(n, v) = println(io, n, " & ", join(v, " & "), " \\\\")
    println(io, "\\multicolumn{$(length(ids) + 1)}{@{}l}{\\textit{Ángulos por \\emph{epoch} de 5\\,s (muestreo completo, en uso)}} \\\\")
    fila("\\emph{Epochs} comparados", string.(ang.epochs))
    fila("Correlación x / y / z", [fmt(a.r_x; d = 4) * " / " * fmt(a.r_y; d = 4) * " / " * fmt(a.r_z; d = 4) for a in eachrow(ang)])
    fila("Diferencia absoluta media x / y / z (grados)", [fmt(a.mae_x) * " / " * fmt(a.mae_y) * " / " * fmt(a.mae_z) for a in eachrow(ang)])
    for (cfg, titulo) in (("auto", "GGIR con el eje que elige él"), ("ejeY", "GGIR con el eje Y forzado"))
        println(io, "\\multicolumn{$(length(ids) + 1)}{@{}l}{\\textit{Periodo de reposo nocturno: $titulo}} \\\\")
        g(id) = (k = findfirst((res.id .== id) .& (res.config .== cfg)); k === nothing ? nothing : res[k, :])
        fila("Eje longitudinal de GGIR / noches", [(x = g(id); x === nothing ? "---" : "$(x.eje_ggir) / $(x.noches)") for id in ids])
        fila("Diferencia absoluta media de inicio / fin (h)", [(x = g(id); x === nothing ? "---" : fmt(x.dif_inicio_h; d = 1) * " / " * fmt(x.dif_fin_h; d = 1)) for id in ids])
        fila("Solapamiento medio (Jaccard)", [(x = g(id); x === nothing ? "---" : fmt(x.jaccard_medio)) for id in ids])
        fila("Duración: nuestra \$-\$ GGIR (h)", [(x = g(id); x === nothing ? "---" : fmt(x.dif_dur_h; d = 1)) for id in ids])
        fila("Sueño en el periodo: nuestro \$-\$ GGIR (h)", [(x = g(id); x === nothing ? "---" : fmt(x.dif_sueno_h; d = 1)) for id in ids])
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

if !("--confirmo" in ARGS)
    println("""Validación externa del periodo de reposo nocturno frente a GGIR (hito H10).
    No se ejecuta sin autorización: añade --confirmo. Lanzaría GGIR (partes 1–4, cadera) dos veces
    sobre data/raw/*.gt3x y compararía ángulos y ventanas con scripts/08_sueno.jl.""")
    exit(0)
end
if !("--sin-ggir" in ARGS && all(isfile(joinpath(VAL, d, "ggir_part4_noches.csv")) for (_, d, _) in CONFIGS))
    for (nombre, d, args) in CONFIGS
        println("GGIR (sueño, $nombre) …")
        correr_ggir(args, joinpath(VAL, d))
    end
else
    println("GGIR: se reutiliza la salida extraída en outputs/validation/ggir_sueno*")
end
p4s = [(nombre, CSV.read(joinpath(VAL, d, "ggir_part4_noches.csv"), DataFrame)) for (nombre, d, _) in CONFIGS]
ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "angulos_5s.arrow")), readdir(INTERIM)))
ids = filter(id -> any(any(p.id .== id) for (_, p) in p4s), ids)
res = map(id -> comparar(id, p4s), ids)
ang = DataFrame(first.(res)); noches = vcat(getindex.(res, 2)...)
rs = resumen(noches)
mkpath(OUT)
CSV.write(joinpath(OUT, "ggir_sueno_angulos.csv"), ang)
CSV.write(joinpath(OUT, "ggir_sueno_noches.csv"), noches)
CSV.write(joinpath(OUT, "ggir_sueno_resumen.csv"), rs)
write(joinpath(OUT, "ggir_sueno.tex"), tabla_tex(ang, rs))
for (a, (_, _, eje)) in zip(eachrow(ang), res)
    @printf("%-12s ángulos: r = %.4f / %.4f / %.4f; |dif| media %.2f / %.2f / %.2f grados (%d epochs); nuestro eje: %s\n",
            a.id, a.r_x, a.r_y, a.r_z, a.mae_x, a.mae_y, a.mae_z, a.epochs, eje)
end
for x in eachrow(rs)
    @printf("%-12s GGIR %-4s (eje %s, %d noches): |Δinicio| %.1f h, |Δfin| %.1f h, Jaccard %.2f, Δduración %+.1f h, Δsueño %+.1f h, noche inválida GGIR %.0f %%\n",
            x.id, x.config, x.eje_ggir, x.noches, x.dif_inicio_h, x.dif_fin_h, x.jaccard_medio, x.dif_dur_h, x.dif_sueno_h, 100x.invalida_ggir)
end
println("Tablas en ", relpath(OUT, ROOT))
