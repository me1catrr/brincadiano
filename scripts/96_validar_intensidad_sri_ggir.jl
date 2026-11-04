# 96 · Tiempo por intensidad y regularidad del sueño frente a GGIR.
#
# Qué hace
#   Lanza GGIR 3.3.9 con partes 1–5 (scripts/R/ggir_sueno.R --eje=2 --parte5): cadera,
#   eje longitudinal Y, ventanas de medianoche a medianoche y umbrales de Hildebrand
#   63,3 / 142,6 / 464,6 mg. Dos contrastes:
#   R2. Minutos por intensidad en el tiempo despierto de cada día
#       (Metrics.intensity_time sobre epochs de 5 s fuera de todo periodo de reposo nocturno)
#       frente a dur_day_total_{IN,LIG,MOD,VIG}_min de la parte 5.
#   R3. Índice de regularidad del sueño (Sleep.sleep_regularity_index sobre los bloques
#       de inactividad sostenida de todo el registro; no uso o sin ángulo = sin dato)
#       frente a SleepRegularityIndex1, y eficiencia dentro de la ventana de postura
#       frente a sleepefficiency. Si GGIR trae la eficiencia en porcentaje, se pasa a fracción.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/96_validar_intensidad_sri_ggir.jl --confirmo
#   julia --project=scripts scripts/96_validar_intensidad_sri_ggir.jl --confirmo --sin-ggir
#   Sin --confirmo termina tras un aviso. --sin-ggir reutiliza
#   outputs/validation/ggir_parte5/ggir_part5_dias.csv. No reutiliza la extracción del
#   93: aquella no tiene parte 5 y vive en ggir_sueno/.
#
# Depende de
#   05  epochs_5s.arrow (bandas de intensidad y tiempo despierto)
#   08  angulos_5s.arrow (ventanas, SIB y eficiencia). Las noches se reconstruyen aquí.
#   01  el registro, para el origen de «días desde el inicio»
#   data/raw/*.gt3x y scripts/R/ggir_sueno.R (el mismo que el 93, con otras opciones)
#
# Lo usan
#   El informe, tras `make sync-assets`. El 10 no lee estas tablas.
#
# Entrada
#   data/raw/*.gt3x
#   data/interim/<ID>/epochs_5s.arrow y angulos_5s.arrow
#
# Salida de la lista blanca de GGIR (no se versiona)
#   outputs/validation/ggir_parte5/ggir_<ID>_angulos5s.csv
#   outputs/validation/ggir_parte5/ggir_part4_noches.csv     incluye SRI y eficiencia
#   outputs/validation/ggir_parte5/ggir_part5_dias.csv      minutos por intensidad
# Resumen (no se versiona)
#   outputs/validation/tables/ggir_intensidad_sri.csv       una fila por registro
#   outputs/validation/tables/ggir_intensidad_sri.tex
#   outputs/validation/tables/ggir_intensidad_dias.csv      una fila por registro y día
#
# Consola: solo el ID y las diferencias en minutos, SRI y eficiencia.

using Brincadiano.GT3X, Brincadiano.Metrics, Brincadiano.Sleep
using Brincadiano.Rutas: Rutas, raw_files
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()          # configuración central (Brincadiano.Rutas)
const INTERIM = joinpath(ROOT, "data", "interim")
const GG = joinpath(ROOT, "outputs", "validation", "ggir_parte5")
const OUT = joinpath(ROOT, "outputs", "validation", "tables")
const RENV = joinpath(ROOT, "scripts", "R")
const BANDAS = (("sedentario", "IN"), ("ligera", "LIG"), ("moderada", "MOD"), ("vigorosa", "VIG"))

fmt(x; d = 1) = (x isa Real && isnan(x)) ? "---" : replace(@sprintf("%.*f", d, round(x; digits = d) == 0 ? 0.0 : x), "." => "{,}", "-" => "−")   # sin «−0,0»
mean_(v) = isempty(v) ? NaN : sum(v) / length(v)
num(v) = [ismissing(x) ? NaN : Float64(x) for x in v]

function correr_ggir()
    files = raw_files("gt3x")
    err = tempname()
    try
        run(pipeline(Cmd(`Rscript ggir_sueno.R --eje=2 --parte5 $GG $files`; dir = RENV); stdout = devnull, stderr = err))
    catch
        foreach(l -> occursin(r"Error"i, l) && println("  ", l), eachline(err))
        rethrow()
    finally
        rm(err; force = true)
    end
end

function comparar(id, p4, p5)
    dir = joinpath(INTERIM, id)
    r = read_interim(dir)
    epf = load_angles(dir); e5 = load_epochs(dir, 5)
    eje = longitudinal_axis(epf; mask = .!epf.nonwear)
    noches = sleep_nights(epf; axis = eje, nonwear = epf.nonwear, rec_start = first(r.acc.time), rec_stop = last(r.acc.time))
    o = Date(first(r.acc.time))
    # R2: despierto = fuera de todo periodo de reposo nocturno (SPT o SPT*).
    iv = [(n.inicio, n.fin) for n in eachrow(noches) if n.metodo != "sin ventana"]
    desp = [!any(a <= t < b for (a, b) in iv) for t in e5.start]
    it = intensity_time(e5, desp)
    it.dia_desde_inicio = [Dates.value(d - o) for d in it.dia]
    g5 = p5[p5.id .== id, :]
    j = innerjoin(it, g5, on = :dia_desde_inicio)
    dias = DataFrame(id = fill(id, nrow(j)), dia = j.dia_desde_inicio .+ 1,
                     despierto_nuestro = j.min_despierto, despierto_ggir = num(j.dur_day_min))
    for (b, g) in BANDAS
        dias[!, b * "_nuestro"] = j[!, b]; dias[!, b * "_ggir"] = num(j[!, "dur_day_total_$(g)_min"])
    end
    # R3: SRI sobre los SIB de todo el registro; no uso o sin ángulo = sin dato.
    sib = sib_vanhees2015(epf.anglez)
    est = [(epf.nonwear[i] || isnan(epf.anglez[i])) ? missing : sib[i] for i in 1:nrow(epf)]
    sri = sleep_regularity_index(epf.start, est)
    sri.dia_desde_inicio = [Dates.value(d - o) for d in sri.dia]
    g4 = p4[p4.id .== id, :]
    js = innerjoin(select(sri, :dia_desde_inicio, :SRI, :frac_valida),
                   select(g4, :dia_desde_inicio, :SleepRegularityIndex1 => :SRI_ggir, :SriFractionValid => :frac_ggir,
                          :sleepefficiency => :efic_ggir), on = :dia_desde_inicio)
    sri_g = num(js.SRI_ggir); okk = isfinite.(sri_g) .& isfinite.(js.SRI)
    noches.dia_desde_inicio = [Dates.value(d - o) for d in noches.noche]
    je = innerjoin(select(noches[noches.metodo .== "postura", :], :dia_desde_inicio, :eficiencia), select(g4, :dia_desde_inicio, :sleepefficiency), on = :dia_desde_inicio)
    ef_g = num(je.sleepefficiency); oke = isfinite.(ef_g) .& isfinite.(je.eficiencia)
    # GGIR expresa la eficiencia en % o en fracción según la versión: se lleva a fracción.
    any(ef_g[oke] .> 1.5) && (ef_g = ef_g ./ 100)
    res = (id = id, dias_intensidad = nrow(j),
           (Symbol("dif_" * b) => mean_(dias[!, b * "_nuestro"] .- dias[!, b * "_ggir"]) for (b, _) in BANDAS)...,
           (Symbol("mae_" * b) => mean_(abs.(dias[!, b * "_nuestro"] .- dias[!, b * "_ggir"])) for (b, _) in BANDAS)...,
           dif_despierto = mean_(dias.despierto_nuestro .- dias.despierto_ggir),
           parejas_sri = count(okk), SRI_nuestro = mean_(js.SRI[okk]), SRI_ggir = mean_(sri_g[okk]),
           mae_SRI = mean_(abs.(js.SRI[okk] .- sri_g[okk])),
           noches_eficiencia = count(oke), efic_nuestra = mean_(je.eficiencia[oke]), efic_ggir = mean_(ef_g[oke]))
    return res, dias
end

function tabla_tex(df)
    io = IOBuffer()
    println(io, "% Generado por scripts/96_validar_intensidad_sri_ggir.jl — tabla individual, no versionar.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}X *{$(nrow(df))}{r}@{}}")
    println(io, "\\toprule")
    println(io, " & ", join(replace.(df.id, "_" => "\\_"), " & "), " \\\\")
    println(io, "\\midrule")
    f(n, g) = println(io, n, " & ", join([g(r) for r in eachrow(df)], " & "), " \\\\")
    sec(t) = println(io, "\\multicolumn{$(nrow(df) + 1)}{@{}l}{\\textit{$t}} \\\\")
    sec("Minutos al día por intensidad, despierto (GGIR parte 5)")
    f("Días comparados", r -> string(r.dias_intensidad))
    f("Tiempo despierto: nuestro \$-\$ GGIR (min)", r -> fmt(r.dif_despierto))
    for (b, _) in BANDAS
        f(uppercasefirst(b) * ": diferencia media / absoluta media (min)", r -> fmt(r[Symbol("dif_" * b)]) * " / " * fmt(r[Symbol("mae_" * b)]))
    end
    sec("Regularidad y eficiencia del sueño (GGIR parte 4)")
    f("SRI: \\Brincadiano{} / GGIR (parejas de días)", r -> fmt(r.SRI_nuestro) * " / " * fmt(r.SRI_ggir) * " (" * string(r.parejas_sri) * ")")
    f("SRI: diferencia absoluta media", r -> fmt(r.mae_SRI))
    f("Eficiencia: \\Brincadiano{} / GGIR (noches)", r -> fmt(r.efic_nuestra; d = 2) * " / " * fmt(r.efic_ggir; d = 2) * " (" * string(r.noches_eficiencia) * ")")
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
    return String(take!(io))
end

if !("--confirmo" in ARGS)
    println("Contrastes R2 y R3 frente a GGIR (partes 1–5). Lanza R: añade --confirmo.")
    exit(0)
end
if !("--sin-ggir" in ARGS && isfile(joinpath(GG, "ggir_part5_dias.csv")))
    println("GGIR (partes 1–5) …"); correr_ggir()
else
    println("GGIR: se reutiliza la salida extraída en outputs/validation/ggir_parte5")
end
p4 = CSV.read(joinpath(GG, "ggir_part4_noches.csv"), DataFrame; missingstring = "NA")
p5 = CSV.read(joinpath(GG, "ggir_part5_dias.csv"), DataFrame; missingstring = "NA")
ids = sort(filter(id -> isfile(joinpath(INTERIM, id, "angulos_5s.arrow")) && id in p5.id, readdir(INTERIM)))
out = [comparar(id, p4, p5) for id in ids]
res = DataFrame(first.(out)); dias = vcat(last.(out)...)
mkpath(OUT)
CSV.write(joinpath(OUT, "ggir_intensidad_sri.csv"), res)
CSV.write(joinpath(OUT, "ggir_intensidad_dias.csv"), dias)
write(joinpath(OUT, "ggir_intensidad_sri.tex"), tabla_tex(res))
for x in eachrow(res)
    @printf("%-12s %d días: Δ despierto %+.0f min; Δ sed %+.1f, lig %+.1f, mod %+.1f, vig %+.2f min (|Δ| %.1f/%.1f/%.1f/%.2f); SRI %.1f vs %.1f (%d parejas, |Δ| %.1f); eficiencia %.2f vs %.2f (%d noches)\n",
            x.id, x.dias_intensidad, x.dif_despierto, x.dif_sedentario, x.dif_ligera, x.dif_moderada, x.dif_vigorosa,
            x.mae_sedentario, x.mae_ligera, x.mae_moderada, x.mae_vigorosa, x.SRI_nuestro, x.SRI_ggir, x.parejas_sri, x.mae_SRI,
            x.efic_nuestra, x.efic_ggir, x.noches_eficiencia)
end
println("Tablas en ", relpath(OUT, ROOT))
