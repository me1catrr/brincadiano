# 26 · Brincadiano frente a la parte 5 de GGIR enviada por el grupo (COHORTE_A V03, mismos .gt3x).
#
# Qué hace
#   Lee del Excel de la parte 5 (`part5_personsummary_MM_L63.3M142.6V464.6`, GGIR 3.2.6 con
#   config_cohorte_a.csv) solo las columnas de la lista blanca y las compara, registro a registro,
#   con lo que calcula el pipeline sobre los mismos .gt3x. Solo se comparan magnitudes con la
#   misma definición:
#     · Días válidos: GGIR exige 16 h (`includedaycrit = 16`); aquí, días naturales con ≥ 16 h de
#       uso según 04 (mismo método de no uso, «2023»). PROVISIONAL: la parte 5 añade además que
#       2/3 del tiempo despierto sea válido.
#     · Minutos/día de actividad moderada y vigorosa (MVPA = moderada + vigorosa) en epochs de 5 s
#       con los mismos umbrales (142,6 y 464,6 mg), media simple de los días válidos (`_pla`). Se
#       compara el día entero: durante el sueño casi no hay MVPA. Diferencias esperadas: GGIR
#       imputa el no uso (`do.imp`) y los huecos del idle sleep (`imputeTimegaps`); aquí no.
#     · Periodo de sueño (`dur_spt_min_pla`, inicio y despertar): solo DESCRIPTIVO. GGIR lo
#       calculó con HDCZA y `sensor.location = wrist`, pero el dispositivo iba en la cadera
#       derecha (PRO-04); el pipeline usa HorAngle de cadera. No vale como validación.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/26_comparar_ggir_cohorte_a.jl     # tras 04, 05 y 08
#
# Salida (no se versiona)
#   outputs/agd/tables/ggir_cohorte_a_comparacion.csv   una fila por registro
#
# Consola: solo el ID y diferencias.

using Brincadiano, Brincadiano.GT3X, Brincadiano.Metrics
using Brincadiano.Rutas: Rutas, raw_files
using CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = Rutas.interim_dir()
const OUT = joinpath(ROOT, "outputs", "agd", "tables")
const UMBRAL_DIA_H = 16
const COLUMNAS = ["filename", "Nvaliddays", "dur_day_total_MOD_min_pla", "dur_day_total_VIG_min_pla",
                  "nonwear_perc_day_pla", "dur_spt_min_pla", "sleeponset_pla", "wakeup_pla"]

# Lector mínimo de .xlsx (una hoja): celdas por referencia, cadenas compartidas o en línea.
function leer_xlsx(path; columnas)
    z = GT3X.ZipReader(read(path))
    ss = String[]
    if "xl/sharedStrings.xml" in GT3X.zip_names(z)
        x = GT3X.zip_readentry(z, "xl/sharedStrings.xml", String)
        ss = [replace(m.captures[1], r"<[^>]+>" => "") for m in eachmatch(r"<si>(.*?)</si>"s, x)]
    end
    x = GT3X.zip_readentry(z, "xl/worksheets/sheet1.xml", String)
    celdas = Dict{Tuple{Int,String},String}()
    for m in eachmatch(r"<c r=\"([A-Z]+)(\d+)\"([^>]*?)(?:/>|>(.*?)</c>)"s, x)
        col, fila, attr, cuerpo = m.captures
        cuerpo === nothing && continue
        v = match(r"<v>(.*?)</v>"s, cuerpo)
        val = v === nothing ? (t = match(r"<t[^>]*>(.*?)</t>"s, cuerpo); t === nothing ? "" : t.captures[1]) :
              occursin("t=\"s\"", attr) ? ss[parse(Int, v.captures[1]) + 1] : v.captures[1]
        celdas[(parse(Int, fila), col)] = val
    end
    cab = Dict(v => c for ((f, c), v) in celdas if f == 1)
    faltan = setdiff(columnas, keys(cab)); isempty(faltan) || error("Faltan columnas en el Excel: $faltan")
    nf = maximum(first, keys(celdas))
    return DataFrame([c => [get(celdas, (f, cab[c]), "") for f in 2:nf] for c in columnas])
end

num(s) = (v = tryparse(Float64, s); v === nothing ? missing : v)
# GGIR da las horas de sueño en horas desde la medianoche del día (pueden pasar de 24).
hora(s) = (v = num(s); ismissing(v) ? missing : round(v; digits = 2))

xs = filter(f -> occursin(r"(?i)ggir.*cohorte_a", basename(f)), raw_files("xlsx"; proyecto = "COHORTE_A"))
length(xs) == 1 || error("Se esperaba un Excel de la parte 5 de GGIR de COHORTE_A en data/raw/ y hay $(length(xs))")
xls = only(xs)
g = leer_xlsx(xls; columnas = COLUMNAS)
g.id = replace.(g.filename, r"\.gt3x$" => "")

sueno = let f = joinpath(ROOT, "outputs", "sleep", "tables", "sueno_resumen.csv")
    isfile(f) ? CSV.read(f, DataFrame) : nothing
end

filas = NamedTuple[]
for r in eachrow(g)
    id = r.id
    fdias = joinpath(ROOT, "outputs", "qc", "tables", "qc_$(id)_dias.csv")
    fep = joinpath(INTERIM, id, "epochs_5s.arrow")
    (isfile(fdias) && isfile(fep)) || (println("  ", id, ": falta 04 o 05; se omite"); continue)
    d = CSV.read(fdias, DataFrame)
    validos = Set(Date.(string.(d.dia[d.uso_h .>= UMBRAL_DIA_H])))
    ep = DataFrame(GT3X.Arrow.Table(fep); copycols = true)
    it = intensity_time(ep, trues(nrow(ep)))
    it = it[in.(it.dia, Ref(validos)), :]
    mod_b = isempty(it) ? missing : mean(it.moderada)
    vig_b = isempty(it) ? missing : mean(it.vigorosa)
    mod_g = num(r.dur_day_total_MOD_min_pla); vig_g = num(r.dur_day_total_VIG_min_pla)
    s = sueno === nothing ? nothing : (k = findfirst(==(id), sueno.id); k === nothing ? nothing : sueno[k, :])
    push!(filas, (; id, dias_validos_ggir = num(r.Nvaliddays), dias_validos_16h = length(validos),
        mod_min_ggir = mod_g, mod_min_brinc = mod_b, vig_min_ggir = vig_g, vig_min_brinc = vig_b,
        mvpa_min_ggir = mod_g + vig_g, mvpa_min_brinc = mod_b + vig_b,
        d_mvpa_min = (d = mod_b + vig_b - mod_g - vig_g; ismissing(d) ? missing : round(d; digits = 1)),
        no_uso_pct_ggir = num(r.nonwear_perc_day_pla),
        spt_min_ggir_hdcza_muneca = num(r.dur_spt_min_pla), inicio_ggir_h = hora(r.sleeponset_pla),
        despertar_ggir_h = hora(r.wakeup_pla),
        spt_min_brinc_horangle = s === nothing ? missing : round(60 * s.duracion_h; digits = 1),
        inicio_brinc_h = s === nothing ? missing : round(s.inicio_h; digits = 2),
        despertar_brinc_h = s === nothing ? missing : round(s.fin_h; digits = 2),
        noches_postura_brinc = s === nothing ? missing : s.postura, spt_validacion = "NO VÁLIDO (HDCZA de muñeca sobre cadera)"))
end
df = DataFrame(filas)
mkpath(OUT)
CSV.write(joinpath(OUT, "ggir_cohorte_a_comparacion.csv"), df)
@printf("%-11s %9s %9s %9s %9s %8s\n", "ID", "días GGIR", "días 16 h", "MVPA GGIR", "MVPA B.", "Δ min")
for r in eachrow(df)
    @printf("%-11s %9s %9d %9.1f %9.1f %8.1f\n", r.id, string(r.dias_validos_ggir), r.dias_validos_16h,
            coalesce(r.mvpa_min_ggir, NaN), coalesce(r.mvpa_min_brinc, NaN), coalesce(r.d_mvpa_min, NaN))
end
println("Salida en ", relpath(OUT, ROOT))
