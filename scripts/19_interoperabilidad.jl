# 19 · Lista de comprobación de interoperabilidad GT3X+ ↔ wGT3X-BT (anexo B, tab:interoperabilidad).
#
# Qué hace
#   Para cada registro de data/interim/ rellena las filas de «Lectura» que se pueden comprobar sin
#   referencia externa y mide qué pasaría si el mismo movimiento se hubiera registrado con un GT3X+
#   (rango ±6 g) en vez de un wGT3X-BT (±8 g): cuenta las muestras con algún eje por encima de ±6 g
#   y recalcula la ENMO por epochs de 5 s con la aceleración recortada a ±6 g por eje. Los
#   registros de GT3X+ propiamente dichos (PRO-06) llegarán con la cohorte: el mismo script se
#   lanza sobre ellos y la tabla por modelo (columna `modelo`) rellena la otra columna.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/19_interoperabilidad.jl            # todos los de data/interim
#   julia --project=scripts scripts/19_interoperabilidad.jl 6021_V00   # solo esos IDs
#
# Depende de
#   01 (o el lote 16). No usa data/raw/ ni herramientas externas.
#
# Salida (no se versiona)
#   outputs/interop/tables/interop_ids.csv       una fila por registro (solo el ID)
#   outputs/interop/tables/interop_modelos.csv   agregado por modelo y firmware
#
# Límites: ENMO sin autocalibrar en las dos versiones (se compara lo mismo con lo mismo);
# recortar por eje es una aproximación a la saturación del sensor de ±6 g, no su simulación.
# Consola: solo el ID y recuentos.

using Brincadiano, Brincadiano.GT3X, Brincadiano.Metrics
using CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT = joinpath(ROOT, "outputs", "interop", "tables")
const LIM = 6.0
const UMBRAL_VIGOROSO_MG = 464.6      # Hildebrand, cadera (ACT-01)

pedidos = filter(a -> !startswith(a, "--"), ARGS)
ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "diagnostics.arrow")), readdir(INTERIM)))
isempty(pedidos) || (ids = filter(in(pedidos), ids))
isempty(ids) && error("No hay registros leídos en data/interim/")

filas = NamedTuple[]
for id in ids
    r = read_interim(joinpath(INTERIM, id))
    fs = round(Int, r.info.sample_rate)
    acc = r.acc
    n = nrow(acc)
    ax = (acc.x, acc.y, acc.z)
    maxabs = maximum(maximum(abs, a) for a in ax)
    sobre = [count(>(LIM), abs.(a)) for a in ax]
    alguno = count(i -> abs(acc.x[i]) > LIM || abs(acc.y[i]) > LIM || abs(acc.z[i]) > LIM, 1:n)
    rec = DataFrame(time = acc.time, x = clamp.(acc.x, -LIM, LIM), y = clamp.(acc.y, -LIM, LIM),
                    z = clamp.(acc.z, -LIM, LIM))
    e0 = epoch_metrics(acc, r.gaps, fs, Second(5))
    e1 = epoch_metrics(rec, r.gaps, fs, Second(5))
    v = .!isnan.(e0.enmo)
    vig0 = count(>=(UMBRAL_VIGOROSO_MG), e0.enmo[v]); vig1 = count(>=(UMBRAL_VIGOROSO_MG), e1.enmo[v])
    ep_cambian = count(abs.(e0.enmo[v] .- e1.enmo[v]) .> 1e-9)
    push!(filas, (; id, formato = String(r.format), modelo = r.info.device_type, firmware = r.info.firmware,
        fs_nominal = fs,         escala_lsb_g = r.scale, rango_declarado_g = r.info.accel_max,
        max_abs_g = round(maxabs; digits = 3), muestras = n,
        sobre6g_x = sobre[1], sobre6g_y = sobre[2], sobre6g_z = sobre[3], sobre6g_alguno = alguno,
        pct_sobre6g = round(100 * alguno / n; digits = 5),
        epochs5s = count(v), epochs_con_cambio = ep_cambian,
        enmo_medio_mg = round(mean(e0.enmo[v]); digits = 3), enmo_medio_recortada_mg = round(mean(e1.enmo[v]); digits = 3),
        vigorosos_5s = vig0, vigorosos_5s_recortada = vig1,
        idle_pct = round(100 * (1 - n / Brincadiano.GT3X.expected_samples(r)); digits = 1),
        luz_registros = nrow(r.lux), capsense_registros = nrow(r.capsense)))
    @printf("%-14s %-8s máx %.2f g   muestras > %.0f g: %d   epochs que cambian: %d/%d   vigorosos %d → %d\n",
            id, r.info.device_type, maxabs, LIM, alguno, ep_cambian, count(v), vig0, vig1)
end
df = DataFrame(filas)
mkpath(OUT)
CSV.write(joinpath(OUT, "interop_ids.csv"), df)
agg = combine(groupby(df, [:modelo, :firmware, :formato]), :id => length => :registros,
    :muestras => sum => :muestras, :sobre6g_alguno => sum => :sobre6g_alguno,
    :vigorosos_5s => sum => :vigorosos_5s, :vigorosos_5s_recortada => sum => :vigorosos_5s_recortada)
CSV.write(joinpath(OUT, "interop_modelos.csv"), agg)
println("Salida en ", relpath(OUT, ROOT))
