# 22 · Siestas candidatas (inactividad sostenida diurna) en niños de 3 años.
#
# Qué hace
#   Aplica `Sleep.nap_candidates` a los registros ya procesados: bloques de inactividad
#   sostenida (SIB de van Hees 2015) fuera del reposo nocturno estimado, con el dispositivo
#   puesto, que empiezan entre las 9 y las 18 h y duran de 15 a 240 min (los mismos parámetros
#   de ejemplo que la ayuda de GGIR 3.3.9 da para `possible_nap_window` y `possible_nap_dur`).
#   Resume por registro cuántas hay, cuántos días tienen alguna y cuánto duran, y repite con
#   otras duraciones mínimas y con unión de bloques para ver cuánto cambia el recuento.
#   No está validada: no hay siestas declaradas ni etiquetadas para estos registros (ficha
#   SUE-02) y no se ha contrastado con GGIR.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/22_siestas.jl [ID ...]   # sin ID, quien tenga angulos_5s.arrow
#
# Depende de
#   08 (data/interim/<ID>/angulos_5s.arrow, con el idle sleep heredado, y el no uso).
#
# Salida (individual; no se versiona; los días van numerados, sin fecha de calendario)
#   outputs/sleep/tables/siestas_<ID>.csv      una fila por candidata
#   outputs/sleep/tables/siestas_resumen.csv   por registro y variante: candidatas y días con alguna
#
# Consola: solo el ID y cifras resumen.

using Brincadiano.Sleep
using CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIM = joinpath(ROOT, "data", "interim")
const OUT_T = joinpath(ROOT, "outputs", "sleep", "tables")
mkpath(OUT_T)

const VARIANTES = [("referencia (15–240 min, sin unir)", (; dur_min = (15, 240), gap_min = 0)),
                   ("≥ 30 min", (; dur_min = (30, 240), gap_min = 0)),
                   ("≥ 60 min", (; dur_min = (60, 240), gap_min = 0)),
                   ("uniendo bloques a menos de 10 min", (; dur_min = (15, 240), gap_min = 10))]

ids = isempty(ARGS) ? sort([d for d in readdir(INTERIM) if isfile(joinpath(INTERIM, d, "angulos_5s.arrow"))]) : String.(ARGS)
resumen = DataFrame()
for id in ids
    dir = joinpath(INTERIM, id)
    isfile(joinpath(dir, "angulos_5s.arrow")) || (@printf("omitido %s: no hay angulos_5s.arrow (¿se ejecutó 08?)\n", id); continue)
    ep = load_angles(dir)
    eje = longitudinal_axis(ep; mask = .!ep.nonwear)
    noches = sleep_nights(ep; axis = eje, nonwear = ep.nonwear)
    if nrow(noches) == 0
        @printf("omitido %s: sin noches completas (un registro corto no tiene reposo nocturno)\n", id)
        continue
    end
    ndias = length(unique(Date.(ep.start[ep.n .> 0])))
    for (nombre, kw) in VARIANTES
        s = nap_candidates(ep, noches; kw...)
        push!(resumen, (id = id, variante = nombre, noches_con_ventana = count(!=("sin ventana"), noches.metodo),
                        candidatas = nrow(s), dias_con_alguna = length(unique(s.dia)), dias_registro = ndias,
                        duracion_mediana_min = isempty(s) ? NaN : median(s.duracion_min),
                        duracion_max_min = isempty(s) ? NaN : maximum(s.duracion_min),
                        pct_idle_mediana = isempty(s) ? NaN : median(s.pct_idle),
                        junto_a_no_uso = count(s.junto_a_no_uso)); promote = true)
        if nombre == VARIANTES[1][1]
            CSV.write(joinpath(OUT_T, "siestas_$(id).csv"), s)
            @printf("%s: %d candidatas en %d de %d días; mediana %.0f min (máx. %.0f); %d junto a no uso; idle sleep mediano %.0f %%\n",
                    id, nrow(s), length(unique(s.dia)), ndias, isempty(s) ? NaN : median(s.duracion_min),
                    isempty(s) ? NaN : maximum(s.duracion_min), count(s.junto_a_no_uso), isempty(s) ? NaN : median(s.pct_idle))
        end
    end
end
isempty(resumen) && error("no se procesó ningún registro")
CSV.write(joinpath(OUT_T, "siestas_resumen.csv"), resumen)
println("escrito: siestas_resumen.csv")
