# 21 · Reposo nocturno con pocas noches de uso.
#
# Qué hace
#   Responde a la pregunta del grupo (15-oct-2026, PRO-02): si muchos niños llevan el dispositivo
#   pocas noches, ¿se puede estimar su reposo nocturno con pocas? Para cada registro, toma las
#   noches con ventana del script 08 y mide, para cada número k de noches, cuánto se aleja la
#   media de k noches de la de todas (inicio, fin, punto medio y duración del periodo). Se
#   enumeran todos los subconjuntos de k noches (`Sleep.night_subsample`).
#   No es la precisión entre niños: con un registro por tipo (uno con el dispositivo puesto, otro
#   con retirada) solo describe la variación entre noches de un mismo niño.
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/21_sueno_pocas_noches.jl [ID ...]   # sin ID, los de sueno_resumen.csv
#
# Depende de
#   08 (outputs/sleep/tables/sueno_<ID>_noches.csv).
#
# Salida (individual; no se versiona; las noches van numeradas, sin fecha de calendario)
#   outputs/sleep/tables/sens_noches.csv              por registro, k y variable
#   outputs/sleep/figures/dos_registros_noches.{pdf,png}   con su ficha de procedencia
#
# Consola: solo el ID y cifras resumen.

using Brincadiano.Sleep, Brincadiano.Viz
using CairoMakie, CSV, DataFrames, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SUENO = joinpath(ROOT, "outputs", "sleep", "tables")
const OUT_T = SUENO
const OUT_F = joinpath(ROOT, "outputs", "sleep", "figures")
mkpath(OUT_T); mkpath(OUT_F)

ids = isempty(ARGS) ? String.(CSV.read(joinpath(SUENO, "sueno_resumen.csv"), DataFrame).id) : String.(ARGS)

todo = DataFrame(); datos = Tuple{String,DataFrame}[]
for id in ids
    f = joinpath(SUENO, "sueno_$(id)_noches.csv")
    if !isfile(f)
        @printf("omitido %s: no hay %s (¿se ejecutó 08?)\n", id, basename(f))
        continue
    end
    n = CSV.read(f, DataFrame)
    n.metodo = String.(n.metodo)
    ok = n[n.metodo .!= "sin ventana", :]
    if nrow(ok) < 2
        @printf("omitido %s: menos de 2 noches con ventana\n", id)
        continue
    end
    s = night_subsample(n)
    insertcols!(s, 1, :id => fill(id, nrow(s)))
    append!(todo, s)
    push!(datos, (id, s))
    nn = nrow(ok)
    m = s[s.variable .== "punto_medio_h", :]; d = s[s.variable .== "duracion_h", :]
    @printf("%s: %d noches (%d por postura, %d con retirada)\n", id, nn,
            count(==("postura"), ok.metodo), count(==("retirada"), ok.metodo))
    for k in (1, 2, 3, nn - 1)
        k in m.k || continue
        a = m[m.k .== k, :][1, :]; b = d[d.k .== k, :][1, :]
        @printf("  k = %d: punto medio, mediana %.2f h y peor %.2f h; duración, mediana %.2f h y peor %.2f h (%d subconjuntos)\n",
                k, a.mediana_abs, a.max_abs, b.mediana_abs, b.max_abs, a.subconjuntos)
    end
end

isempty(todo) && error("no hay ningún registro con al menos 2 noches")
CSV.write(joinpath(OUT_T, "sens_noches.csv"), todo)
guardar(joinpath(OUT_F, "dos_registros_noches"), plot_nights_subsample(datos);
        procedencia = (; tipo = :datos, registros = [d[1] for d in datos], funcion = "plot_nights_subsample",
                       entradas = ["outputs/sleep/tables/sueno_<ID>_noches.csv"]))
println("escrito: sens_noches.csv y dos_registros_noches")
