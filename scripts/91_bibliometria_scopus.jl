# 91 · Mapa bibliométrico de los colaboradores a partir de exportaciones de Scopus.
#
# Qué hace
#   Cuenta, para dos perfiles de autor exportados de Scopus, artículos, citas, periodo
#   y el reparto por temas (obesidad, HOMA-IR, actividad física, acelerometría, sueño,
#   ritmo circadiano, luz, GGIR/ENMO, etc.). No toca acelerometría ni data/. Los temas
#   se buscan con expresiones regulares, sin distinguir mayúsculas, sobre título,
#   resumen y palabras clave de autor e indexadas. Es un cribado: un artículo cuenta
#   en un tema si menciona el término, aunque no sea su objeto. «Luz» no busca la
#   palabra suelta «light», que en estos textos es casi siempre «light physical
#   activity».
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/91_bibliometria_scopus.jl <Gil-Campos.csv> <Jurado-Castro.csv>
#   Los dos argumentos son obligatorios y van en ese orden.
#
# Depende de
#   De ningún script del pipeline. Los CSV se pasan desde fuera del repositorio
#   (exportación de Scopus del 2026-10-02, todos los campos). No se versionan:
#   la licencia de Scopus no lo permite. Son metadatos bibliográficos públicos,
#   no datos de participantes.
#
# Lo usan
#   Solo el capítulo 1. Las salidas ya quedan en assets/.../report/ y se versionan;
#   no pasan por `make sync-assets`. Hay que recompilar el TikZ de publicaciones
#   anuales si se regenera el .tikz (`make tikz-preview` en report/).
#
# Entrada
#   Dos CSV de Scopus, rutas absolutas o relativas al directorio de trabajo.
#
# Salida versionada (agregados; material editorial del capítulo 1)
#   report/assets/01_contexto_clinico/report/tables/bibliometria-perfil.tex
#   report/assets/01_contexto_clinico/report/tables/bibliometria-temas.tex
#   report/assets/01_contexto_clinico/report/tikz/publicaciones-anuales.tikz
#
# Consola: recuentos por autor y por tema. Ningún identificador de participante.

using CSV, DataFrames, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ASSETS = joinpath(ROOT, "report", "assets", "01_contexto_clinico", "report")

const TEMAS = [
    ("Obesidad o sobrepeso",                 "obes-, overweight",                          r"obes|overweight"i),
    ("HOMA-IR o resistencia a la insulina",  "HOMA, insulin resist-",                      r"homa|insulin resist"i),
    ("Actividad física o sedentarismo",      "physical activity, sedentar-, exercise, fitness", r"physical activity|sedentar|exercise|fitness"i),
    ("Metabolómica o microbiota",            "metabolom-, microbio-",                      r"metabolom|microbio"i),
    ("Cohorte COHORTE_A",                       "COHORTE_A",                                     r"cohorte_a"i),
    ("Acelerometría",                        "acceleromet-, ActiGraph",                    r"acceleromet|actigraph"i),
    ("Sueño",                                "sleep",                                      r"\bsleep"i),
    ("Cohorte GENOBOX",                      "GENOBOX",                                    r"genobox"i),
    ("Aprendizaje automático o IA",          "machine/deep learning, artificial intelligence, random forest",
                                             r"machine learning|artificial intelligence|deep learning|random forest"i),
    ("Ritmo circadiano o cronotipo",         "circadian, chrono-, rest-activity",          r"circadian|chrono|rest-activity"i),
    ("Exposición a la luz",                  "light exposure, lux, illuminance, light-dark", r"\blight exposure|\blux\b|illuminance|light-dark"i),
    ("Señal RAW, GGIR o ENMO",               "GGIR, ENMO, raw acceleration/data",          r"\bggir\b|\benmo\b|raw acceleration|raw data"i),
]

leer(path) = CSV.read(path, DataFrame; stringtype = String)

texto(r) = join((coalesce(r[c], "") for c in ("Title", "Abstract", "Author Keywords", "Index Keywords")), ' ')

citas(df) = [ismissing(c) ? 0 : Int(c) for c in df[!, "Cited by"]]

function indice_h(c)
    s = sort(c; rev = true)
    return count(i -> s[i] >= i, eachindex(s))
end

acceso_abierto(df) = count(x -> !ismissing(x) && !isempty(x), df[!, "Open Access"])

function revistas(df, n)
    cuenta = combine(groupby(df, "Source title"), nrow => :n)
    sort!(cuenta, [order(:n; rev = true), order("Source title")])
    return [(r["Source title"], r.n) for r in eachrow(first(cuenta, n))]
end

tex(s) = replace(s, "&" => "\\&", "%" => "\\%", "_" => "\\_")

function tabla_perfil(io, M, J, comunes)
    cm, cj = citas(M), citas(J)
    años(df) = extrema(df.Year)
    fila(nombre, a, b) = println(io, nombre, " & ", a, " & ", b, " \\\\")
    println(io, "% Generado por scripts/91_bibliometria_scopus.jl — no editar a mano.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}p{0.24\\linewidth} X X@{}}")
    println(io, "\\toprule")
    println(io, " & Dra.~M.~Gil Campos & Dr.~J.~M.~Jurado Castro \\\\")
    println(io, "\\midrule")
    fila("Registros en Scopus", nrow(M), nrow(J))
    fila("Periodo", join(años(M), "--"), join(años(J), "--"))
    fila("Artículos originales / revisiones",
         "$(count(==("Article"), M[!, "Document Type"])) / $(count(==("Review"), M[!, "Document Type"]))",
         "$(count(==("Article"), J[!, "Document Type"])) / $(count(==("Review"), J[!, "Document Type"]))")
    fila("Citas totales", sum(cm), sum(cj))
    fila("Índice \$h\$ (sobre esta exportación)", indice_h(cm), indice_h(cj))
    fila("En acceso abierto", "$(acceso_abierto(M)) ($(round(Int, 100acceso_abierto(M)/nrow(M)))\\,\\%)",
                              "$(acceso_abierto(J)) ($(round(Int, 100acceso_abierto(J)/nrow(J)))\\,\\%)")
    fila("Publicaciones conjuntas", "$comunes ($(round(Int, 100comunes/nrow(M)))\\,\\%)",
                                    "$comunes ($(round(Int, 100comunes/nrow(J)))\\,\\%)")
    fila("Revistas más frecuentes",
         join(("$(tex(t)) ($n)" for (t, n) in revistas(M, 4)), "; "),
         join(("$(tex(t)) ($n)" for (t, n) in revistas(J, 4)), "; "))
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
end

function tabla_temas(io, M, J, eids_comunes)
    tm, tj = texto.(eachrow(M)), texto.(eachrow(J))
    comun = M.EID .∈ Ref(eids_comunes)
    println(io, "% Generado por scripts/91_bibliometria_scopus.jl — no editar a mano.")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}p{0.27\\linewidth} r r r >{\\raggedright\\arraybackslash}X@{}}")
    println(io, "\\toprule")
    println(io, "Tema & Gil Campos & Jurado Castro & Conjuntas & Términos buscados \\\\")
    println(io, "\\midrule")
    for (nombre, terminos, re) in TEMAS
        m = count(t -> occursin(re, t), tm)
        j = count(t -> occursin(re, t), tj)
        c = count(i -> comun[i] && occursin(re, tm[i]), eachindex(tm))
        resaltar = m + j == 0 || (nombre == "Ritmo circadiano o cronotipo" && c == 0)
        celda(x) = resaltar ? "\\textbf{$x}" : string(x)
        println(io, nombre, " & ", celda(m), " & ", celda(j), " & ", celda(c),
                " & {\\footnotesize\\emph{", terminos, "}} \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabularx}")
end

function figura_anual(io, M, J, eids_comunes)
    años = minimum(vcat(M.Year, J.Year)):maximum(vcat(M.Year, J.Year))
    n(df, a) = count(==(a), df.Year)
    nc(a) = count(i -> M.Year[i] == a && M.EID[i] in eids_comunes, 1:nrow(M))
    coords(f) = join(("($a,$(f(a)))" for a in años), " ")
    println(io, "% Generado por scripts/91_bibliometria_scopus.jl — no editar a mano.")
    println(io, "% Publicaciones por año en Scopus (exportación del 2026-10-02; 2026 incompleto).")
    println(io, """
    \\begin{tikzpicture}
    \\begin{axis}[
      width=\\linewidth, height=5.2cm,
      ybar=0pt, bar width=3.2pt,
      xmin=$(first(años) - 0.7), xmax=$(last(años) + 0.7),
      xtick={$(join(unique(vcat(collect(first(años):5:last(años)), last(años))), ","))},
      xticklabel style={/pgf/number format/1000 sep={}},
      ymin=0, ylabel={Publicaciones},
      axis lines*=left, ymajorgrids, grid style={black!10},
      legend style={draw=none, at={(0.02,0.98)}, anchor=north west, font=\\sffamily\\footnotesize},
      legend cell align=left,
      tick label style={font=\\sffamily\\footnotesize}, label style={font=\\sffamily\\footnotesize},
    ]
    \\addplot[fill=brinBlue, draw=none] coordinates {$(coords(a -> n(M, a)))};
    \\addplot[fill=brinTeal, draw=none] coordinates {$(coords(a -> n(J, a)))};
    \\addplot[fill=brinAmber, draw=none] coordinates {$(coords(nc))};
    \\legend{Gil Campos, Jurado Castro, Conjuntas}
    \\end{axis}
    \\end{tikzpicture}""")
end

function main(args)
    length(args) == 2 || error("Uso: scripts/91_bibliometria_scopus.jl <Gil-Campos.csv> <Jurado-Castro.csv>")
    M, J = leer(args[1]), leer(args[2])
    eids_comunes = Set(M.EID) ∩ Set(J.EID)
    mkpath(joinpath(ASSETS, "tables")); mkpath(joinpath(ASSETS, "tikz"))
    open(io -> tabla_perfil(io, M, J, length(eids_comunes)), joinpath(ASSETS, "tables", "bibliometria-perfil.tex"), "w")
    open(io -> tabla_temas(io, M, J, eids_comunes), joinpath(ASSETS, "tables", "bibliometria-temas.tex"), "w")
    open(io -> figura_anual(io, M, J, eids_comunes), joinpath(ASSETS, "tikz", "publicaciones-anuales.tikz"), "w")
    @printf("Gil Campos: %d registros · Jurado Castro: %d · conjuntos: %d\n", nrow(M), nrow(J), length(eids_comunes))
    println("Tablas y figura en ", relpath(ASSETS, ROOT))
end

main(ARGS)
