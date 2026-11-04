"""
    Rutas

Configuración central de rutas y de la estructura de `data/raw/` (`docs/datos_raw.md`).

- `data_dir()`: `ENV["BRINCADIANO_DATA"]` si existe; si no, `<repositorio>/data`. Así los datos
  pueden vivir en otro volumen (p. ej. un disco cifrado) sin tocar el código.
- `raw_dir()`, `interim_dir()`, `outputs_dir()`.
- `raw_files(formato; proyecto)`: archivos originales de un formato, en cualquier subcarpeta de
  `raw_dir()` (la estructura por proyecto y, durante la transición, la raíz).
- `find_raw(nombre)`: ruta de un archivo original por su nombre; `nothing` si no está.
- `clasificar(nombre)`: proyecto, estado de la atribución, sujeto, visita, formato y carpeta de
  destino, deducidos **solo del nombre** y de lo declarado por el grupo. Lo que el nombre no
  permite atribuir va a `pendientes/`.

No lee el contenido de ningún archivo.
"""
module Rutas

export data_dir, raw_dir, interim_dir, outputs_dir, raw_files, find_raw, clasificar, id_registro, PROYECTOS

const REPO = normpath(joinpath(@__DIR__, ".."))

data_dir() = get(ENV, "BRINCADIANO_DATA", joinpath(REPO, "data"))
raw_dir() = joinpath(data_dir(), "raw")
interim_dir() = joinpath(data_dir(), "interim")
outputs_dir() = get(ENV, "BRINCADIANO_OUTPUTS", joinpath(REPO, "outputs"))

"Proyecto → carpeta bajo `data/raw/`."
const PROYECTOS = (COHORTE_A = "cohorte_a", MELIPOP = "melipop", MiopiaUNAV = "miopia-unav",
                   pruebas = "pruebas", pendiente = "pendientes")

"Formato → subcarpeta dentro del proyecto."
const SUBCARPETAS = Dict("gt3x" => "gt3x", "agd" => "agd", "actilife_raw_csv" => "actilife_raw_csv",
                         "ggir_config" => "config", "xlsx" => "exportaciones", "otro" => "otros")

"""
    formato(nombre) -> String

`gt3x`, `agd`, `actilife_raw_csv`, `ggir_config`, `xlsx` u `otro`, por la extensión y el nombre.
"""
function formato(nombre::AbstractString)
    l = lowercase(nombre)
    endswith(l, ".gt3x") && return "gt3x"
    endswith(l, ".agd") && return "agd"
    endswith(l, "raw.csv") && return "actilife_raw_csv"
    startswith(l, "config") && endswith(l, ".csv") && return "ggir_config"
    endswith(l, ".xlsx") && return "xlsx"
    return "otro"
end

"""
    id_registro(nombre) -> String

ID pseudonimizado del registro: el nombre sin extensión, sin el sufijo `RAW` de las
exportaciones CSV de ActiLife ni la duración del *epoch* que ActiLife añade a los AGD
(`H-01-0011sec.agd` → `H-01-001`; `C6005_V00.agd` → `C6005_V00`).
"""
id_registro(nombre::AbstractString) =
    replace(splitext(basename(nombre))[1], r"RAW$" => "", r"(1|5|10|15|30|60)sec$" => "")

"""
    clasificar(nombre) -> NamedTuple

Atribución de un archivo original por su nombre: `proyecto` (clave de [`PROYECTOS`](@ref)),
`estado` («declarado (…)» con su fuente, o «hipótesis»), `sujeto`, `visita`, `formato`,
`colocacion`, `procedencia` y `destino` (ruta relativa a `raw_dir()`).

Los `C60xx_V00` no se atribuyen: el prefijo `C` sugiere COHORTE_A, pero el grupo advirtió que los
códigos coinciden por azar entre proyectos; van a `pendientes/` hasta que haya evidencia.
"""
function clasificar(nombre::AbstractString)
    fmt = formato(nombre)
    base = id_registro(nombre)
    p, estado, sujeto, visita, coloc, proc = "pendiente", "sin atribuir", "", "", "", ""
    if occursin(r"^\d{4}_V00$", base)
        p, estado, sujeto, visita = "MELIPOP", "declarado (DAT-03)", base[1:4], "V00"
        coloc, proc = "cadera derecha (declarado, PRO-04)", "registros de prueba, 2026-10-02"
    elseif occursin(r"^C\d{4}_V03$", base)
        p, estado, sujeto, visita = "COHORTE_A", "declarado (envío JMJC, 2026-10-30)", base[1:5], "V03"
        coloc, proc = "cadera derecha (declarado, PRO-04)", "JMJC, 2026-10-30"
    elseif occursin(r"^C\d{4}_V00$", base)
        estado, sujeto, visita = "hipótesis: solo el nombre (¿COHORTE_A V00?); pendiente de JMJC", base[1:5], "V00"
        coloc, proc = "por confirmar", "JMJC, 2026-11-03"
    elseif startswith(base, "H-01-")
        p, estado, sujeto = "MiopiaUNAV", "declarado (JMJC, 2026-11-03)", base
        coloc = "por fuera de la ropa (declarado); sitio corporal por confirmar"
        proc = fmt == "agd" ? "JMJC, 2026-11-03" : "JMJC, 2026-10-30"
    elseif startswith(base, "VTL")
        p, estado, sujeto = "pruebas", "declarado (DAT-03: adulto, dispositivo en un casco)", base
        coloc, proc = "casco (declarado, DAT-03)", "registros de prueba, 2026-10-02"
    elseif fmt == "ggir_config"
        occursin(r"(?i)cohorte_a", base) && (p = "COHORTE_A")
        occursin(r"(?i)unav", base) && (p = "MiopiaUNAV")
        p == "pendiente" || (estado = "declarado (nombre del archivo)")
        proc = "JMJC, 2026-11-03"
    elseif fmt == "xlsx"
        # «EJEMPLO_BBDD_Datos acel»: ActiLife + GGIR + luz de los cinco H-01 (CLAUDE.md, 2026-10-30).
        occursin(r"(?i)cohorte_a", base) && (p = "COHORTE_A")
        occursin(r"(?i)datos acel", base) && (p = "MiopiaUNAV")
        p == "pendiente" || (estado = "declarado (envío JMJC, 2026-10-30)")
        proc = "JMJC, 2026-10-30"
    end
    sub = SUBCARPETAS[fmt]
    fmt == "xlsx" && (sub = joinpath(sub, "2026-10-30_JMJC"))
    p == "pendiente" && fmt == "agd" && occursin(r"^C\d{4}_V00$", base) && (sub = "agd_C60xx_V00")
    destino = joinpath(getfield(PROYECTOS, Symbol(p)), sub, basename(nombre))
    return (; proyecto = p, estado, sujeto, visita, formato = fmt, colocacion = coloc, procedencia = proc, destino)
end

"""
    raw_files(formato; proyecto = nothing, raiz = raw_dir()) -> Vector{String}

Rutas de los archivos originales de `formato` (véase [`clasificar`](@ref)), ordenadas por
nombre, en cualquier subcarpeta. Con `proyecto` (`"COHORTE_A"`, `"MELIPOP"`…), solo los que
[`clasificar`](@ref) atribuye a ese proyecto. Se omiten los archivos y carpetas ocultos.
"""
function raw_files(fmt::AbstractString; proyecto = nothing, raiz::AbstractString = raw_dir())
    isdir(raiz) || return String[]
    out = String[]
    for (dir, dirs, fs) in walkdir(raiz)
        filter!(d -> !startswith(d, "."), dirs)
        for f in fs
            startswith(f, ".") && continue
            formato(f) == fmt || continue
            proyecto === nothing || clasificar(f).proyecto == proyecto || continue
            push!(out, joinpath(dir, f))
        end
    end
    return sort(out; by = basename)
end

"""
    find_raw(nombre; raiz = raw_dir()) -> Union{String,Nothing}

Ruta del archivo original `nombre` en cualquier subcarpeta. Error si hay más de uno con ese
nombre (los nombres de los originales son únicos en `data/raw/`).
"""
function find_raw(nombre::AbstractString; raiz::AbstractString = raw_dir())
    isdir(raiz) || return nothing
    hits = String[]
    for (dir, dirs, fs) in walkdir(raiz)
        filter!(d -> !startswith(d, "."), dirs)
        nombre in fs && push!(hits, joinpath(dir, nombre))
    end
    length(hits) > 1 && error("Hay $(length(hits)) archivos llamados $nombre en data/raw/")
    return isempty(hits) ? nothing : only(hits)
end

end # module
