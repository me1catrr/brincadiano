# 23 · AGD de ActiLife: inventario, lectura con lista blanca y correspondencia con los .gt3x.
#
# Qué hace
#   1. Inventario: busca todos los .agd bajo data/raw/ (también en subcarpetas), con tamaño y
#      SHA-256; marca las copias idénticas (mismo hash) y lee solo la primera.
#   2. Lectura (Brincadiano.AGD, sqlite3 en solo lectura): metadatos de la lista blanca,
#      epochs (counts, pasos, lux, inclinómetro), filtros de tiempo y recuento de tablas.
#   3. Correspondencia con cada .gt3x leído en data/interim/: número de serie, inicio, fin y
#      descarga. Un par se da por CONFIRMADO (metadatos) si coinciden la serie y el inicio y la
#      descarga a ±1 min; «mismo dispositivo, otra grabación» si coincide solo la serie. Los pares
#      por nombre (mismo código numérico) se listan siempre, coincidan o no: así se ve la
#      discrepancia C60xx_V00 / 6005_V00 / C60xx_V03 sin asociarlos de antemano.
#   4. Señal, solo en los pares confirmados: counts por minuto del AGD frente a la ENMO por minuto
#      del pipeline (Spearman, desfase de −3 a +3 min), y lux del AGD frente a lux_raw del
#      log.bin segundo a segundo, solo donde el pipeline conoce la luz (fuera de huecos).
#   El proyecto de cada ID sale del nombre (Brincadiano.Rutas.clasificar: COHORTE_A, MELIPOP,
#   MiopiaUNAV, pruebas o pendiente) y se escribe en outputs/inventario/ids_cohorte.csv para
#   `17 --cohortes`. Es lo declarado por el grupo, no una verificación; los C60xx_V00, «pendiente».
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/23_agd.jl            # tras 01 (o el lote 16) y 05
#
# Depende de
#   01/16 (meta, gaps y lux) y 05 (epochs_60s) para los pares; sqlite3 del sistema.
#
# Salida (no se versiona)
#   data/interim/agd/<archivo>/{epochs,filters,meta}.arrow
#   outputs/agd/tables/agd_archivos.csv          inventario: ruta, bytes, SHA-256, duplicados
#   outputs/agd/tables/agd_meta.csv              metadatos (lista blanca) y contenido por archivo
#   outputs/agd/tables/agd_correspondencia.csv   cada AGD frente a cada .gt3x candidato
#   outputs/agd/tables/agd_senal.csv             counts frente a ENMO y lux frente a lux_raw
#   outputs/inventario/ids_cohorte.csv           id, cohorte (por nombre; hipótesis)
#
# Consola: solo el nombre del archivo (ID pseudonimizado), recuentos y el veredicto.

using Brincadiano, Brincadiano.AGD, Brincadiano.GT3X, Brincadiano.Light
using Brincadiano.Rutas: Rutas, raw_files, clasificar, id_registro
using CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = Rutas.raw_dir()
const INTERIM = Rutas.interim_dir()
const INTERIM_AGD = joinpath(INTERIM, "agd")
const OUT = joinpath(ROOT, "outputs", "agd", "tables")
const TOL = Minute(1)

# Proyecto por el nombre (Brincadiano.Rutas.clasificar): los C60xx_V00 quedan como «pendiente».
conjunto(id) = clasificar(id).proyecto
# «H-01-0011sec» → «H-01-001»: ActiLife añade la duración del epoch al nombre.
id_agd(archivo) = id_registro(archivo)
codigo(id) = (m = match(r"(\d{3,4})(?:_V\d+)?$", id); m === nothing ? "" : m.captures[1])

sha256(path) = first(split(read(`shasum -a 256 $path`, String)))
minutos(a, b) = (a === nothing || b === nothing || ismissing(a) || ismissing(b)) ? missing :
    round(Dates.value(a - b) / 60_000; digits = 1)

function rangos(x)
    p = sortperm(x); r = similar(x, Float64); i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[p[j + 1]] == x[p[i]]; j += 1; end
        r[p[i:j]] .= (i + j) / 2; i = j + 1
    end
    return r
end
spearman(a, b) = length(a) < 3 ? missing : cor(rangos(a), rangos(b))

# ---------------------------------------------------------------- 1 · inventario
archivos = raw_files("agd")
isempty(archivos) && error("No hay archivos .agd en data/raw/")
inv = DataFrame(archivo = basename.(archivos), ruta = relpath.(archivos, ROOT),
                bytes = filesize.(archivos), sha256 = sha256.(archivos))
inv.duplicado_de = Vector{Union{String,Missing}}(missing, nrow(inv))
for (i, h) in enumerate(inv.sha256)
    j = findfirst(==(h), inv.sha256)
    j < i && (inv.duplicado_de[i] = inv.ruta[j])
end
mkpath(OUT)
CSV.write(joinpath(OUT, "agd_archivos.csv"), inv)
@printf("AGD: %d archivos, %d únicos por SHA-256\n", nrow(inv), count(ismissing, inv.duplicado_de))

# ---------------------------------------------------------------- 2 · lectura
gt3x_ids = sort(filter(d -> isfile(joinpath(INTERIM, d, "meta.arrow")), readdir(INTERIM)))
leer(id, n) = DataFrame(GT3X.Arrow.Table(joinpath(INTERIM, id, n * ".arrow")); copycols = true)
gmeta = Dict(id => leer(id, "meta")[1, :] for id in gt3x_ids)

metas = NamedTuple[]; corr = NamedTuple[]; senal = NamedTuple[]
for r in eachrow(inv[ismissing.(inv.duplicado_de), :])
    a = read_agd(joinpath(ROOT, r.ruta)); i = a.info; e = a.epochs
    id = id_agd(r.archivo)
    dst = joinpath(INTERIM_AGD, splitext(r.archivo)[1]); mkpath(dst)
    GT3X.Arrow.write(joinpath(dst, "epochs.arrow"), e)
    GT3X.Arrow.write(joinpath(dst, "filters.arrow"), a.filters)
    fila_meta = (; archivo = r.archivo, id_agd = id, conjunto = conjunto(id), serial = i.serial,
        device = i.device, firmware = i.firmware, software = i.software, software_version = i.software_version,
        filtro = i.filter, modo = something(i.mode, missing), epoch_s = i.epoch_s,
        epochs = nrow(e), epochs_declarados = something(i.epoch_count, missing),
        fs_original = something(i.sample_rate, missing), inicio = i.start, fin = something(i.stop, missing),
        descarga = something(i.download, missing), bateria_v = something(i.battery_voltage, missing),
        terminado = something(i.finished, missing), columnas = join(names(e)[2:end], " "),
        frac_axis1_pos = "axis1" in names(e) ? round(mean(coalesce.(e.axis1, 0) .> 0); digits = 4) : missing,
        frac_lux_pos = "lux" in names(e) ? round(mean(coalesce.(e.lux, 0) .> 0); digits = 4) : missing,
        lux_max = "lux" in names(e) ? maximum(skipmissing(e.lux); init = 0.0) : missing,
        lux_nulos = "lux" in names(e) ? count(ismissing, e.lux) : missing,
        filtros = nrow(a.filters),
        filtros_min = round(sum(Dates.value.(a.filters.stop .- a.filters.start); init = 0) / 60_000; digits = 1),
        filas_capsense = a.rows["capsense"], filas_sleep = a.rows["sleep"],
        filas_diario = a.rows["logDiaryTimes"], filas_proximity = a.rows["proximity"])
    push!(metas, fila_meta)
    GT3X.Arrow.write(joinpath(dst, "meta.arrow"), DataFrame([fila_meta]))

    # ------------------------------------------------------------ 3 · correspondencia
    cand = filter(g -> gmeta[g].serial == i.serial || codigo(g) == codigo(id), gt3x_ids)
    isempty(cand) && push!(corr, (; archivo = r.archivo, id_agd = id, conjunto_agd = conjunto(id), gt3x = missing,
        conjunto_gt3x = missing, por_nombre = false, misma_serie = false, d_inicio_min = missing, d_fin_min = missing,
        d_descarga_min = missing, veredicto = "sin .gt3x candidato"))
    for g in cand
        m = gmeta[g]
        serie = m.serial == i.serial
        di = minutos(i.start, m.start); dd = minutos(i.download, m.download)
        ok = serie && !ismissing(di) && abs(di) <= 1 && !ismissing(dd) && abs(dd) <= 1
        ver = ok ? "CONFIRMADO (metadatos)" : serie ? "mismo dispositivo, otra grabación" :
              "DISCREPANCIA: mismo código, otro dispositivo"
        push!(corr, (; archivo = r.archivo, id_agd = id, conjunto_agd = conjunto(id), gt3x = g,
            conjunto_gt3x = conjunto(g), por_nombre = codigo(g) == codigo(id), misma_serie = serie,
            d_inicio_min = di, d_fin_min = minutos(i.stop, m.stop), d_descarga_min = dd, veredicto = ver))

        # -------------------------------------------------------- 4 · señal (solo confirmados)
        ok || continue
        ep = leer(g, "epochs_60s")
        e.minuto = floor.(e.time, Minute)
        e.vm = sqrt.(coalesce.(e.axis1, 0.0) .^ 2 .+ coalesce.(e.axis2, 0.0) .^ 2 .+ coalesce.(e.axis3, 0.0) .^ 2)
        cpm = combine(groupby(e, :minuto), :axis1 => (x -> sum(coalesce.(x, 0.0))) => :cpm_axis1, :vm => sum => :cpm_vm)
        res = Dict{Int,Any}()
        for lag in -3:3
            j = innerjoin(transform(cpm, :minuto => ByRow(t -> t + Minute(lag)) => :start), ep; on = :start)
            j = j[.!ismissing.(j.enmo_idle0), :]
            res[lag] = (n = nrow(j), rho = spearman(j.cpm_vm, Float64.(j.enmo_idle0)),
                        rho_axis1 = spearman(j.cpm_axis1, Float64.(j.enmo_idle0)))
        end
        mejor = argmax(l -> coalesce(res[l].rho, -Inf), -3:3)
        # Luz: segundo a segundo, solo donde el pipeline la conoce (fuera de huecos).
        lx = leer(g, "lux"); gp = leer(g, "gaps")
        s1 = light_epochs(lx, gp; origin = i.start, stop = i.start + Second(nrow(e) - 1), epoch = Second(1))
        lagd = Dict(t => v for (t, v) in zip(e.time, e.lux))
        x = Float64[]; y = Float64[]
        for q in eachrow(s1)
            isnan(q.lux_media) && continue
            v = get(lagd, q.start, missing); ismissing(v) && continue
            push!(x, q.lux_media); push!(y, v)
        end
        pos = (x .> 0) .& (y .> 0)
        push!(senal, (; archivo = r.archivo, gt3x = g, minutos = res[0].n, rho_vm_enmo = res[0].rho,
            rho_axis1_enmo = res[0].rho_axis1, desfase_mejor_min = mejor, rho_mejor = res[mejor].rho,
            seg_luz_conocida = length(x), frac_seg_luz_conocida = round(length(x) / nrow(e); digits = 4),
            frac_ambos_cero = isempty(x) ? missing : round(mean((x .== 0) .& (y .== 0)); digits = 4),
            frac_cero_discrepante = isempty(x) ? missing : round(mean((x .== 0) .⊻ (y .== 0)); digits = 4),
            rho_lux = spearman(x, y),
            razon_lux_mediana = any(pos) ? round(median(y[pos] ./ x[pos]); digits = 4) : missing,
            razon_lux_p05 = any(pos) ? round(quantile(y[pos] ./ x[pos], 0.05); digits = 4) : missing,
            razon_lux_p95 = any(pos) ? round(quantile(y[pos] ./ x[pos], 0.95); digits = 4) : missing,
            pendiente_origen = any(pos) ? round(sum(x[pos] .* y[pos]) / sum(x[pos] .^ 2); digits = 4) : missing,
            lux_raw_max = isempty(x) ? missing : maximum(x), lux_agd_max = isempty(y) ? missing : maximum(y)))
    end
    v = [c.veredicto for c in corr if c.archivo == r.archivo]
    @printf("  %-22s %-20s epochs %6d  %s\n", r.archivo, conjunto(id), nrow(e), join(unique(v), " · "))
end

CSV.write(joinpath(OUT, "agd_meta.csv"), DataFrame(metas))
CSV.write(joinpath(OUT, "agd_correspondencia.csv"), DataFrame(corr))
isempty(senal) || CSV.write(joinpath(OUT, "agd_senal.csv"), DataFrame(senal))

ids = sort(unique([gt3x_ids; [m.id_agd for m in metas]]))
mkpath(joinpath(ROOT, "outputs", "inventario"))
CSV.write(joinpath(ROOT, "outputs", "inventario", "ids_cohorte.csv"), DataFrame(id = ids, cohorte = conjunto.(ids)))

c = DataFrame(corr)
for v in unique(c.veredicto)
    @printf("%-45s %d\n", v, count(==(v), c.veredicto))
end
println("Salida en ", relpath(OUT, ROOT))
