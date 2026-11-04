"""
    Light

Señal de luz del ActiGraph (hito H11, cap. 12): serie por *epoch*, control de calidad por
estados y métricas de exposición, también en ventanas ancladas al sueño.

La luz se trabaja en **unidades brutas** (`lux_raw`): la conversión a lux del wGT3X-BT no
está resuelta (cap. 5). Los umbrales por defecto (100 y 1000 brutos) serían unos 125 y
1250 lux si valiera el factor 1,25 del GT3X+, lo que no está garantizado.

Cómo se reconstruye la serie:
- El dispositivo solo escribe la luz cuando cambia, así que cada valor vale hasta el
  registro siguiente (función escalonada, resolución de 1 s).
- Durante los huecos (*idle sleep* y sin datos) la luz es **desconocida**: casi no se
  registra mientras el dispositivo duerme (cap. 6). No se rellena.

Estados por *epoch* (`light_states`), en este orden de prioridad:
1. `retirado`: no uso probable (la luz no es la del niño);
2. `sueño`: dentro del periodo de reposo nocturno (SPT o SPT*, cap. 9);
3. `sin dato`: menos de la mitad del *epoch* con luz conocida;
4. `oscuro y activo` / `oscuro e inactivo`: luz media < `dark` estando despierto, con o
   sin movimiento (ENMO ≥ `active_mg`). «Oscuro y activo» reúne la oclusión del sensor por
   la ropa y los interiores por debajo del umbral de detección (unos 50 lux): con estos
   datos no se pueden separar;
5. `válido`: el resto.
"""
module Light

using Arrow, DataFrames, Dates

export light_epochs, light_states, light_daily, light_summary, save_light, load_light, ESTADOS_LUZ

const ESTADOS_LUZ = ("válido", "oscuro y activo", "oscuro e inactivo", "sin dato", "sueño", "retirado")
const MS_H = 3_600_000

"""
    light_epochs(lux, gaps; origin, stop, epoch = Minute(1), thr = (100, 1000)) -> DataFrame

Serie de luz por *epoch* desde `origin` hasta `stop`: `start`, `lux_obs` (fracción del
*epoch* con luz conocida), `lux_media` y `lux_max` (brutos, sobre los segundos conocidos;
`NaN` si ninguno) y `s_bajo`, `s_alto` (segundos conocidos con luz ≥ `thr[1]` y ≥ `thr[2]`:
en la cadera las exposiciones intensas duran segundos y la media del minuto las diluye). Cada registro vale hasta el siguiente; los segundos que tocan un hueco
o anteriores al primer registro son desconocidos.
"""
function light_epochs(lux::AbstractDataFrame, gaps::AbstractDataFrame; origin::DateTime, stop::DateTime,
                      epoch::Period = Minute(1), thr = (100, 1000))
    S = Dates.value(stop - origin) ÷ 1000 + 1
    val = fill(NaN, S)
    sec(x) = fld(Dates.value(x - origin), 1000) + 1
    t = lux.time; v = Float64.(lux.lux_raw)
    for i in eachindex(t)
        a = max(sec(t[i]), 1)
        b = min(i < length(t) ? sec(t[i+1]) - 1 : S, S)
        a <= b && (val[a:b] .= v[i])
    end
    for g in eachrow(gaps)
        a = max(sec(g.start), 1); b = min(cld(Dates.value(g.stop - origin), 1000), S)
        a <= b && (val[a:b] .= NaN)
    end
    w = Dates.value(convert(Second, epoch))
    nb = cld(S, w)
    obs = zeros(nb); media = fill(NaN, nb); mx = fill(NaN, nb); sb = zeros(Int, nb); sa = zeros(Int, nb)
    for k in 1:nb
        seg = view(val, (k-1)*w+1:min(k*w, S))
        n = 0; s = 0.0; m = -Inf
        for x in seg
            isnan(x) && continue
            n += 1; s += x; m = max(m, x)
            x >= thr[1] && (sb[k] += 1); x >= thr[2] && (sa[k] += 1)
        end
        obs[k] = n / w
        n > 0 && (media[k] = s / n; mx[k] = m)
    end
    return DataFrame(start = [origin + Second((k - 1) * w) for k in 1:nb], lux_obs = obs, lux_media = media, lux_max = mx,
                     s_bajo = sb, s_alto = sa)
end

# ¿Cae t en alguno de los intervalos [a, b)?
_dentro(t, iv) = any(a <= t < b for (a, b) in iv)

"""
    light_states(le, ep, noches; dark = 1, active_mg = 63.3, min_obs = 0.5) -> DataFrame

Une la luz por *epoch* (`le`) con los *epochs* de aceleración del mismo tamaño (`ep`, con
`enmo_idle0` y `nonwear`) y los periodos de reposo nocturno (`noches` de `Sleep.sleep_nights`), y
asigna a cada *epoch* uno de `ESTADOS_LUZ`. `active_mg`: umbral de actividad ligera en
cadera (Hildebrand, 63,3 mg).
"""
function light_states(le::AbstractDataFrame, ep::AbstractDataFrame, noches::AbstractDataFrame;
                      dark::Real = 1, active_mg::Real = 63.3, min_obs::Real = 0.5)
    d = innerjoin(le, select(ep, :start, :enmo_idle0, :nonwear), on = :start)
    sort!(d, :start)
    iv = [(n.inicio, n.fin) for n in eachrow(noches) if n.metodo != "sin ventana"]
    d.estado = map(eachrow(d)) do r
        r.nonwear && return "retirado"
        _dentro(r.start, iv) && return "sueño"
        r.lux_obs < min_obs && return "sin dato"
        if r.lux_media < dark
            return (isfinite(r.enmo_idle0) && r.enmo_idle0 >= active_mg) ? "oscuro y activo" : "oscuro e inactivo"
        end
        return "válido"
    end
    return d
end

_mean(v) = isempty(v) ? NaN : sum(v) / length(v)
_h(t, d) = Dates.value(t - DateTime(d)) / MS_H

# Media circular de horas del día (0–24).
function _hora_media(h)
    isempty(h) && return NaN
    s = sum(sin.(2π .* h ./ 24)); c = sum(cos.(2π .* h ./ 24))
    return mod(atan(s, c) * 24 / 2π, 24)
end

# Número de tramos de al menos `minlen` epochs consecutivos con v = true.
function _tramos(v, minlen)
    n = 0; run = 0
    for x in v
        if x
            run += 1
        else
            run >= minlen && (n += 1); run = 0
        end
    end
    run >= minlen && (n += 1)
    return n
end

# Resumen de luz en un subconjunto de epochs (índices `k` de `d`).
function _ventana(d, k, thr)
    v = k[d.estado[k] .== "válido"]
    return (min_valido = length(v), pct_valido = isempty(k) ? NaN : 100 * length(v) / length(k),
            media_bruta = _mean(d.lux_media[v]), min_sobre_bajo = count(d.lux_media[v] .>= thr[1]))
end

"""
    light_daily(d, noches; thr = (100, 1000), mvpa_mg = 142.6, bout_min = 5,
                tras_despertar = Hour(2), antes_dormir = Hour(3)) -> DataFrame

Una fila por día natural con tiempo despierto y con el dispositivo puesto (`d` de
`light_states`, *epochs* de 60 s):
- calidad: minutos despierto con el dispositivo y % en cada estado;
- exposición (*epochs* válidos): minutos con media ≥ `thr[1]` y ≥ `thr[2]` brutos; además,
  en todo el tiempo despierto con luz conocida, minutos acumulados segundo a segundo
  (`seg_*`, si `d` trae `s_bajo` y `s_alto`); primera, última y
  hora media (circular) con luz ≥ `thr[2]`; entradas en luz ≥ `thr[2]` de al menos
  `bout_min` minutos (transiciones de baja a alta);
- integración: minutos de AFMV (ENMO ≥ `mvpa_mg`) según la luz: oscuro, interior
  (`dark` a `thr[2]`), compatible con exterior (≥ `thr[2]`) o sin luz válida;
- ventanas ancladas: `tras_despertar` desde el fin del periodo de reposo nocturno anterior y
  `antes_dormir` hasta el inicio de la siguiente (minutos válidos, % válido, media bruta y
  minutos ≥ `thr[1]`). `NaN` si falta el ancla.
"""
function light_daily(d::AbstractDataFrame, noches::AbstractDataFrame; thr = (100, 1000), mvpa_mg::Real = 142.6,
                     bout_min::Integer = 5, tras_despertar::Period = Hour(2), antes_dormir::Period = Hour(3))
    dias = sort(unique(Date.(d.start)))
    ok_n = noches[noches.metodo .!= "sin ventana", :]
    filas = NamedTuple[]
    for dia in dias
        k = findall(t -> Date(t) == dia, d.start)
        desp = k[(d.estado[k] .!= "sueño") .& (d.estado[k] .!= "retirado")]
        isempty(desp) && continue
        st = d.estado[desp]
        v = desp[st .== "válido"]
        alto = v[d.lux_media[v] .>= thr[2]]
        halto = _h.(d.start[alto], dia)
        mv = desp[isfinite.(d.enmo_idle0[desp]) .& (d.enmo_idle0[desp] .>= mvpa_mg)]
        ml = d.lux_media[mv]; mst = d.estado[mv]
        # Anclas: fin de la ventana de la noche anterior e inicio de la de esta noche.
        kf = findfirst(==(dia - Day(1)), ok_n.noche); ki = findfirst(==(dia), ok_n.noche)
        despierta = kf === nothing ? nothing : ok_n.fin[kf]
        duerme = ki === nothing ? nothing : ok_n.inicio[ki]
        vd = despierta === nothing ? nothing : _ventana(d, findall(t -> despierta <= t < despierta + tras_despertar, d.start), thr)
        va = duerme === nothing ? nothing : _ventana(d, findall(t -> duerme - antes_dormir <= t < duerme, d.start), thr)
        nanv = (min_valido = NaN, pct_valido = NaN, media_bruta = NaN, min_sobre_bajo = NaN)
        vd === nothing && (vd = nanv); va === nothing && (va = nanv)
        alto_mask = [d.estado[i] == "válido" && d.lux_media[i] >= thr[2] for i in k]
        push!(filas, (dia = dia, dia_semana = dayofweek(dia), fin_de_semana = dayofweek(dia) >= 6,
                      min_despierto = length(desp),
                      pct_valido = 100 * count(==("válido"), st) / length(desp),
                      pct_oscuro_activo = 100 * count(==("oscuro y activo"), st) / length(desp),
                      pct_oscuro_inactivo = 100 * count(==("oscuro e inactivo"), st) / length(desp),
                      pct_sin_dato = 100 * count(==("sin dato"), st) / length(desp),
                      min_valido = length(v),
                      min_sobre_bajo = count(d.lux_media[v] .>= thr[1]), min_sobre_alto = length(alto),
                      seg_sobre_bajo = hasproperty(d, :s_bajo) ? sum(d.s_bajo[desp]) / 60 : NaN,
                      seg_sobre_alto = hasproperty(d, :s_alto) ? sum(d.s_alto[desp]) / 60 : NaN,
                      primera_alta_h = isempty(halto) ? NaN : minimum(halto),
                      ultima_alta_h = isempty(halto) ? NaN : maximum(halto),
                      hora_media_alta_h = _hora_media(halto),
                      entradas_alta = _tramos(alto_mask, bout_min),
                      afmv_min = length(mv),
                      afmv_oscuro = count((mst .!= "válido") .& (mst .!= "sin dato")),
                      afmv_interior = count((mst .== "válido") .& (ml .< thr[2])),
                      afmv_exterior = count((mst .== "válido") .& (ml .>= thr[2])),
                      afmv_sin_luz = count(mst .== "sin dato"),
                      despertar_h = despierta === nothing ? NaN : _h(despierta, dia),
                      tras_despertar_pct_valido = vd.pct_valido, tras_despertar_media = vd.media_bruta,
                      tras_despertar_min_bajo = vd.min_sobre_bajo,
                      dormir_h = duerme === nothing ? NaN : _h(duerme, dia),
                      antes_dormir_pct_valido = va.pct_valido, antes_dormir_media = va.media_bruta,
                      antes_dormir_min_bajo = va.min_sobre_bajo))
    end
    return DataFrame(filas)
end

_sd(v) = length(v) < 2 ? NaN : sqrt(sum((v .- _mean(v)) .^ 2) / (length(v) - 1))
_finite(v) = filter(isfinite, collect(Float64, v))

"""
    light_summary(dd; min_despierto = 600) -> NamedTuple

Resumen por registro de `light_daily` en los días con al menos `min_despierto` minutos
despierto con el dispositivo puesto: calidad media, exposición, regularidad (desviación
estándar entre días de la hora media con luz alta), diario frente a fin de semana, AFMV
por luz y ventanas ancladas.
"""
function light_summary(dd::AbstractDataFrame; min_despierto::Integer = 600)
    x = dd[dd.min_despierto .>= min_despierto, :]
    m(c) = _mean(_finite(x[!, c]))
    lab = x[.!x.fin_de_semana, :]; fds = x[x.fin_de_semana, :]
    hm = _finite(x.hora_media_alta_h)
    return (dias = nrow(x), pct_valido = m(:pct_valido), pct_oscuro_activo = m(:pct_oscuro_activo),
            pct_oscuro_inactivo = m(:pct_oscuro_inactivo), pct_sin_dato = m(:pct_sin_dato),
            min_sobre_bajo = m(:min_sobre_bajo), min_sobre_alto = m(:min_sobre_alto),
            seg_sobre_bajo = m(:seg_sobre_bajo), seg_sobre_alto = m(:seg_sobre_alto),
            min_sobre_alto_diario = _mean(_finite(lab.min_sobre_alto)), min_sobre_alto_finde = _mean(_finite(fds.min_sobre_alto)),
            dias_con_alta = length(hm), hora_media_alta_h = _hora_media(hm),
            hora_media_alta_de_h = length(hm) < 2 ? NaN : _sd([mod(h - _hora_media(hm) + 12, 24) - 12 for h in hm]),
            entradas_alta = m(:entradas_alta), afmv_min = m(:afmv_min), afmv_exterior = m(:afmv_exterior),
            afmv_oscuro = m(:afmv_oscuro), afmv_sin_luz = m(:afmv_sin_luz),
            tras_despertar_pct_valido = m(:tras_despertar_pct_valido), tras_despertar_min_bajo = m(:tras_despertar_min_bajo),
            antes_dormir_pct_valido = m(:antes_dormir_pct_valido), antes_dormir_min_bajo = m(:antes_dormir_min_bajo))
end

"""
    save_light(dir, d) / load_light(dir) -> DataFrame

Guarda o lee la luz por *epoch* de 60 s con su estado en `dir/luz_60s.arrow`.
"""
save_light(dir::AbstractString, d::AbstractDataFrame) = Arrow.write(joinpath(dir, "luz_60s.arrow"), d; compress = :zstd)
load_light(dir::AbstractString) = DataFrame(Arrow.Table(joinpath(dir, "luz_60s.arrow")); copycols = true)

end # module
