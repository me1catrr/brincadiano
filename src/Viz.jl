"""
    Viz

Gráficas del proyecto. Aquí viven las funciones **puras** que preparan los datos para
dibujar (agregados por intervalo, ejes en días, rupturas en los huecos), que se pueden
testear sin Makie. Las funciones de dibujo solo se declaran; su implementación está en
la extensión `BrincadianoCairoMakieExt` y se activa con `using CairoMakie`.

Los resúmenes por intervalo de este módulo son **de visualización**, no métricas del
estudio: las métricas (ENMO, MAD, epochs) son de `Metrics` (hito H5).

Privacidad: los ejes temporales se rotulan como «día k (día de la semana)», nunca con
la fecha de calendario del registro.
"""
module Viz

using DataFrames, Dates

export vector_norm, binned_summary, days_since, day_labels, break_at_gaps
export ANCHO_TEXTO, ANCHO_COMPLETO, ALTO_MAX, figure_size, tema_informe, guardar, procedencia_tex
export plot_overview, plot_raw_window, plot_days24h, plot_histograms, plot_metrics_window
export plot_profile24h, plot_timing_circular, plot_sleep_nights, plot_light_days
export plot_actogram, plot_profiles_compare, plot_enmo_ecdf, plot_timing_compare, plot_sleep_raster,
       plot_vigorous_summary, plot_circular_plane, plot_bland_altman, plot_agreement_timeline,
       plot_timing_sensitivity, plot_rhythm_spectrum, plot_harmonic_profile, plot_psd_classes,
       plot_cadence, plot_episode_spectral_qc, plot_qc_days, plot_capsense_runs, plot_anchor_nights,
       plot_timing_days, plot_clock_coverage, plot_threshold_sweep, plot_light_anchored,
       plot_mvpa_light, plot_light_enmo_density, plot_light_states_by_hour, plot_light_profile,
       plot_imputation_dumbbell, plot_intensity_ggir, plot_nights_subsample

"""
    vector_norm(x, y, z) -> Vector{Float32}

Longitud del vector de aceleración, `√(x² + y² + z²)`, muestra a muestra (en g).
"""
vector_norm(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    Float32[sqrt(Float64(a)^2 + Float64(b)^2 + Float64(c)^2) for (a, b, c) in zip(x, y, z)]
vector_norm(acc::AbstractDataFrame) = vector_norm(acc.x, acc.y, acc.z)

"""
    binned_summary(t, v, period; origin = floor(first(t), Day), stop = last(t)) -> DataFrame

Resume `v` (muestras en los instantes `t`, crecientes) en intervalos consecutivos de
duración `period` a partir de `origin`. Devuelve **todos** los intervalos entre `origin` y
`stop`, también los vacíos, con columnas `start`, `n` (muestras), `vmin`, `vmax`, `vmean`.
En los intervalos sin muestras los tres últimos son `NaN`: así Makie corta las líneas en
los huecos en lugar de unir los extremos (los huecos nunca se rellenan).
"""
function binned_summary(t::AbstractVector{DateTime}, v::AbstractVector{<:Real}, period::Period;
                        origin::DateTime = floor(first(t), Day), stop::DateTime = last(t))
    length(t) == length(v) || throw(DimensionMismatch("t y v deben tener la misma longitud"))
    w = Dates.value(convert(Millisecond, period))
    w > 0 || throw(ArgumentError("period debe ser positivo"))
    nb = fld(Dates.value(stop - origin), w) + 1
    n = zeros(Int, nb)
    vmin = fill(Inf, nb); vmax = fill(-Inf, nb); vsum = zeros(nb)
    @inbounds for (ti, vi) in zip(t, v)
        k = fld(Dates.value(ti - origin), w) + 1
        1 <= k <= nb || continue
        n[k] += 1
        x = Float64(vi)
        vmin[k] = min(vmin[k], x); vmax[k] = max(vmax[k], x); vsum[k] += x
    end
    empty = n .== 0
    vmin[empty] .= NaN; vmax[empty] .= NaN
    vmean = ifelse.(empty, NaN, vsum ./ max.(n, 1))
    return DataFrame(start = origin .+ Millisecond(w) .* (0:nb-1), n = n,
                     vmin = vmin, vmax = vmax, vmean = vmean)
end

"""
    days_since(t, origin) -> Float64 o Vector{Float64}

Tiempo en días (fraccionarios) desde `origin`. Es el eje x de las vistas de varios días.
"""
days_since(t::DateTime, origin::DateTime) = Dates.value(t - origin) / 86_400_000
days_since(t::AbstractVector{DateTime}, origin::DateTime) = [days_since(ti, origin) for ti in t]

const DIAS_SEMANA = ("lun", "mar", "mié", "jue", "vie", "sáb", "dom")

"""
    day_labels(origin, ndays) -> Vector{String}

Rótulos «día k (lun)» para los días `0:ndays-1` desde `origin` (medianoche del primer
día). No incluye la fecha de calendario.
"""
day_labels(origin::DateTime, ndays::Integer) =
    ["día $(k + 1) ($(DIAS_SEMANA[dayofweek(origin + Day(k))]))" for k in 0:ndays-1]

"""
    break_at_gaps(t, ys...; maxstep) -> (t′, ys′...)

Inserta un `NaN` (y su instante) allí donde dos muestras consecutivas están separadas más
de `maxstep`, para que las líneas no crucen los huecos. `t` en segundos (`Real`).
"""
function break_at_gaps(t::AbstractVector{<:Real}, ys::AbstractVector{<:Real}...; maxstep::Real)
    cut = findall(>(maxstep), diff(t))
    tt = Float64.(t); yy = [Float64.(y) for y in ys]
    for i in reverse(cut)
        insert!(tt, i + 1, (tt[i] + tt[i + 1]) / 2)
        for y in yy
            insert!(y, i + 1, NaN)
        end
    end
    return (tt, yy...)
end

# ─── Tamaños de figura en puntos (página A4 apaisada del informe, report/main.tex) ──────
#
# Las figuras se dibujan a su tamaño final: 1 unidad de Makie = 1 pt del PDF (`guardar` usa
# `pt_per_unit = 1`), así el texto de la figura mide lo mismo que el del informe.

"Ancho del texto principal: 160 mm, en puntos."
const ANCHO_TEXTO = 160 * 72 / 25.4
"Ancho completo (texto + separación + margen, 160 + 8 + 90 mm), para `figure*`."
const ANCHO_COMPLETO = 258 * 72 / 25.4
"Alto máximo de una figura para que quepa con su pie y su fuente en la caja de 180 mm (510 pt)."
const ALTO_MAX = 440.0

"""
    figure_size(ancho, alto) -> (w, h)

Tamaño en puntos: `ancho = :texto` (`ANCHO_TEXTO`) o `:completo` (`ANCHO_COMPLETO`); el alto
se limita a `ALTO_MAX` para que la figura y su pie quepan en la página.
"""
function figure_size(ancho::Symbol, alto::Real)
    w = ancho === :texto ? ANCHO_TEXTO : ancho === :completo ? ANCHO_COMPLETO :
        throw(ArgumentError("ancho debe ser :texto o :completo"))
    return (w, min(Float64(alto), ALTO_MAX))
end

"""
    tema_informe() -> Makie.Theme

Tema común de las figuras: la sans del informe (TeX Gyre Heros), cuerpo de 8,5 pt, líneas
finas y rejilla tenue. Todas las funciones `plot_*` lo aplican.
"""
function tema_informe end

"""
    guardar(base, fig; png = true, procedencia = nothing)

Guarda `base.pdf` a tamaño real (1 unidad = 1 pt) y, si `png`, `base.png` a 3 px por punto.

Con `procedencia = (; tipo, registros, funcion, entradas)` escribe además
`base.procedencia.tex` (ver [`procedencia_tex`](@ref)): el informe lo lee con
`\\figprocedencia` y lo muestra bajo la figura y en el anexo de procedencia. El script, la
ruta de salida, el commit y la fecha los completa `guardar`.
"""
function guardar end

const TIPOS_PROCEDENCIA = (:datos, :sintetico, :fixture)
# IDs pseudonimizados (nombre del archivo sin extensión, p. ej. 6005_V00) y rutas relativas.
const _RE_ID = r"^[A-Za-z0-9][A-Za-z0-9_-]*$"
const _RE_RUTA = r"^[A-Za-z0-9_.<>*{},/ -]+$"

_tex(s::AbstractString) = replace(s, "_" => "\\_", "%" => "\\%", "#" => "\\#", "&" => "\\&")

"""
    procedencia_tex(; tipo, registros, funcion, script, entradas, salida, commit, fecha) -> String

Ficha de procedencia de una figura para el informe: una sola llamada
`\\BRprocedencia{tipo}{registros}{función}{script}{entradas}{salida}{commit}{fecha}`.

- `tipo`: `:datos` (registros de participantes), `:sintetico` o `:fixture` (archivo público).
- `registros`: IDs pseudonimizados (vacío si no hay); `funcion`: la función de `Viz` que dibuja.
- `script`, `entradas`, `salida`: rutas **relativas** a la raíz del repositorio.

Solo admite IDs y rutas (lanza `ArgumentError` con cualquier otra cosa), de modo que la ficha
no puede llevar campos de `info.txt` ni rutas absolutas del equipo.
"""
function procedencia_tex(; tipo::Symbol, registros = String[], funcion::AbstractString,
                         script::AbstractString, entradas = String[], salida::AbstractString,
                         commit::AbstractString, fecha::Date)
    tipo in TIPOS_PROCEDENCIA || throw(ArgumentError("tipo debe ser uno de $(TIPOS_PROCEDENCIA)"))
    tipo === :datos && isempty(registros) && throw(ArgumentError("una figura con datos necesita sus registros"))
    for id in registros
        occursin(_RE_ID, id) || throw(ArgumentError("registro no válido: $(repr(id))"))
    end
    occursin(r"^[A-Za-z_][A-Za-z0-9_.!]*$", funcion) || throw(ArgumentError("función no válida: $(repr(funcion))"))
    for r in (script, salida, entradas...)
        (occursin(_RE_RUTA, r) && !startswith(r, "/") && !occursin("..", r)) ||
            throw(ArgumentError("ruta no válida (debe ser relativa al repositorio): $(repr(r))"))
    end
    occursin(r"^[0-9a-f]{7,40}\+?$|^sin git$", commit) || throw(ArgumentError("commit no válido: $(repr(commit))"))
    campos = (String(tipo), join(registros, ", "), funcion, script, join(entradas, ", "), salida,
              commit, Dates.format(fecha, dateformat"yyyy-mm-dd"))
    return "% Generado por Viz.guardar: no editar a mano.\n\\BRprocedencia" *
           join("{" * _tex(c) * "}" for c in campos) * "\n"
end

# ─── Funciones de dibujo (implementadas en BrincadianoCairoMakieExt) ───────────

"""
    plot_overview(rec; bin = Minute(1)) -> Figure

Vista de todo el registro: envolvente de ‖a‖ por intervalo con los huecos sombreados,
cobertura horaria del muestreo, luz bruta, contacto con la piel (CAPSENSE) y batería.
"""
function plot_overview end

"""
    plot_raw_window(rec, t0, dur; title = "") -> Figure

Señal a la frecuencia de muestreo (x, y, z y ‖a‖) entre `t0` y `t0 + dur`, con los
huecos sombreados.
"""
function plot_raw_window end

"""
    plot_days24h(rec; bin = Minute(1), double = false, runs = nothing) -> Figure

Una fila por día (0–24 h, o 0–48 h con `double = true`): media por intervalo de
|‖a‖ − 1 g| con los huecos sombreados y la luz bruta superpuesta. Con `runs` (tabla de
tramos inmóviles de `QC.qc_report`, con `start`, `stop` y `etiqueta`) añade una franja
inferior coloreada por etiqueta.
"""
function plot_days24h end

"""
    plot_histograms(rec) -> Figure

Distribución de ‖a‖ (todas las muestras) y de la luz bruta ponderada por el tiempo que
se mantiene cada valor.
"""
function plot_histograms end

"""
    plot_metrics_window(rec, t0, dur; epochs = (Second(1), Second(5), Second(60)),
                        offset = (0, 0, 0), scale = (1, 1, 1)) -> Figure

Señal ‖a‖ a la frecuencia de muestreo, con el área por encima de 1 g sombreada (lo que
suma ENMO), y debajo ENMO (`enmo_idle0`) agregado en cada longitud de *epoch*, todas en
la misma escala, para ver cómo los *epochs* largos diluyen los esfuerzos breves. `offset`
y `scale` aplican la calibración de `Metrics.autocalibrate`.
"""
function plot_metrics_window end

"""
    plot_profile24h(perfiles, rar; smooth = 15) -> Figure

Perfil medio de 24 h (media móvil de `smooth` min) de uno o varios perfiles de
`Circadian.profile24h`, dados como pares `(etiqueta, perfil)`, con las ventanas L5 y M10
de `rar` (de `Circadian.rar_metrics`) sombreadas y el punto medio de L5 marcado.
"""
function plot_profile24h end

"""
    plot_timing_circular(fase, dosis, dia, θ, R; ancla = "punto medio de L5") -> Figure

Diana circular del *timing*: 0 h arriba = el ancla, sentido horario = horas después.
Cada episodio es un punto en el anillo de su día (`dia` = 1, 2, …), de tamaño ∝ √dosis;
la flecha es el vector medio ponderado (dirección θ, longitud R).
"""
function plot_timing_circular end

"""
    plot_sleep_nights(ep, noches; axis, threshold = 60, title = "") -> Figure

Una fila por noche (de 12:00 a 12:00): |ángulo del eje longitudinal| por minuto (0–90°,
0 = horizontal), el umbral de «tumbado» (`threshold`), el no uso probable sombreado en
gris, el *idle sleep* heredado en ámbar claro, el periodo de reposo nocturno (SPT o SPT*) como
banda y los bloques de inactividad sostenida (SIB) dentro de ella como franja inferior.
`ep`: *epochs* de 5 s de `Sleep.fill_idle` con `nonwear`; `noches`: de `Sleep.sleep_nights`.
"""
function plot_sleep_nights end

"""
    plot_light_days(d; thr = (100, 1000), title = "") -> Figure

Una fila por día (de 00:00 a 24:00) con la luz por *epoch* de `Light.light_states`:
log₁₀(1 + luz bruta) como barras, coloreadas por estado (válido, oscuro y activo, oscuro e
inactivo), los umbrales `thr` como líneas discontinuas y, de fondo, el periodo de reposo nocturno, la
retirada y los tramos sin dato.
"""
function plot_light_days end

# ─── Figuras de síntesis (H13; datos ya procesados, scripts/10_figuras.jl) ────────────

"""
    plot_actogram(e60, luz, noches, episodios, capsense; ancho = :completo) -> Figure

Actograma multimodal: una fila por día (00–24 h) con ENMO por minuto (`e60`, barras), la
luz (`luz` de `Light.light_states`) y el contacto con la piel (`capsense`) como franjas
inferiores, el periodo de reposo nocturno (`noches`) y el no uso de fondo, y los episodios vigorosos
(`episodios`, columna `medio`) como marcas superiores.
"""
function plot_actogram end

"""
    plot_profiles_compare(perfiles; smooth = 15) -> Figure

Perfiles medios de 24 h de varios registros superpuestos; `perfiles` es un vector de
`(etiqueta, perfil, L5_medio_h)`, con `perfil` de `Circadian.profile24h`. El punto medio
de L5 (el ancla) se marca con una línea del mismo color.
"""
function plot_profiles_compare end

"""
    plot_enmo_ecdf(series; umbrales = (63.3, 142.6, 464.6)) -> Figure

Distribución acumulada de ENMO (mg, escala logarítmica) de varios registros (`series`:
vector de `(etiqueta, valores)`), con los umbrales de intensidad de Hildebrand en cadera.
"""
function plot_enmo_ecdf end

"""
    plot_timing_compare(datos) -> Figure

Diana conjunta del timing: `datos` es un vector de `(etiqueta, fase, dosis, θ, R)` (fase y
θ en radianes desde el ancla). Cada registro ocupa un anillo y tiene su vector medio.
"""
function plot_timing_compare end

"""
    plot_sleep_raster(datos; ggir = false) -> Figure

Periodos de reposo nocturno por noche de varios registros (`datos`: vector de `(etiqueta, noches)` con
`inicio_h`, `fin_h`, `libre`, `metodo` y, opcional, `L5_noche_h`). Diario en azul y noches
libres en ámbar; SPT* (retirada) en trazo claro. Con `ggir = true`, `noches` trae también
`inicio_ggir` y `fin_ggir` y se dibujan debajo, en gris.
"""
function plot_sleep_raster end

"""
    plot_vigorous_summary(r, ventana, datos; offset, scale) -> Figure

Arriba, ‖a‖ a la frecuencia de muestreo alrededor de un episodio vigoroso (`ventana =
(t0, t1, ep_inicio, ep_fin)`); abajo, para cada registro de `datos` (vector de
`(etiqueta, episodios, dias)`), la duración de los episodios y su dosis por hora del día.
"""
function plot_vigorous_summary end

"""
    plot_circular_plane(A, θ0_h, puntos) -> Figure

Plano (R cos θ, R sin θ) del modelo metabólico: de fondo, el HOMA-IR relativo predicho
exp(A·R·cos(θ − θ0)) de una cohorte sintética; encima, los `puntos` (vector de
`(etiqueta, R, θ_h)`). 0 h arriba y sentido horario, como la diana.
"""
function plot_circular_plane end

"""
    plot_bland_altman(datos; unidad = "mg") -> Figure

Bland–Altman por registro (`datos`: vector de `(etiqueta, a, b)`): media frente a diferencia
a − b como densidad hexagonal, con el sesgo y los límites de acuerdo (±1,96 DE).
"""
function plot_bland_altman end

"""
    plot_imputation_dumbbell(filas) -> Figure

Efecto de imputar el *idle sleep* (`filas`: `id`, `variante`, `ancla_h`, `θ_h`).
Cada registro es una fila; las tres variantes (huecos, ceros, última muestra) se unen
por el arco corto, para que un salto de medianoche no se dibuje como 20 h.
"""
function plot_imputation_dumbbell end

"""
    plot_intensity_ggir(filas) -> Figure

Minutos diarios por intensidad, Brincadiano frente a GGIR (`filas`: `id`, `banda`,
`nuestro`, `ggir`). Un panel por registro; cada banda es un dumbbell.
"""
function plot_intensity_ggir end

"""
    plot_agreement_timeline(datos; paso = Minute(15)) -> Figure

Acuerdo en el tiempo de dos clasificaciones binarias (p. ej., no uso nuestro y de GGIR):
`datos` es un vector de `(etiqueta, inicio, a, b)` con las ventanas de duración `paso`; cada
registro es un par de franjas a lo largo de los días y las discrepancias se marcan en rojo.
"""
function plot_agreement_timeline end

"""
    plot_timing_sensitivity(filas; ancho = :texto) -> Figure

Sensibilidad del *timing* (*forest plot*): `filas` es una tabla con `variante`, `id`,
`θ_h`, `θ_ic_desde_h`, `θ_ic_hasta_h`, `R`, `R_ic_inf` y `R_ic_sup`, con las variantes en
el orden en que se dibujan (de arriba abajo; la primera es la de referencia). A la
izquierda, θ (horas tras el ancla) con su intervalo *bootstrap*; a la derecha, R. Cada
registro tiene su color y una línea vertical en su θ de referencia; si hay columna `grupo`,
una línea fina separa los grupos de variantes.
"""
function plot_timing_sensitivity end

"""
    plot_rhythm_spectrum(datos; fmax = 8) -> Figure

Espectro del ritmo de varios registros (`datos`: vector de `(etiqueta, gls, ventana, chi2)`,
con `gls` de `Spectral.gls_periodogram`, `ventana` de `Spectral.spectral_window` en las
mismas frecuencias y `chi2` de `Spectral.chi2_periodogram`). Tres paneles: varianza
explicada por frecuencia (líneas finas en los armónicos de 24 h), ventana espectral del
muestreo y Qp/Qp_crit por periodo.
"""
function plot_rhythm_spectrum end

"""
    plot_harmonic_profile(datos; K = (1, 4)) -> Figure

Perfil horario medio (puntos) y su reconstrucción con los primeros armónicos de 24 h, un
panel por registro (`datos`: vector de `(etiqueta, h)` con `h` de `Spectral.harmonic_is`).
La leyenda da la fracción de IS que recoge cada reconstrucción.
"""
function plot_harmonic_profile end

"""
    plot_psd_classes(datos, frecuencias; suelo = nothing, cuantizacion = nothing) -> Figure

PSD media del movimiento por clase de intensidad, un panel por registro (`datos`: vector de
`(etiqueta, clases)`, con `clases` un vector de `(nombre, psd_media)`), en escala
logarítmica. `suelo` (vector de PSD, p. ej. la de las ventanas más quietas) y
`cuantizacion` (nivel en m𝑔²/Hz) se dibujan como referencias.
"""
function plot_psd_classes end

"""
    plot_cadence(datos, frecuencias) -> Figure

Distribución de la frecuencia dominante por clase de intensidad, un panel por registro
(`datos`: vector de `(etiqueta, clases)`, con `clases` un vector de `(nombre, f_dom)`).
"""
function plot_cadence end

"""
    plot_episode_spectral_qc(datos; umbral = 0.5) -> Figure

QC espectral de la actividad vigorosa, un panel por registro (`datos`: vector de
`(etiqueta, episodios, aislados)`, cada uno un `DataFrame` con `f_dom` y `frac_3`):
frecuencia dominante frente a fracción de potencia en 10–15 Hz de la ventana de cada
episodio y de cada segundo vigoroso aislado; la línea marca el umbral de vibración.
"""
function plot_episode_spectral_qc end

"""
    plot_qc_days(datos; minimo_h = 10) -> Figure

Composición de las horas de cada día natural, un panel por registro (`datos`: vector de
`(etiqueta, dias)`, con `dias` como la tabla de `QC.wear_by_day`: `dia`, `cubiertas_h`,
`no_uso_h`, `inmovil_contacto_h`, `sin_datos_h`, `uso_h`, `valido`). Cada día es una barra
de 24 h: uso (con lo inmóvil con contacto aparte, que cuenta como uso), no uso probable,
sin datos y fuera del registro; la línea marca el mínimo de uso de un día válido.
"""
function plot_qc_days end

"""
    plot_capsense_runs(datos; ratio = 0.25) -> Figure

Tramos inmóviles del control de calidad, un panel por registro (`datos`: vector de
`(etiqueta, tramos, contacto_ref)`, con `tramos` como los de `QC.qc_report` —`horas`,
`contacto` y `nocturno` en fracción, `etiqueta`— y `contacto_ref` en fracción): duración
frente a fracción de contacto con la piel (CAPSENSE), con la regla de `QC.label_runs`
(`ratio`·`contacto_ref`) y el contacto de referencia del propio registro.
"""
function plot_capsense_runs end

"""
    plot_anchor_nights(datos) -> Figure

Coherencia del ancla noche a noche (`datos`: vector de `(etiqueta, noches)`, con `noches`
como la tabla de `Sleep.sleep_nights` del script 08: `dia_semana`, `metodo`,
`spt_menos_L5_h` y `pct_spt_idle`). Izquierda: punto medio del periodo de reposo nocturno
menos punto medio de L5 de esa noche (h); derecha: porcentaje del periodo en *idle sleep*
(ángulo heredado).
"""
function plot_anchor_nights end

"""
    plot_timing_days(datos; ancla = "punto medio de L5") -> Figure

Estabilidad del *timing* entre días, una diana por registro (`datos`: vector de
`(etiqueta, por_dia, θ_h, R, ic)`, con `por_dia` como la tabla de
`Circadian.timing_by_day` —`dia`, `θ_h`, `R`—, θ y R de la semana y `ic` el intervalo
*bootstrap* de θ como `(desde_h, hasta_h)`). Cada día es un vector (diario o fin de
semana); el de la semana y su intervalo, encima. 0 h arriba (el ancla), sentido horario.
"""
function plot_timing_days end

"""
    plot_clock_coverage(datos) -> Figure

Cobertura del reloj y su corrección, una columna por registro (`datos`: vector de
`(etiqueta, cobertura, episodios, dias_validos, ancla_h)`, con `cobertura` el vector de
86 400 s de `Circadian.clock_coverage` y `episodios` con `medio`, `dosis` y `peso`).
Arriba: días observados por hora del reloj. Abajo: dosis vigorosa por hora del reloj, por
día válido (sin corregir) y por día observado a esa hora (corregida, `coverage_weights`).
"""
function plot_clock_coverage end

"""
    plot_threshold_sweep(datos; referencia = 464.6, marcas = ()) -> Figure

*Timing* según el umbral de actividad vigorosa (`datos`: vector de `(etiqueta, barrido)`,
con `barrido` una tabla de `umbral_mg`, `episodios`, `θ_h` y `R`). De arriba abajo:
desplazamiento de θ respecto al del umbral de `referencia` (h, en [−12, 12)), R y número de
episodios; líneas verticales en la referencia y en `marcas`.
"""
function plot_threshold_sweep end

"""
    plot_nights_subsample(datos; variables) -> Figure

Error de estimar el reposo nocturno con pocas noches (`datos`: vector de `(etiqueta, resumen)`, con
`resumen` la tabla de `Sleep.night_subsample`: `k`, `variable`, `mediana_abs`, `max_abs`). Un panel por
variable de `variables` (pares `(columna, rótulo)`; por defecto el punto medio y la duración del
periodo): error absoluto de la media de `k` noches frente a la de todas, en horas; línea continua, la
mediana de todos los subconjuntos de `k` noches; discontinua, el peor.
"""
function plot_nights_subsample end

"""
    plot_light_anchored(datos) -> Figure

Luz anclada al sueño, día a día (`datos`: vector de `(etiqueta, dias)`, con `dias` como la
tabla de `Light.light_daily`). Izquierda: del despertar al inicio del reposo nocturno
(`despertar_h`, `dormir_h`), con la primera, la última y la hora media de luz alta
(`primera_alta_h`, `ultima_alta_h`, `hora_media_alta_h`) y las ventanas de 2 h tras
despertar y 3 h antes de dormir sombreadas. Derecha: % de luz válida en esas ventanas.
"""
function plot_light_anchored end

"""
    plot_mvpa_light(datos) -> Figure

Minutos de AFMV de cada día según la luz (`datos`: vector de `(etiqueta, dias)` de
`Light.light_daily`): a oscuras, interior, compatible con exterior y sin luz válida,
apilados; laborables y fin de semana se distinguen en el rótulo.
"""
function plot_mvpa_light end

"""
    plot_light_enmo_density(datos; dark = 1, active_mg = 63.3) -> Figure

Densidad conjunta de luz y movimiento por minuto, un panel por registro (`datos`: vector
de `(etiqueta, d)`, con `d` la serie de `Light.light_states`): log₁₀(1 + luz media bruta)
frente a log₁₀(1 + ENMO), solo minutos despierto con el dispositivo puesto y luz conocida.
Las líneas son los umbrales de `light_states`: la esquina inferior derecha es «oscuro y
activo» (oclusión o interior por debajo del umbral de detección).
"""
function plot_light_enmo_density end

"""
    plot_light_states_by_hour(datos) -> Figure

Estado de la luz según la hora del reloj, un panel por registro (`datos`: vector de
`(etiqueta, d)` de `Light.light_states`): proporción de minutos de cada hora en cada uno
de `Light.ESTADOS_LUZ`, apilada.
"""
function plot_light_states_by_hour end

"""
    plot_light_profile(datos; min_validos = 10) -> Figure

Perfil de 24 h de la luz válida y del movimiento, un panel por registro (`datos`: vector de
`(etiqueta, d)` de `Light.light_states`): mediana y rango intercuartílico de la luz media
bruta de los minutos válidos de cada hora (solo horas con al menos `min_validos`), en
escala logarítmica, y ENMO medio de los minutos despierto con el dispositivo puesto.
"""
function plot_light_profile end

end # module
