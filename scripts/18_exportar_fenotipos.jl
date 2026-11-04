# 18 · Exportación de los fenotipos por «bloques de análisis».
#
# Qué hace
#   Reúne en una tabla por bloque (una fila por registro, solo el ID pseudonimizado) las
#   variables que el pipeline ya calcula, y escribe un diccionario con su significado, su
#   unidad y su estado (CONFIRMADO / PROVISIONAL / PENDIENTE / NO VÁLIDO). No recalcula nada:
#   lee las tablas `*_resumen.csv` de los scripts 04, 05, 06, 08, 09 y 12. Es el formato en que
#   se entregarían las variables de una cohorte a un análisis posterior, y el punto de partida
#   para añadir bloques nuevos (por ejemplo, la actividad al aire libre por intensidad).
#
# Cómo ejecutarlo
#   julia --project=scripts scripts/18_exportar_fenotipos.jl                  # todos los bloques
#   julia --project=scripts scripts/18_exportar_fenotipos.jl calidad timing   # solo esos
#
# Bloques: calidad, volumen, ritmo, espectro, timing, sueno, luz.
#
# Depende de
#   04, 05, 06, 08, 09 y 12, que escriben las tablas de `outputs/*/tables/`. Un bloque cuya
#   tabla falte se omite y se dice; no se rellena.
#
# Salida (individual; no se versiona)
#   outputs/fenotipos/fenotipos_<bloque>.csv   una fila por registro, con `id` primero
#   outputs/fenotipos/fenotipos_ancho.csv      todos los bloques unidos por `id`, columnas `bloque__variable`
#   outputs/fenotipos/diccionario.csv          bloque, variable, descripción, unidad, estado
#
# Consola: solo el ID, el número de variables y las omisiones. Las unidades de la luz son brutas
# del sensor: la conversión a lux del wGT3X-BT no está resuelta.

using CSV, DataFrames, Printf

const RAIZ = normpath(joinpath(@__DIR__, ".."))
const ENT = joinpath(RAIZ, "outputs")
const SAL = joinpath(ENT, "fenotipos")

# bloque => (archivo de origen, [(variable, descripción, unidad, estado), ...])
const BLOQUES = Dict(
    "calidad" => ("qc/tables/qc_resumen.csv", [
        ("horas", "Duración del registro", "h", "confirmado"),
        ("pct_muestreado", "Tiempo realmente muestreado (el resto es idle sleep)", "%", "confirmado"),
        ("fs_efectiva", "Frecuencia de muestreo efectiva", "Hz", "confirmado"),
        ("idle_n", "Número de huecos de idle sleep", "n", "confirmado"),
        ("idle_h", "Tiempo total en idle sleep", "h", "confirmado"),
        ("idle_max_h", "Hueco de idle sleep más largo", "h", "confirmado"),
        ("pct_saturadas", "Muestras saturadas", "%", "confirmado"),
        ("no_uso_n", "Tramos de no uso probable", "n", "confirmado"),
        ("no_uso_h", "Tiempo de no uso probable", "h", "confirmado"),
        ("contacto_ref_pct", "Contacto CAPSENSE durante el uso", "%", "provisional"),
        ("dias", "Días del registro", "n", "confirmado"),
        ("dias_validos", "Días válidos (≥ 10 h con el dispositivo puesto)", "n", "confirmado"),
    ]),
    "volumen" => ("metrics/tables/metricas_resumen.csv", [
        ("cal_aplicada", "¿Se aplicó la autocalibración (van Hees 2014)?", "lógico", "confirmado"),
        ("sesgo_mg", "Sesgo en reposo antes de calibrar", "mg", "confirmado"),
        ("enmo_medio_mg", "ENMO medio de los días válidos", "mg", "confirmado"),
        ("mad5_mediana_mg", "Mediana de la MAD en epochs de 5 s", "mg", "confirmado"),
        ("enmo5_p99_mg", "Percentil 99 de la ENMO en epochs de 5 s", "mg", "confirmado"),
        ("dias_validos", "Días válidos", "n", "confirmado"),
    ]),
    "ritmo" => ("circadian/tables/ritmo_resumen.csv", [
        ("L5", "ENMO medio de las 5 h menos activas", "mg", "confirmado"),
        ("L5_medio_h", "Punto medio de L5", "h", "confirmado"),
        ("L5_estimable", "¿L5 estimable (cubre ≥ 90 % del día)?", "lógico", "confirmado"),
        ("M10", "ENMO medio de las 10 h más activas", "mg", "confirmado"),
        ("M10_medio_h", "Punto medio de M10", "h", "confirmado"),
        ("RA", "Amplitud relativa", "", "confirmado"),
        ("IS", "Estabilidad interdiaria", "", "confirmado"),
        ("IV", "Variabilidad intradiaria", "", "confirmado"),
        ("ancla", "Ancla del ritmo (L5, o L5* con el no uso como reposo)", "texto", "confirmado"),
    ]),
    "espectro" => ("circadian/tables/espectro_resumen.csv", [
        ("pico_periodo_h", "Periodo del pico del periodograma", "h", "provisional"),
        ("potencia_24h", "Varianza explicada por el ritmo de 24 h", "", "provisional"),
        ("amplitud_24h", "Amplitud del primer armónico", "mg", "provisional"),
        ("acrofase_24h", "Acrofase del primer armónico", "h", "provisional"),
        ("IS_armonicos", "IS por armónicos", "", "provisional"),
    ]),
    "timing" => ("timing/tables/timing_resumen.csv", [
        ("ancla", "Ancla del timing", "texto", "confirmado"),
        ("episodios", "Episodios vigorosos", "n", "confirmado"),
        ("seg_vig_dia", "Segundos vigorosos por día", "s", "confirmado"),
        ("θ_h", "Hora media de la actividad vigorosa tras el ancla (corregida por cobertura)", "h", "confirmado"),
        ("R", "Concentración de la hora (longitud resultante)", "", "confirmado"),
        ("R_cos", "R·cos θ (predictor del plano circular)", "", "confirmado"),
        ("R_sin", "R·sin θ (predictor del plano circular)", "", "confirmado"),
        ("θ_ic_amplitud_h", "Amplitud del intervalo de θ (bootstrap de días)", "h", "provisional"),
        ("θ_informable", "¿R permite interpretar θ?", "lógico", "provisional"),
    ]),
    "sueno" => ("sleep/tables/sueno_resumen.csv", [
        ("noches", "Noches analizadas", "n", "provisional"),
        ("postura", "Noches con periodo por la postura", "n", "provisional"),
        ("retirada", "Noches con retirada (SPT*)", "n", "provisional"),
        ("inicio_h", "Hora media de inicio del reposo nocturno estimado", "h", "provisional"),
        ("fin_h", "Hora media de fin del reposo nocturno estimado", "h", "provisional"),
        ("punto_medio_h", "Punto medio del reposo nocturno estimado", "h", "provisional"),
        ("duracion_h", "Duración del periodo", "h", "provisional"),
        ("sueno_h", "Sueño dentro del periodo", "h", "provisional"),
        ("eficiencia", "Eficiencia (exploratoria en cadera)", "", "provisional"),
        ("jetlag_social_h", "Jet lag social", "h", "provisional"),
    ]),
    "luz" => ("light/tables/luz_resumen.csv", [
        ("dias", "Días con ≥ 10 h despierto con el dispositivo", "n", "provisional"),
        ("pct_valido", "Tiempo despierto con luz válida (distinta de cero)", "%", "confirmado"),
        ("pct_oscuro_activo", "Tiempo despierto a oscuras y en movimiento", "%", "confirmado"),
        ("pct_oscuro_inactivo", "Tiempo despierto a oscuras y quieto", "%", "confirmado"),
        ("pct_sin_dato", "Tiempo despierto sin luz conocida", "%", "confirmado"),
        ("min_sobre_bajo", "Minutos con luz ≥ 100 (bruto)", "min", "provisional"),
        ("min_sobre_alto", "Minutos con luz ≥ 1000 (bruto)", "min", "provisional"),
        ("hora_media_alta_h", "Hora media (circular) de luz alta", "h", "provisional"),
        ("afmv_exterior", "Minutos de AFMV con luz compatible con exterior", "min", "provisional"),
        ("antes_dormir_pct_valido", "Luz válida en las 3 h antes de dormir", "%", "novalido"),
    ]),
)
const ORDEN = ["calidad", "volumen", "ritmo", "espectro", "timing", "sueno", "luz"]

pedidos = isempty(ARGS) ? ORDEN : ARGS
for b in pedidos
    haskey(BLOQUES, b) || error("bloque desconocido: $b (válidos: $(join(ORDEN, ", ")))")
end
mkpath(SAL)

dicc = DataFrame(bloque = String[], variable = String[], descripcion = String[], unidad = String[], estado = String[])
ancho = nothing
for b in pedidos
    archivo, vars = BLOQUES[b]
    ruta = joinpath(ENT, archivo)
    if !isfile(ruta)
        @printf("omitido  %-9s no existe %s (¿se ejecutó su script?)\n", b, archivo)
        continue
    end
    t = CSV.read(ruta, DataFrame)
    "id" in names(t) || error("$archivo no tiene columna id")
    presentes = [v for v in vars if v[1] in names(t)]
    faltan = [v[1] for v in vars if !(v[1] in names(t))]
    isempty(faltan) || @printf("aviso    %-9s faltan en la tabla: %s\n", b, join(faltan, ", "))
    sal = select(t, "id", [v[1] for v in presentes]...)
    CSV.write(joinpath(SAL, "fenotipos_$(b).csv"), sal)
    for v in presentes
        push!(dicc, (b, v[1], v[2], v[3], v[4]))
    end
    ren = rename(sal, [n => (n == "id" ? "id" : "$(b)__$(n)") for n in names(sal)])
    global ancho = ancho === nothing ? ren : outerjoin(ancho, ren, on = "id")
    @printf("exportado %-9s %d registros, %d variables\n", b, nrow(sal), length(presentes))
end

if ancho !== nothing
    CSV.write(joinpath(SAL, "fenotipos_ancho.csv"), sort(ancho, "id"))
    CSV.write(joinpath(SAL, "diccionario.csv"), dicc)
    @printf("ancho    %d registros, %d columnas; diccionario con %d variables\n", nrow(ancho), ncol(ancho), nrow(dicc))
end
