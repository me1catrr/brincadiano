"""
    Brincadiano

Análisis de acelerometría RAW ActiGraph (GT3X+/wGT3X-BT) en niños: actividad física,
sueño, ritmos circadianos y luz, y su relación con la salud metabólica infantil.

Submódulos:
- [`Rutas`](@ref): rutas centrales y estructura de `data/raw/` por proyecto.
- [`GT3X`](@ref): lectura nativa de archivos `.gt3x`.
- [`QC`](@ref): control de calidad por archivo.
- [`Metrics`](@ref): ENMO, MAD y agregación en epochs.
- [`Circadian`](@ref): L5/M10, IS/IV y estadística circular.
- [`Sleep`](@ref): periodo de reposo nocturno sin diario (postura en cadera) e inactividad sostenida.
- [`Model`](@ref): modelo metabólico, cohortes sintéticas y potencia.
- [`Spectral`](@ref): análisis espectral del ritmo (Lomb–Scargle, χ², armónicos de 24 h).
- [`Light`](@ref): procesamiento de la señal de luz (lux).
- [`AGD`](@ref): lectura de los `.agd` de ActiLife (solo como referencia, lista blanca).
- [`Viz`](@ref): gráficas (implementadas en la extensión de CairoMakie).
"""
module Brincadiano

include("Rutas.jl")
include("GT3X.jl")
include("AGD.jl")
include("QC.jl")
include("Metrics.jl")
include("Circadian.jl")
include("Sleep.jl")
include("Model.jl")
include("Spectral.jl")
include("Light.jl")
include("Viz.jl")

export Rutas, GT3X, AGD, QC, Metrics, Circadian, Sleep, Model, Spectral, Light, Viz

end # module
