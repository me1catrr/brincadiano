# Brincadiano

Análisis en Julia de acelerometría RAW ActiGraph (GT3X+ / wGT3X-BT) en niños: actividad
física, sueño, ritmos circadianos y luz.
Colaboración con el grupo de pediatría del IMIBIC/UCO.

Paquete `Brincadiano` 0.1.0. Requiere Julia ≥ 1.11; el `Manifest.toml` de la raíz está
resuelto con Julia 1.12.3.

Este repositorio contiene solo el código de procesado de los archivos RAW: el paquete, los
scripts del pipeline y de validación, y los tests. Ningún dato original ni tratado.

## Datos

`data/` y `outputs/` contienen información de salud de menores, pseudonimizados, y **no se
versionan**. `data/raw/` es de solo lectura. Los intermedios van a `data/interim/<ID>/`
(Arrow) y los resultados a `outputs/`. En consola, logs, figuras y tablas solo aparece el
ID del archivo (por ejemplo `6005_V00`). De `info.txt` solo se leen los campos de una lista
blanca (`src/GT3X/info.jl`) y el registro METADATA de `log.bin`, con los datos del sujeto,
nunca se decodifica (`src/GT3X/logbin.jl`).

La única excepción versionada es el fixture público de ActiGraph
`test/fixtures/legacy_sample.gt3x` (licencia MIT).

## Entornos de Julia

Un workspace (`[workspace]` en el `Project.toml` de la raíz) y un solo `Manifest.toml`.
`Pkg.instantiate()` en la raíz resuelve también `scripts/`.

| Entorno | Proyecto | Para qué |
|---------|----------|----------|
| Paquete | `--project=.` | núcleo, tests |
| Análisis | `--project=scripts` | scripts (incluye CairoMakie) |

Dependencias del núcleo: Arrow, DataFrames, Dates, TimeZones, ZipArchives, StatsBase, DSP
y FFTW. CSV y Test solo entran en los tests. CairoMakie es *weak dependency*: `using
Brincadiano` no lo carga; `using Brincadiano, CairoMakie` activa la extensión
`BrincadianoCairoMakieExt`.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
BRINCADIANO_SLOW_TESTS=1 julia --project=. -e 'using Pkg; Pkg.test()'  # + 7 días frente a ActiLife
```

Los scripts se lanzan con el entorno de análisis:

```bash
julia --project=scripts scripts/01_read_gt3x.jl
```

## Módulos

| Módulo | Archivo | Qué hace |
|--------|---------|----------|
| `GT3X` | `src/GT3X.jl` | lectura nativa de `.gt3x` (`log.bin` y formato antiguo) y exportación Arrow |
| `QC` | `src/QC.jl` | muestreo, huecos, saturación, luz, no uso y días válidos |
| `Metrics` | `src/Metrics.jl` | ENMO, MAD, epochs y autocalibración |
| `Circadian` | `src/Circadian.jl` | L5/M10, IS/IV y estadística circular |
| `Sleep` | `src/Sleep.jl` | ventana de sueño en cadera e inactividad sostenida |
| `Light` | `src/Light.jl` | serie de luz, estados de calidad y métricas de exposición |
| `Model` | `src/Model.jl` | modelo metabólico, cohortes sintéticas y potencia |
| `Spectral` | `src/Spectral.jl` | espectro del ritmo (Lomb–Scargle, χ², IS por armónicos) y del movimiento (PSD por ventana) |
| `Viz` | `src/Viz.jl` | figuras; la implementación está en la extensión de CairoMakie |

Los huecos (*idle sleep*, muestras ausentes) se guardan explícitos. La imputación, cuando
se usa, es una función aparte (`GT3X.impute_idle_sleep`) y no forma parte de la lectura.

## Pipeline

Orden de los scripts de análisis. Cada uno documenta entradas y salidas en su cabecera.

| Script | Paso |
|--------|------|
| `01_read_gt3x.jl` | `.gt3x` → `data/interim/<ID>/*.arrow` |
| `02_gt3x_summary.jl` | resumen de cobertura y huecos |
| `03_explorar_sujeto.jl` | figuras de un registro (`--fixture` usa el archivo público) |
| `04_qc.jl` | control de calidad y tramos de no uso |
| `05_metricas.jl` | epochs de 1, 5 y 60 s |
| `06_circadiano.jl` | L5/M10, IS/IV y timing de la actividad vigorosa |
| `07_modelo_sintetico.jl` | potencia del modelo con cohortes sintéticas (sí se versiona) |
| `08_sueno.jl` | ventana de sueño |
| `09_luz.jl` | estados y métricas de luz |
| `10_figuras.jl` | figuras de síntesis a partir de lo ya calculado |
| `11_sensibilidad_timing.jl` | sensibilidad del timing al ancla, al umbral y al tipo de día |
| `12_espectro.jl` | espectro del ritmo y timing por armónicos |
| `13_espectro_movimiento.jl` | espectro del movimiento: PSD por intensidad, cadencia, QC espectral |
| `15_figuras_complementarias.jl` | figuras complementarias y cifras de sus pies |
| `16_lote.jl` | lectura por lotes reanudable (script 01 por archivo; `--analisis` lanza el resto) |
| `17_inventario.jl` | filas del embudo de registros y recuento por modelo y firmware |
| `18_exportar_fenotipos.jl` | tablas de fenotipos por bloque de análisis, con diccionario de variables |
| `19_interoperabilidad.jl` | comprobaciones GT3X+ ↔ wGT3X-BT y efecto de recortar a ±6 g |
| `20_sensibilidades.jl` | sensibilidad del sueño, de los episodios vigorosos y del submuestreo de días |
| `21_sueno_pocas_noches.jl` | error del reposo nocturno estimado con pocas noches |
| `22_siestas.jl` | siestas candidatas (inactividad sostenida diurna) |
| `23_agd.jl` | `.agd` de ActiLife: inventario, lectura con lista blanca y correspondencia con los `.gt3x` |
| `24_config_ggir.jl` | auditoría de los `config_*.csv` de GGIR del grupo |
| `25_inventario_raw.jl` | inventario trazable de `data/raw/` y manifiesto de migración |
| `26_comparar_ggir_cohorte_a.jl` | Brincadiano frente a la parte 5 de GGIR (COHORTE_A V03) |
| `27_migrar_raw.jl` | migración reversible de `data/raw/` a la estructura por proyecto (por defecto, en seco) |

Las tablas y figuras de un participante viven en `outputs/` y no se versionan.

## Benchmarks externos

GGIR, read.gt3x, LightLogR y pygt3x comparan resultados; no son dependencias del paquete.
Los scripts `90`–`96` lanzan R o Python como procesos y, salvo que indiquen lo contrario,
exigen `--confirmo`.

| Entorno | Dónde | Restaurar |
|---------|-------|-----------|
| R 4.6.1 + renv (GGIR 3.3.9, read.gt3x 1.2.0, LightLogR 0.10.6) | `scripts/R/` | `cd scripts/R && Rscript -e 'renv::restore()'` |
| Python 3.12 (pygt3x 0.7.1) | `.venv/` | ver `scripts/python/requirements.txt` |

Se versionan `scripts/R/renv.lock` y `scripts/python/requirements.txt`. No se versionan
`.venv/` ni la librería de renv.

| Script | Contraste |
|--------|-----------|
| `90_validate_external.jl` | lector frente a read.gt3x y pygt3x |
| `91_bibliometria_scopus.jl` | tablas del capítulo 1 a partir de dos CSV de Scopus (fuera del repo) |
| `92_validar_ggir.jl` | calibración, ENMO, no uso, L5/M10 e IS/IV frente a GGIR; media circular frente a numpy |
| `93_validar_sueno_ggir.jl` | ventana de sueño frente a GGIR |
| `94_validar_luz_lightlogr.jl` | métricas de luz frente a LightLogR |
| `95_imputacion.jl` | efecto de imputar el *idle sleep* |
| `96_validar_intensidad_sri_ggir.jl` | tiempo por intensidad e índice de regularidad del sueño frente a GGIR |

## Informe y libro

El informe LaTeX y el libro Quarto se mantienen fuera de este repositorio. La edición web
del libro está en <https://me1catrr.github.io/brincadiano/>.

## Estructura

```
src/Brincadiano.jl          módulo raíz
src/GT3X/                   lector .gt3x
ext/BrincadianoCairoMakieExt.jl
test/                       tests (Test + CSV); fixture público en test/fixtures/
scripts/                    análisis (--project=scripts)
scripts/R/                  renv de benchmark
scripts/python/             requirements de pygt3x
data/raw/                   .gt3x originales (solo lectura, no versionado)
data/interim/<ID>/          Arrow (no versionado)
outputs/                    tablas y figuras (no versionado)
```
