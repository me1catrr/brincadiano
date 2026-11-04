# Fixtures de test

Solo se versionan aquí archivos **públicos que no proceden de ningún participante**.
Los datos del estudio viven en `data/raw/` y nunca se versionan.

| Archivo | Procedencia | Licencia | Verificación |
|---|---|---|---|
| `legacy_sample.gt3x` | `src/GT3X.Parsing.Examples/sample.gt3x` de [actigraph/NHANES-GT3X-File-Format](https://github.com/actigraph/NHANES-GT3X-File-Format), descargado el 2026-10-02 | MIT (ActiGraph) | blob git `48cce9c7a4aa50e4b2c5b5991a098e7b82f8569c` (igual al publicado por la API de GitHub) |

`legacy_sample.gt3x` es un registro de demostración del fabricante en el formato
antiguo (GT3X+ `NEO`, firmware 1.4.0): `activity.bin` + `lux.bin`, sin `log.bin`.
