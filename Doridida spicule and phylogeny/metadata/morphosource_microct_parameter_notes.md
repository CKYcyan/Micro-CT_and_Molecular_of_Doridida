# MorphoSource microCT parameter extraction

Source folder: `1.Doridida spicule and phylogeny\MicroCT_data\voi`

## Output files

- morphosource_microct_parameters_by_sample.csv: one row per sample/log, with columns named to match the requested MorphoSource fields.
- morphosource_microct_parameter_notes.md: extraction notes and field mapping.

## Event date note

The logs indicate multiple scan dates. The CSV uses the per-sample `Study Date and Time` value from the log as `Event date`.

| Event date | Count |
|---|---:|
|  | 1 |
| 2022/02/09 | 3 |
| 2022/06/13 | 3 |
| 2022/06/14 | 3 |
| 2022/06/15 | 2 |
| 2022/06/20 | 3 |
| 2022/11/23 | 1 |
| 2022/11/24 | 8 |
| 2022/11/25 | 4 |
| 2022/11/28 | 2 |
| 2023/03/23 | 2 |
| 2023/03/24 | 6 |
| 2023/03/25 | 2 |
| 2023/03/27 | 2 |
| 2023/03/28 | 3 |
| 2023/06/28 | 5 |
| 2023/06/29 | 4 |

## Field mapping

| MorphoSource field | Source in log / handling |
|---|---|
| Scanner Modality | Fixed value: X-Ray Computed Tomography (CT/microCT) |
| Creator | Fixed value: Kuan Yu Cho |
| Event date | Parsed from `Study Date and Time` |
| Software | `Reconstruction Program` |
| Filter | `Filter` |
| Exposure time | `Exposure (ms)` |
| Flux normalization | `Reference Intensity` plus `Flat Field Correction`; this suggests the relevant normalization metadata, but MorphoSource wording may require review |
| Pixel spacing calibration | `Image Pixel Size (um)` plus `Scaled Image Pixel Size (um)` |
| Shading correction | `Flat Field Correction`; when present, `Heel Effect Correction On` is appended |
| Frame averaging | `Frame Averaging` |
| Projections | `Number Of Files` |
| Voltage | `Source Voltage (kV)` |
| Power | Calculated as kV x mA from voltage and amperage |
| Amperage | `Source Current (uA)` |
| X-ray tube type | `Source Type` plus `Source spot size` |
| Detector type | `Camera Type` |
| Detector pixels X/Y | `Number Of Columns` and `Number Of Rows` |
| Detector pixels size X/Y | `Camera Pixel Size (um)` |
| Detector configuration | `Camera binning`, `Image Format`, `Depth (bits)`, rows, and columns |
| Source object distance | `Object to Source (mm)` |
| Source detector distance | `Camera to Source (mm)` |
| Rotation number | `Use 360 Rotation`, `Rotation Step (deg)`, and calculated angular coverage |
| Optical magnification | Calculated as camera pixel size / image pixel size, which indicates geometric magnification |
| Acquisition type | `Type of Detector Motion` plus `Scanning Trajectory` |

## Fields not recorded in logs

- Surrounding material
- Target type
- Target material
- Phase contrast

These fields likely need specimen handling notes, scanner documentation, or manual confirmation before submission.
