# Micro-CT and Molecular Data for Doridida

Data and code associated with **Micro-CT and molecular phylogenetics suggest evolutionary patterns in spicule and shell architecture within Doridida (Gastropoda: Heterobranchia)**.

The manuscript's statistical analysis and figure workflow uses **R only**. The 2026-10-08 release uses 51 specimens, with mantle and foot paired within specimens. See the [R workflow instructions](Doridida%20spicule%20and%20phylogeny/README.md).

## Repository contents

- Doridida spicule and phylogeny/: canonical morphometric input, R code, tests, dated R results and provenance.
- Molecular_Phylogeny/: historical IQ-TREE/MrBayes inputs, trees and summaries.
- Trace_History/: historical Mesquite and ancestral-state artifacts.
- Morphological_Analysis/: historical morphology files; its CSV matches the canonical SpA-coded input.

Raw Micro-CT volumes and specimen images are archived separately in MorphoSource. Large sampled MrBayes trees and local session histories remain excluded by existing Git ignore rules.

## Validation status

R morphometric analyses and regression checks completed on 2026-10-08. Some permutation p-values differ from the manuscript snapshot and require reconciliation. Tree parsing is narrower than inference validation: historical MCC support labels and several Mesquite mappings remain unresolved.

See the [validation report](Doridida%20spicule%20and%20phylogeny/metadata/VALIDATION_20261008.md) and [provenance](Doridida%20spicule%20and%20phylogeny/metadata/ANALYSIS_PROVENANCE.md).

Earlier auxiliary scripts/results remain as historical records. The historical Python script and its outputs are excluded from the manuscript's current workflow and validation evidence.
