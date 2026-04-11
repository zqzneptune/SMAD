# Statistical Modelling of AP-MS Data (SMAD)

[![Bioc Release](https://www.bioconductor.org/shields/build/release/bioc/SMAD.svg)](https://www.bioconductor.org/checkResults/release/bioc-LATEST/SMAD/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `SMAD` package implements statistical modelling of affinity purification–mass spectrometry (AP-MS) data to compute confidence scores for identifying *bona fide* protein-protein interactions (PPI). By assigning probability scores, `SMAD` facilitates the removal of non-specific background contaminants commonly found in proteomics data.

## Key Features

- **Multiple Scoring Models**: Support for validated algorithms including CompPASS, HGScore, SAINTexpress, PE, DICE, and Hart.
- **Spoke & Matrix Models**: Flexibility to work with different topological interpretations of AP-MS data.
- **Bioconductor Integration**: Built to standard proteomics community requirements.
- **User-Friendly API**: Consistent interface across different scoring methods.

## Installation

### From Bioconductor (Recommended)

To install the stable version from Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SMAD")
```

### From GitHub (Development Version)

To install the latest development version:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("zqzneptune/SMAD")
```

## Input Data Format

`SMAD` requires input data as a standard `data.frame`. A built-in dataset `TestDatInput` is provided to demonstrate the required structure:

```r
library(SMAD)
data(TestDatInput)
head(TestDatInput)
```

| Column | Description | Required For |
|:---:|:---|:---|
| **idRun** | Unique identifier for each AP-MS run | All methods |
| **idBait** | Identifier for the bait protein | CompPASS, PE, SAINT |
| **idPrey** | Identifier for the identified prey protein | All methods |
| **countPrey** | Prey quantitative measure (e.g., spectral counts) | CompPASS, HG, SAINT |
| **lenPrey** | Protein sequence length of the prey | HG, SAINT |

*Tip: For replicates, ensure **idRun** is unique (e.g., by adding a suffix like "_A", "_B") to uniquely identify each experiment.*

## Quick Usage

### 1. CompPASS
The Comparative Proteomic Analysis Software Suite (CompPASS) identifies high-confidence interactors by comparing occurrences across multiple experiments.

```r
# Returns scores including WD-score, Z-score, S-score, and D-score
datScore <- CompPASS(TestDatInput)
```

### 2. HGScore
Based on a hypergeometric distribution error model incorporating NSAF for length normalization.

```r
datScore <- HG(TestDatInput)
```

### 3. SAINTexpress
An integrated version of the widely used SAINT algorithm for significance analysis of interactomes.

```r
# Supporting both spectral counts (spc) and intensities (int)
# See ?SAINTexpress_spc for detailed parameter documentation
```

## Documentation

For detailed information on the scoring algorithms and more examples, please refer to the package vignettes:

- **Introduction to SMAD**: `vignette("quickstart", package = "SMAD")`
- **Detailed Scoring Functions**: `vignette("scoring_functions", package = "SMAD")`

## License

`SMAD` is released under the **MIT License**.

© 2018-2024 Qingzhou Zhang

---

### References

- **CompPASS**: Sowa et al. (2009) [Cell 138(2):389-403](https://doi.org/10.1016/j.cell.2009.04.042)
- **HGScore**: Hart et al. (2007) [BMC Bioinformatics 8:236](https://doi.org/10.1186/1471-2105-8-236)
- **DICE**: Zhang et al. (2008) [Bioinformatics 24(7):979–986](https://doi.org/10.1093/bioinformatics/btn036)
- **PE**: Collins et al. (2007) [Mol Syst Biol 3:88](https://doi.org/10.1038/msb4100128)
- **SAINTexpress**: Teo et al. (2014) [J Proteomics 100:37-43](https://doi.org/10.1016/j.jprot.2013.11.004)