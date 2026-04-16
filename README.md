# SMAD: Statistical Modelling of AP-MS and BioID Data

[![Bioc Release](https://www.bioconductor.org/shields/build/release/bioc/SMAD.svg)](https://www.bioconductor.org/checkResults/release/bioc-LATEST/SMAD/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`SMAD` is a high-performance Bioconductor R package designed for the statistical analysis of proteomics data, specifically focusing on **Affinity Purification Mass Spectrometry (AP-MS)** and **BioID / Proximity-dependent Labeling** (e.g., TurboID, APEX2). 

The package implements a comprehensive suite of validated algorithms to assign confidence scores to identified **Protein-Protein Interactions (PPI)**, enabling robust background contaminant removal and the identification of *bona fide* interactors from high-throughput mass spectrometry experiments.

## Why use SMAD?

- **Optimized for BioID & AP-MS**: Tailored workflows for spectral count and intensity-based proteomics data.
- **Robust Background Removal**: Efficiently filters non-specific contaminants common in affinity purification.
- **Validated Algorithms**: Access to world-class scoring models like SAINTexpress, CompPASS, and HGScore.
- **Bioconductor Standards**: Fully integrated with the Bioconductor ecosystem for seamless proteomics pipelines.
- **Speed & Efficiency**: Core algorithms implemented in C++ for fast processing of large datasets.

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

## BioID & Proximity Labeling Workflow

For BioID, TurboID, or APEX2 experiments, `SAINTexpress` is the recommended scoring method. Below is a typical workflow using spectral counts:

```r
library(SMAD)

# 1. Prepare your data (Interactome, Prey, and Bait tables)
# 2. Run SAINTexpress for significance analysis
# results <- SAINTexpress_spc(your_inter_df, your_prey_df, your_bait_df)

# The output includes:
# - SaintScore: Probability of interaction
# - BFDR: Bayesian False Discovery Rate
# - FoldChange: Enrichment over controls
```

## Core Scoring Algorithms

`SMAD` provides a unified interface to the most widely used algorithms in the proteomics community:

### 1. SAINTexpress (The Gold Standard for BioID)
Integrated R version of the [SAINTexpress](https://doi.org/10.1016/j.jprot.2013.11.004) algorithm. Ideal for **BioID**, **TurboID**, and standard **AP-MS** data with control experiments. It supports both spectral counts (spc) and intensities (int).

```r
# Ideal for BioID and Proximity Labeling data
# result <- SAINTexpress_spc(inter, prey, bait)
```

### 2. CompPASS
The **Comparative Proteomic Analysis Software Suite**. Efficiently identifies high-confidence interactors by comparing bait-prey occurrences across multiple diverse experiments without requiring explicit controls.

```r
# Returns scores including WD-score, Z-score, S-score, and D-score
datScore <- CompPASS(TestDatInput)
```

### 3. HGScore (Hypergeometric Scoring)
A distribution-based model that incorporates NSAF (Normalized Spectral Abundance Factor) for length normalization, suitable for high-throughput PPI mapping.

```r
datScore <- HG(TestDatInput)
```

### 4. Additional Models
- **DICE**: Specifically designed for quantifying interaction significance in large networks.
- **PE (Prey Equivalence)**: Scoring based on the probability of observing a prey across baits.
- **Hart**: Likelihood-based scoring for topological analysis.

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