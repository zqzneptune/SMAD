# SMAD: Statistical Modelling of AP-MS Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SMAD** is an R package designed for the statistical analysis of Affinity Purification–Mass Spectrometry (AP-MS) data. It implements several widely-used algorithms to compute confidence scores, helping researchers distinguish *bona fide* protein-protein interactions (PPI) from non-specific background noise.

## Table of Contents
- [Installation](#installation)
- [Input Data Format](#input-data-format)
- [Available Scoring Algorithms](#available-scoring-algorithms)
- [Quick Start](#quick-start)
- [References](#references)

## Installation

You can install the development version of **SMAD** from GitHub using `devtools` or `remotes`:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install SMAD
devtools::install_github("zqzneptune/SMAD")

# Load the package
library(SMAD)
```

## Input Data Format

To ensure compatibility with SMAD scoring functions, your input `data.frame` should follow a specific structure. You can explore the built-in demo dataset:

```r
data(TestDataInput)
head(TestDataInput)
```

### Required Columns

| Column | Description |
|:-------|:------------|
| **idRun** | Unique identifier for a specific AP-MS run (replicate). |
| **idBait** | Unique identifier for the bait protein. |
| **idPrey** | Unique identifier for the prey protein. |
| **countPrey** | Spectral counts (or peptide counts) for the prey protein. |
| **lenPrey** | Protein sequence length of the prey (required for HGScore/NSAF). |

> **Note on Replicates:** If you have biological or technical duplicates, ensure the `idRun` is unique (e.g., append suffixes like `_A`, `_B`). The combination of `idRun` and `idBait` should uniquely identify a single purification.

---

## Available Scoring Algorithms

### 1. CompPASS
The **Comparative Proteomic Analysis Software Suite (CompPASS)** uses a "spoke model" to identify high-confidence interactions. It outputs several metrics, including **Z-score, S-score, D-score, and WD-score**. This implementation is optimized for performance based on the BioPlex pipeline.

*   **Required columns:** `idRun`, `idBait`, `idPrey`, `countPrey`.
*   **Usage:** `datScore <- CompPASS(datInput)`

### 2. DICE
The **Dice coefficient** measures the similarity between prey pair-wise combinations. It is particularly useful for identifying preys that frequently co-occur across different runs.

*   **Required columns:** `idRun`, `idPrey`.
*   **Usage:** `datScore <- DICE(datInput)`

### 3. Hart
Based on a **hypergeometric distribution error model**, this algorithm calculates the probability of finding a prey across different bait purifications by chance.

*   **Required columns:** `idRun`, `idPrey`.
*   **Usage:** `datScore <- Hart(datInput)`

### 4. HGScore
**HGScore** enhances the Hart hypergeometric model by incorporating **NSAF** (Normalized Spectral Abundance Factor), which accounts for protein length and abundance. It is designed based on a "matrix model."

*   **Required columns:** `idRun`, `idPrey`, `countPrey`, `lenPrey`.
*   **Usage:** `datScore <- HG(datInput)`

### 5. PE (Purification Enrichment)
The **PE score** incorporates both spoke and matrix models to calculate an interaction score based on the frequency and exclusivity of prey identifications.

*   **Required columns:** `idRun`, `idBait`, `idPrey`.
*   **Usage:** `datScore <- PE(datInput)`

---

## Quick Start

```r
library(SMAD)

# Load example data
data(TestDataInput)

# Run CompPASS scoring
results <- CompPASS(TestDataInput)

# View top-scoring interactions based on WD-score
head(results[order(results$WD, decreasing = TRUE), ])
```

## License
Released under the [MIT License](https://opensource.org/licenses/MIT).  
© [Qingzhou Zhang](https://github.com/zqzneptune)

## References
1.  **CompPASS:** [Sowa et al., 2009 (Cell)](https://doi.org/10.1016/j.cell.2009.04.042)
2.  **BioPlex:** [Huttlin et al., 2015 (Cell)](https://doi.org/10.1016/j.cell.2015.06.043)
3.  **Hart Scoring:** [Hart et al., 2007 (BMC Bioinformatics)](https://doi.org/10.1186/1471-2105-8-236)
4.  **HGScore:** [Guruharsha et al., 2011 (Cell)](https://doi.org/10.1016/j.cell.2011.08.047)
5.  **DICE:** [Zhang et al., 2008 (Bioinformatics)](https://doi.org/10.1093/bioinformatics/btn036)
6.  **PE Score:** [Collins et al., 2007 (MCP)](https://doi.org/10.1074/mcp.M600381-MCP200)