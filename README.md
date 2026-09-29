# heart-capillary-sex-differences

Code for: **Sex differences in human heart capillaries: a multi-cohort analysis of capillary endothelium and pericytes** (Adelus, in preparation).

Coronary microvascular dysfunction and heart failure with preserved ejection fraction are more common in women. This project asks whether healthy women's and men's heart capillaries differ in the genes that carry out capillary function. It pools every healthy adult heart ventricle in the CELLxGENE Census and analyzes each donor as one unit. Nulls are tested with equivalence bounds, so a null result is informative rather than just a failure to find a difference.

The pre-specified plan is in [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md) and was committed before any expression data were analyzed.

## Steps

| Script | What it does |
|---|---|
| `R/00_setup.R` | Packages, Census connection, pre-specified settings and gene sets |
| `R/01_inventory.R` | Selects eligible donors and datasets (metadata only) |
| `R/02_download.R` | Downloads ventricular cells one donor at a time and keeps pseudobulks |
| `R/03_ec_subtypes_pseudobulk.R` | Assigns capillary, arterial and venous ECs with one rule across datasets |
| `R/04_aim1_capillary_programs.R` | Primary analysis: six capillary-function programs, meta-analysis, equivalence, calibration |
| `R/05_aim2_structure.R` | Pericyte-to-capillary ratio, composition, capillary–pericyte signaling |
| `R/06_aim3_age.R` | Sex × age band (menopause proxy), descriptive |
| `R/07_aim4_sex_chromosomes.R` | X–Y paralogs across cell types in the same donors |
| `R/08_figures.R` | Figures and tables |

Run everything with `source("run_all.R")` from this folder. The first run downloads several GB of Census data and takes a few hours; downloads resume where they stopped.

## Data

All data are open: the [CZ CELLxGENE Discover Census](https://chanzuckerberg.github.io/cellxgene-census/), release 2025-11-08. No login or data-access application is needed.

## Author

Maria Adelus · [ORCID 0000-0002-9676-9214](https://orcid.org/0000-0002-9676-9214) · MIT License
