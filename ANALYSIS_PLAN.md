# Analysis plan: sex differences in human heart capillaries

**Author:** Maria Adelus · **Fixed before analysis of expression data.** The commit date of this file on GitHub is the time stamp. Any later change is listed under *Deviations* at the end, with the reason.

## Question

Do healthy women's and men's ventricular capillary endothelial cells differ in the gene programs that carry out capillary function? Does any difference change across a menopause proxy? Is the X–Y paralog pattern seen in endothelium specific to capillaries?

## Data (all open; CELLxGENE Census release 2025-11-08)

- **Tissue:** healthy (`disease == normal`) adult (18 years or older) heart ventricles only: left ventricle (including labels for its anterior wall and apical region), apex and interventricular septum. Right ventricle, atria, nodes, coronary arteries and unspecified "heart" or "cardiac ventricle" labels are excluded. Only primary data (`is_primary_data == TRUE`), so no cell is counted twice.
- **Eligible donors:** at least 100 endothelial cells (atlas labels) in ventricular tissue, and known sex.
- **Datasets:** a dataset is included when it has at least 3 eligible women and 3 eligible men.
- **Unit of analysis:** the donor. There is one pseudobulk per donor and cell group, summed across ventricular regions.

## Cell definitions

- **Blood endothelial cells:** every cell labeled as an endothelial cell. Lymphatic and endocardial cells are excluded.
- **Capillary, arterial and venous subtypes:** assigned with one rule in all datasets. Each cell gets a marker score (mean log(1 + CP10k)), z-scored within the dataset, and takes the subtype with the highest z.
  - Capillary markers: CA4, RGCC, BTNL9
  - Arterial markers: GJA5, HEY1, SEMA3G, GJA4, DKK2, FBLN5
  - Venous markers: ACKR1, NR2F2, CPE
  - The capillary markers exclude every outcome gene.
- **Pseudobulk minimums:** a capillary pseudobulk needs at least 30 capillary cells. Pseudobulks for other cell groups (pericytes, cardiomyocytes, fibroblasts, myeloid cells, smooth muscle) also need at least 30 cells.

## Primary outcome (Aim 1)

Six capillary-function programs, fixed here:

| Program | Genes |
|---|---|
| NO / eNOS | NOS3, KLF2, KLF4, CAV1, GCH1, SLC7A1, DDAH1, DDAH2 |
| Endothelin / ACE | EDN1, ECE1, EDNRB, ACE |
| Prostacyclin | PTGS1, PTGS2, PTGIS, PLA2G4A |
| Barrier / junction | CLDN5, CDH5, ESAM, OCLN, TJP1, JAM2, PECAM1 |
| Fatty-acid transport | CD36, FABP4, FABP5, LPL, GPIHBP1 |
| Angiogenic tip | ESM1, APLN, DLL4, ANGPT2, PGF, KCNE3 |

**Score:** the capillary pseudobulk is converted to TMM log-CPM after filterByExpr. Each gene is z-scored across donors within the dataset, and the program score is the mean z. The score is then standardized within the dataset.

**Estimand:** the women-minus-men difference in the program score, in donor-level SD units, adjusted for age (and assay when assay varies within a dataset). It is estimated per dataset by linear regression and pooled by random-effects meta-analysis (REML; fixed effect if fewer than 3 datasets).

**Inference:**

- Holm correction across the six programs.
- Equivalence (TOST): a program is called equivalent when the 90% CI lies within ±0.8 SD, the smallest effect of interest (a large effect). Whether the 90% CI also lies within the stricter ±0.5 SD is reported for each program.
- Why ±0.8: the metadata inventory (run before any expression data) found 45 eligible donors (23 women, 22 men) in 2 datasets. With this sample the minimal detectable effect is about 0.85 SD, so equivalence within ±0.5 SD is not attainable and the study is designed to rule out large differences. This choice was made from the donor count only, before any outcome was seen.
- Each program gets one of three labels:
  - *difference* (Holm p < 0.05);
  - *equivalent* (not different, and the 90% CI lies within ±0.8 SD);
  - *inconclusive* (neither).
- Reported for each program: the pooled estimate with 95% and 90% CIs, the prediction interval, I², the number of datasets agreeing in direction, and the minimal detectable effect at 80% power.

**Checks:**

1. **Positive control:** X-inactivation-escape program (KDM6A, KDM5C, DDX3X, EIF1AX, ZFX, USP9X, JPX). It must be higher in women; otherwise the pipeline is considered to have failed.
2. **Calibration:** sex labels are shuffled within dataset × assay 1,000 times. The false-positive rate at α = 0.05 and permutation p values are reported.
3. **Sensitivity analyses:**
   - adding a dissociation/ischemia stress score as a covariate;
   - adding log capillary cell number as a covariate;
   - leaving each dataset out in turn.
4. **Genome-wide:** limma-voom per dataset, and inverse-variance meta-analysis of gene-level sex effects. XIST and Y-linked genes serve as positive controls. Between-dataset correlation of autosomal effects is reported.

## Secondary outcomes

**Aim 2: structure and signaling** (Holm correction within each family)

- *Composition:*
  - pericytes per capillary EC, as log((pericytes + 0.5) / (capillary ECs + 0.5));
  - capillary share of endothelial cells (logit);
  - endothelial share of all cells (logit).

  These are also reported by suspension type (cells vs nuclei).
- *Capillary–pericyte signaling:* a donor score of z(ligand) + z(receptor) for five pairs: PDGFB–PDGFRB, JAG1–NOTCH3, EDN1–EDNRA, ANGPT1–TEK and ANGPT2–TEK.
- Same estimand and equivalence rules as Aim 1.

**Aim 3: menopause proxy (descriptive)**

- Donors under 50 compared with donors 50 or older (primary; the inventory found only 2 eligible women under 45). Sensitivity: under 45 vs 55 or older, excluding 45–54.
- The estimate is the sex × age-band interaction for the six programs and the pericyte ratio. 
- A dataset contributes only with at least 2 women and 2 men in each band. Results are reported as estimates with CIs only.

**Aim 4: sex-chromosome specificity**

- Eight X–Y pairs (KDM6A/UTY, KDM5C/KDM5D, USP9X/USP9Y, DDX3X/DDX3Y, EIF1AX/EIF1AY, ZFX/ZFY, RPS4X/RPS4Y1, NLGN4X/NLGN4Y), plus XIST and JPX.
- Measured in capillary ECs, arterial ECs, cardiomyocytes, fibroblasts, pericytes, myeloid cells and smooth muscle from the same donors.
- Reported for each pair: the X copy alone, the Y copy alone and X + Y (log2 CPM), and the Y share in men.
- *Specificity test:* a mixed model with a donor random effect; the sex × cell-group interaction uses capillary ECs as the reference. Results are pooled across datasets with BH correction.

## What counts as a result

- **Positive:** Holm p < 0.05, with the same direction in most datasets and passed calibration.
- **Informative null:** equivalence shown (90% CI within ±0.8 SD). The claim is then that large baseline differences are unlikely, not that there are none.
- **Inconclusive:** anything else, reported as such.

## Not available in these data (stated as limitations)

Cause of death, ischemic time, menopausal status, hormone therapy, medications and ancestry.

## Deviations

*(none yet)*
