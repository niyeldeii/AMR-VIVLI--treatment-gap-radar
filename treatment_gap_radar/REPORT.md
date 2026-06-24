# Treatment Gap Radar: Aligning Antimicrobial R&D with Where Resistance Is Worsening

*Vivli AMR Surveillance Open Data Re-use Data Challenge — AMR ID 00013367*

---

## Abstract

Antimicrobial resistance (AMR) is rising faster than new therapies reach the clinic, and research
and development (R&D) is finite. We built the **Treatment Gap Radar**, a reproducible framework that
integrates 13 industry AMR surveillance datasets from the Vivli AMR Register (~1.0 million isolates;
17.9 million isolate–antibiotic results, 2004–2025) with the Global AMR R&D Hub investment and
pipeline data ($18.9 billion across 19,023 projects). For each pathogen we compute a six-indicator
**Resistance Need Index (RNI)** and compare it to an **R&D Attention Index (RAI)** to identify
"treatment gaps" — pathogens where resistance pressure is high but R&D attention is low. We validate
the framework with bootstrap confidence intervals, logistic time-trend models, PCA-justified
weighting, weight-sensitivity analysis, and a machine-learning model that predicts resistance in
surveillance blind spots (held-out-country R² = 0.73). The clearest gaps are Gram-positive
(*Enterococcus faecium*/VRE, *Staphylococcus epidermidis*); broad-spectrum Gram-negative R&D already
covers the Enterobacterales; and the largest blind spot is geographic — Sub-Saharan Africa shows
>80% third-generation-cephalosporin resistance yet is barely surveilled.

## 1. Background and objective

Surveillance tells us *where resistance is*; R&D databases tell us *where investment goes*. They are
rarely examined together. Our objective: **systematically test whether global antimicrobial R&D is
aimed at the pathogens whose resistance is actually escalating**, and to flag misalignments useful
for stewardship and funding prioritization.

## 2. Data

**Resistance (Vivli AMR Register, 13 datasets):** ATLAS (antibiotics), GEARS, SIDERO-WT, KEYSTONE,
SOAR (×3), Innoviva *Acinetobacter*, J&J DREAM (MDR-TB), Venus PLEA I/II & GASAR, and the SPIDAAR
real-world-evidence cohort (Kenya, Ghana, Uganda, Malawi). ATLAS dominates (1,011,168 isolates,
83 countries, 2004–2024). **R&D (Global AMR R&D Hub):** 19,023 projects, $18.9B, with funder,
amount, product/pipeline, research area, and infectious-agent tags.

## 3. Methods

**Harmonization.** Each dataset was normalized into one long table
(isolate × antibiotic × susceptibility) with canonical pathogen, antibiotic (and drug class), and
ISO-3 country. ATLAS carries native S/I/R; for the MIC-only datasets we derived S/I/R with a curated
CLSI M100 breakpoint table (eight organism groups) and tuberculosis critical concentrations. The
SPIDAAR phenotype flags (3GC-resistance, MRSA) were mapped to representative agents.

**Resistance Need Index (RNI).** Per pathogen we computed six indicators — resistance prevalence,
time trend (MIC drift), multidrug-resistance frequency, geographic spread, therapeutic scarcity,
and pediatric involvement — each min–max normalized and combined (default equal weights).

**R&D Attention Index (RAI).** R&D Hub records resolve mainly to Gram class, so each pathogen
inherits its Gram-class broad-spectrum investment/projects/pipeline (Gram-negative $4.9B,
Gram-positive $3.0B, tuberculosis $4.6B) plus a bonus for work naming the species. Clinical
(therapeutics/diagnostics) work is up-weighted.

**Gap.** RNI vs RAI defines four quadrants; the priority quadrant is high need / low attention.

**Statistical rigor.** Bootstrap 95% CIs for prevalence; logistic regression of resistance on year
(odds ratio per year + 95% CI); PCA on the indicators to test whether equal weighting is justified;
and a weight-sensitivity analysis (Spearman rank correlation across weighting schemes).

**Blind-spot prediction (ML).** A gradient-boosted model predicts resistance from *generalizable*
features (pathogen, drug, drug class, Gram, continent, income tier, year) — deliberately excluding
country identity — validated by holding out entire countries (GroupKFold), so its accuracy on unseen
countries measures its ability to estimate resistance where surveillance is absent.

**Forecasting / early warning.** For each pathogen–drug pair we project the logistic resistance–year
fit forward and solve for the year %R crosses 50%, propagating the slope's 95% CI into a CI on that
crossing year.

**R&D attribution robustness.** Because R&D attribution is the framework's most contestable step, we
recompute the priority gaps under three schemes — species-named only, Gram-class only, and class +
named — and report rank stability.

**External validation.** We compare our resistance estimates to published WHO GLASS 2022 global
figures.

## 4. Results

- **Validity.** Independent datasets reproduce ATLAS: *A. baumannii* meropenem resistance 60.9%
  (95% CI 60.4–61.3) in ATLAS vs 62.3% in Innoviva; *P. aeruginosa* meropenem 20–21% across ATLAS
  and GEARS. Known epidemiology holds (carbapenem resistance *A. baumannii* ≫ *E. coli*).
- **Resistance is worsening.** *A. baumannii* meropenem resistance rises with an odds ratio of ~1.07
  per year (p < 0.001).
- **Treatment gaps.** The clearest high-need/low-attention pathogens are **Gram-positive**:
  *Enterococcus faecium* (VRE) and *Staphylococcus epidermidis*. Broad-spectrum Gram-negative R&D
  already covers the Enterobacterales (*Providencia*, *Proteus*, *Serratia*, *Citrobacter*), which a
  naïve name-only attribution would wrongly flag. WHO-critical *A. baumannii*/*K. pneumoniae* and TB
  are well-served (high need, high attention).
- **Methodology is robust.** PCA's first component explains ~72% of indicator variance with near-equal
  loadings (equal weighting justified); RNI rankings are stable across weighting schemes
  (Spearman 0.98–0.99).
- **Geographic blind spots.** The prediction model attains R² = 0.73 on held-out countries
  (weighted MAE 8.1% vs 18.4% for a global-mean baseline). Sub-Saharan Africa is barely surveilled,
  yet SPIDAAR shows 82% (*E. coli*) and 90% (*K. pneumoniae*) ceftriaxone resistance, and the model
  predicts high resistance across under-surveilled regions.
- **Early warning.** *Klebsiella pneumoniae* carbapenem (meropenem) resistance is projected to cross
  50% by ≈2036 (95% CI 2035–2038), rising from ~21% today; *A. baumannii* is already above 50% and
  projected at ~86% by 2035. Declining pathogens (MRSA, VRE, pneumococcal penicillin) show no crossing.
- **The gap is robust to R&D attribution.** Under both defensible attribution schemes (class, class +
  named) the priority gaps are identical (*E. faecium*, TB, *S. epidermidis*; Spearman 1.0/0.87);
  only naive name-only attribution differs (0.33) and is rejected.
- **External validation.** MRSA matches GLASS almost exactly (36.6% vs 35%). Our 3GC estimates sit
  below GLASS medians (ATLAS high-income sampling bias) while the Africa subset brackets the high
  end — corroborating both our pipeline and the blind-spot argument.

## 5. Stewardship and policy implications

1. **Targeted R&D** for the Gram-positive gaps (VRE, coagulase-negative staphylococci).
2. **Sustained investment** in the well-served critical Gram-negatives, where resistance is still
   climbing — high attention is necessary but not yet winning.
3. **Expand surveillance in LMICs**, where the model flags the highest *unmeasured* resistance; the
   gap there is a double failure of monitoring *and* R&D reach.

## 6. Limitations

- Only ATLAS provides native S/I/R; other datasets rely on a curated, simplified CLSI/TB breakpoint
  table (documented in `src/breakpoints.py`), not for clinical use.
- R&D Hub records are tagged mostly at Gram-class level and by funder geography; R&D attention is
  therefore attributed at Gram-class level (plus a species bonus), and we do not over-claim
  species-level R&D differences.
- Continents and income tiers used by the prediction model are coarse, documented approximations.
- Surveillance datasets over-represent high-income hospital settings, biasing global estimates.

## 7. Reproducibility and open science

All code is open at the project repository. Raw Vivli data is **not** redistributed (per the Data
Use Agreement); only aggregate results (minimum cell count n ≥ 30) are shared. Rebuild with
`python -m src.pipeline` (framework) and `python -m src.build_rigor` (rigor + ML); explore via the
Streamlit dashboard and the analysis notebook.

## Acknowledgement

This publication or presentation is based on research using data from GSK, Innoviva Specialty
Therapeutics, Johnson & Johnson, Paratek, Pfizer, Shionogi, Venatorx, Venus Remedies Limited,
obtained through https://amr.vivli.org
