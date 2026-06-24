# Treatment Gap Radar: a reproducible framework aligning antimicrobial R&D attention with where resistance is worsening

*Formatted for submission to the Wellcome Open Research — Vivli AMR Open Data Re-use Data Challenge collection.*

**Authors:** Olakulehin Adebusuyi¹, *[co-authors]* ; on behalf of the Treatment Gap Radar team
**Affiliations:** ¹ Faculty of Pharmaceutical Sciences, University of Ibadan, Nigeria
**Corresponding author:** olakulehin13@gmail.com
**Article type:** Research Article / Method Article
**Keywords:** antimicrobial resistance; surveillance; research and development; treatment gap; machine learning; forecasting; health equity; Vivli; AMR R&D Hub

---

## Abstract

**Background.** Antimicrobial resistance (AMR) is rising while antibiotic research and development
(R&D) is limited and concentrated on a few high-profile pathogens. Surveillance data and R&D
investment data are rarely analysed together, so there is no systematic, data-driven check of
whether R&D is directed at the pathogens whose resistance is escalating fastest.

**Methods.** We integrated 13 industry AMR surveillance datasets from the Vivli AMR Register
(≈1.0 million isolates; 17.9 million isolate–antibiotic results, 2004–2025) with the Global AMR R&D
Hub investment and pipeline data ($18.9 billion; 19,023 projects). For each pathogen we computed a
six-indicator **Resistance Need Index (RNI)** and an **R&D Attention Index (RAI)**, and defined the
treatment gap as their difference. We added bootstrap confidence intervals, logistic time-trend
models, principal-component and weight-sensitivity analyses, a logistic early-warning forecast of
threshold crossing, a gradient-boosted model to predict resistance in surveillance blind spots
(validated on held-out countries), an R&D-attribution robustness analysis, and external validation
against WHO GLASS 2022.

**Results.** Resistance estimates reproduced across independent datasets and matched WHO GLASS for
MRSA (36.6% vs 35%). The clearest treatment gaps were Gram-positive — *Enterococcus faecium* (VRE)
and *Staphylococcus epidermidis* — high resistance need with limited targeted R&D; broad-spectrum
Gram-negative R&D already covers the Enterobacterales; WHO-critical *Acinetobacter baumannii*,
*Klebsiella pneumoniae* and tuberculosis are well-served. The priority gaps were identical under both
defensible R&D-attribution schemes (Spearman 1.0/0.87). *K. pneumoniae* carbapenem resistance is
projected to cross 50% by ≈2036 (95% CI 2035–2038). The blind-spot model reached R² = 0.73 on unseen
countries; Sub-Saharan Africa is barely surveilled yet shows >80% third-generation-cephalosporin
resistance.

**Conclusions.** A reproducible "radar" can flag where R&D is misaligned with resistance need and
where surveillance is blind. It points to under-addressed Gram-positive pathogens, sustained need in
critical Gram-negatives, and an urgent equity gap in low- and middle-income-country surveillance.

---

## 1. Introduction

AMR is among the leading global health threats, yet the antibiotic pipeline is thin and economically
fragile. Investment and development effort concentrate on a small set of pathogens, while resistance
evolves across many. Surveillance programmes (e.g. the Vivli AMR Register) and R&D trackers (e.g. the
Global AMR R&D Hub) exist in parallel but are seldom analysed jointly. We ask a simple, policy-relevant
question: **is global antimicrobial R&D attention aligned with where resistance is actually getting
worse?** We answer it by constructing two comparable per-pathogen indices and overlaying them, then
stress-testing the result and extending it with forecasting and prediction.

## 2. Methods

### 2.1 Data
Thirteen Vivli AMR Register datasets were used: ATLAS (Pfizer), GEARS (Venatorx), SIDERO-WT
(Shionogi), KEYSTONE (Paratek), SOAR 201818/201910/207965 (GSK), Innoviva *Acinetobacter*
surveillance, DREAM MDR-tuberculosis (Johnson & Johnson), Venus PLEA I/II and GASAR III, and the
SPIDAAR real-world-evidence cohort (Pfizer; Kenya, Ghana, Uganda, Malawi). R&D data were the Global
AMR R&D Hub export (19,023 projects; $18.9 billion).

### 2.2 Harmonization
Datasets were normalized into a single long table (isolate × antibiotic × susceptibility) with
canonical pathogen, antibiotic and drug class, and ISO-3166 country. ATLAS provided native
susceptible/intermediate/resistant (S/I/R) interpretation; for the 12 minimum-inhibitory-concentration
(MIC)-only datasets, S/I/R was derived using a curated CLSI M100 breakpoint table covering eight
organism groups, plus tuberculosis critical concentrations.

### 2.3 Resistance Need Index (RNI)
Six indicators per pathogen — resistance prevalence, time trend (MIC drift), multidrug-resistance
frequency, geographic spread, therapeutic scarcity, and pediatric involvement — were each min–max
normalized to 0–1 and combined (equal weights by default). Equal weighting was tested by principal
component analysis (PCA) and a weight-sensitivity analysis.

### 2.4 R&D Attention Index (RAI)
R&D Hub records resolve mainly to Gram class, so each pathogen inherited its Gram-class
broad-spectrum pool of investment, project count and pipeline products, plus a bonus for projects
naming the species. Clinical (therapeutics/diagnostics) work was up-weighted. Robustness was assessed
by recomputing the gap under three attribution schemes (species-named only; Gram-class only;
class + named).

### 2.5 Gap, rigor, forecasting, prediction, validation
The gap (RNI − RAI) defined four quadrants. We computed bootstrap 95% CIs for prevalence and logistic
regressions of resistance on calendar year (odds ratio per year, 95% CI). An early-warning forecast
projected each logistic trajectory forward and solved for the year %R crosses 50%, propagating the
slope CI. A gradient-boosted regressor predicted %R from generalizable features (pathogen, drug, drug
class, Gram, continent, World Bank income tier, year), validated by holding out whole countries
(GroupKFold). Estimates were compared to WHO GLASS 2022 figures. Finally, a causal-inference-style
panel model (pathogen × year resistance; year × R&D-attention-tertile interaction; pathogen fixed
effects) tested whether R&D attention relates to the *rate* of resistance growth, with a transparent
counterfactual scenario.

## 3. Results

### 3.1 Validity
Independent datasets reproduced ATLAS (*A. baumannii* meropenem 60.9% vs 62.3% in Innoviva;
*P. aeruginosa* meropenem 20–21% across ATLAS and GEARS). Against WHO GLASS 2022, MRSA matched
(36.6% vs 35%); our third-generation-cephalosporin estimates sat below GLASS medians, consistent with
ATLAS's high-income hospital sampling, while the Sub-Saharan-Africa subset (82–90%) bracketed the
high end.

### 3.2 The treatment gap
The clearest priority gaps were Gram-positive — *Enterococcus faecium* (VRE) and *Staphylococcus
epidermidis*. Broad-spectrum Gram-negative R&D already covers the Enterobacterales; WHO-critical
*A. baumannii*, *K. pneumoniae* and tuberculosis were well-served. Priority gaps were identical under
the two defensible attribution schemes (Spearman 1.0 and 0.87 versus the combined scheme); naive
name-only attribution differed (0.33) and was rejected.

### 3.3 Rigor
PCA's first component explained ≈72% of indicator variance with near-equal loadings, justifying equal
weighting; rankings were stable across weighting schemes (Spearman 0.98–0.99). Resistance trends were
significant for critical Gram-negatives (e.g. *A. baumannii* meropenem OR 1.07/year, p < 0.001).

### 3.4 Early warning
*K. pneumoniae* carbapenem (meropenem) resistance is projected to cross 50% by ≈2036 (95% CI
2035–2038), from ~21% currently; *A. baumannii* is already above 50% and projected to ≈86% by 2035.
Declining agents (MRSA, VRE, pneumococcal penicillin) showed no crossing.

### 3.5 Predicting blind spots
The gradient-boosted model attained R² = 0.73 on held-out countries (weighted MAE 8.1% vs 18.4% for a
global-mean baseline), enabling resistance estimates for under-surveilled regions.

### 3.6 R&D is reactive, not pre-emptive
In the panel model, pathogens with the **highest** R&D attention had **rising** resistance
(+0.18 %R/year, p < 0.001), while lower-attention pathogens were stable or declining, after adjusting
for each pathogen's baseline (fixed effects). This is consistent with R&D attention being allocated
*reactively* — concentrated where resistance is already worsening — rather than pre-emptively. The
illustrative counterfactual indicates substantial avertable resistance (e.g. ≈31 percentage points
for *K. pneumoniae* carbapenem by 2035 under best-case growth), underscoring the value of acting
earlier and in the blind spots.

## 4. Discussion

The framework converts heterogeneous surveillance and R&D data into a single, defensible answer to a
policy question. Three implications follow: (1) targeted R&D for under-addressed Gram-positive
pathogens (VRE); (2) sustained investment in critical Gram-negatives, where resistance continues to
climb; and (3) urgent expansion of surveillance in low- and middle-income countries, where the model
predicts the highest *unmeasured* resistance. The forecasting and blind-spot components give the radar
a forward-looking, equity-oriented dimension beyond description.

### 4.1 Limitations
Only ATLAS provided native S/I/R; other datasets used a curated, simplified CLSI/TB breakpoint table
(not for clinical use); resistance can be sensitive to breakpoint vs epidemiological-cutoff (ECOFF)
choice. R&D Hub records are coarsely tagged (mostly Gram class) and by funder geography, so R&D
attention is attributed at Gram-class level (plus a species bonus); we do not over-claim species-level
R&D differences. Continents and income tiers in the prediction model are coarse approximations.
Surveillance over-represents high-income hospital settings, biasing global estimates downward —
quantified here against GLASS.

### 4.2 Conclusions
The Treatment Gap Radar is an open, reproducible tool that flags misalignment between resistance need
and R&D attention and highlights surveillance blind spots, supporting more strategic AMR R&D
prioritization and stewardship.

## 5. Data availability
Source data are third-party datasets accessed through the Vivli AMR Register (https://amr.vivli.org)
under the Vivli Data Use Agreement and cannot be redistributed. Aggregate results underlying the
figures (minimum cell count n ≥ 30) are shared in the project repository. The SPIDAAR dataset is
licensed CC BY 4.0.

## 6. Software/code availability
All analysis code is openly available: https://github.com/niyeldeii/AMR-VIVLI--treatment-gap-radar
(reproduce with `python -m src.pipeline`). Built with Python (pandas, scikit-learn, statsmodels,
plotly, Streamlit).

## 7. Grant information / Funding
The 2025 Vivli AMR Surveillance Data Challenge was supported by Johnson & Johnson, Paratek, Pfizer and
a U.S. National Institutes of Health award. *[Add any team-specific funding.]*

## 8. Acknowledgements
This publication or presentation is based on research using data from GSK, Innoviva Specialty
Therapeutics, Johnson & Johnson, Paratek, Pfizer, Shionogi, Venatorx, Venus Remedies Limited, obtained
through https://amr.vivli.org.

## References (indicative)
1. World Health Organization. *Global Antimicrobial Resistance and Use Surveillance System (GLASS)
   Report 2022.* Geneva: WHO; 2022.
2. World Health Organization. *WHO Bacterial Priority Pathogens List, 2024.* Geneva: WHO; 2024.
3. Clinical and Laboratory Standards Institute. *M100 Performance Standards for Antimicrobial
   Susceptibility Testing.* CLSI.
4. Global AMR R&D Hub. *Dynamic Dashboard of AMR R&D investments and pipeline.* https://globalamrhub.org
5. Vivli AMR Register. https://amr.vivli.org
