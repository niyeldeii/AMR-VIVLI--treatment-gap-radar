# 🛰️ Treatment Gap Radar

**Does global antimicrobial R&D go where resistance is actually getting worse?**

A reproducible framework for the **Vivli AMR Surveillance Open Data Re-use Data Challenge**
(AMR ID 00013367) that integrates **13 antimicrobial-resistance (AMR) surveillance datasets** with
**Global AMR R&D Hub** investment/pipeline data to measure, for every pathogen, whether research and
development attention is aligned with resistance need — and to flag the gaps.

> 📄 Full write-up: [`treatment_gap_radar/REPORT.md`](treatment_gap_radar/REPORT.md)  ·
> 📓 Notebook: [`treatment_gap_radar/notebooks/Treatment_Gap_Radar.ipynb`](treatment_gap_radar/notebooks/Treatment_Gap_Radar.ipynb)  ·
> 🛰️ Dashboard: `streamlit run treatment_gap_radar/app/dashboard.py`

---

## 1. The problem

Resistance keeps rising; drug development is slow, expensive, and finite. Funders and companies
concentrate on a handful of "famous" pathogens. Nobody routinely checks, in a data-driven way,
whether R&D money flows to the bacteria whose resistance is genuinely escalating. The Treatment Gap
Radar answers that question by putting **two scores per pathogen** on the same map.

## 2. The framework

| Index | Question | Built from |
|---|---|---|
| **Resistance Need Index (RNI)** | *How bad is this bug getting?* | Six indicators from Vivli surveillance |
| **R&D Attention Index (RAI)** | *How much is the world investing against it?* | Global AMR R&D Hub |
| **Gap = RNI − RAI** | *Is attention aligned with need?* | The overlay |

**The six RNI indicators:** resistance prevalence · time trend (MIC drift) · multidrug-resistance
frequency · geographic spread · therapeutic scarcity · pediatric involvement. Each is min–max
normalized to 0–1 and combined (equal weights by default — *validated by PCA*, see §6).

**Gap quadrants:**

| Quadrant | Interpretation |
|---|---|
| 🔴 High need / low attention | **Treatment gaps** (priority) |
| 🟢 High need / high attention | Well-served |
| 🔵 Low need / high attention | Possible over-investment |
| ⚪ Low need / low attention | Lower priority |

## 3. Data

**Resistance — 13 Vivli AMR Register datasets (~1.0 M isolates; 17.9 M isolate–antibiotic results, 2004–2025):**

| Dataset | Contributor | Notes |
|---|---|---|
| ATLAS Antibiotics | Pfizer | Primary; 1,011,168 isolates, 83 countries, native S/I/R |
| GEARS | Venatorx | Gram-negative |
| SIDERO-WT | Shionogi | Cefiderocol-era Gram-negatives |
| KEYSTONE | Paratek | Omadacycline; Gram-positive heavy |
| SOAR ×3 (201818/201910/207965) | GSK | Respiratory (*S. pneumoniae*, *H. influenzae*) |
| Innoviva *Acinetobacter* surveillance | Innoviva | *A. baumannii* focus |
| DREAM | Johnson & Johnson | MDR-tuberculosis (bedaquiline) |
| PLEA I/II, GASAR III | Venus Remedies | Polymyxin studies |
| SPIDAAR RWE | Pfizer | Real-world cohort, Sub-Saharan Africa |

**R&D — Global AMR R&D Hub:** 19,023 projects, **$18.9 billion**, with funder, amount, product
pipeline, research area, and infectious-agent tags.

## 4. What we built (pipeline)

```
raw data ─▶ harmonize ─▶ isolates_long (17.9M rows) ─▶ indicators ─▶ RNI ─┐
Projects.xlsx ─▶ R&D Attention Index (RAI) ───────────────────────────────┼─▶ GAP ─▶ deliverables
                                                                          │
              rigor (CIs, trends, PCA, sensitivity) + ML (blind-spot model)┘
```

**Harmonization.** Every dataset has a different schema, drug-naming convention, date format and
header layout. We normalize all of them into one long table
(`isolate × antibiotic × susceptibility`) with canonical pathogen, antibiotic + drug class, and
ISO-3 country (`src/canon.py`, `src/loaders.py`, `src/harmonize.py`).

**Susceptibility (S/I/R).** ATLAS ships native interpretation. For the 12 MIC-only datasets we
derived S/I/R with a **curated CLSI M100 breakpoint table** (8 organism groups) plus
**tuberculosis critical concentrations** (`src/breakpoints.py`). Validated by cross-dataset
agreement (see §6).

**RNI** (`src/indicators.py`, `src/rni.py`) — the six indicators above per pathogen.

**RAI** (`src/rai.py`) — R&D Hub records resolve mostly to Gram class, so each pathogen inherits its
**Gram-class broad-spectrum pool** (Gram-negative $4.9 B · Gram-positive $3.0 B · TB $4.6 B) **plus a
bonus for work naming the species**. This is the defensible attribution: a broad-spectrum
Gram-negative antibiotic genuinely treats *Providencia* even when that genus is never named.

**Gap** (`src/gap.py`) — RNI vs RAI → quadrants and a gap score.

## 5. Statistical rigor (so it's defensible, not arbitrary)

`src/models.py`, precomputed by the pipeline (`src/pipeline.py`):

- **Bootstrap 95% CIs** for resistance prevalence (e.g. *A. baumannii* meropenem 60.9%, CI 60.4–61.3).
- **Logistic time-trend models** — odds ratio of resistance per year + 95% CI + p-value
  (*A. baumannii* meropenem **OR ≈ 1.07/yr, p < 0.001** → resistance rising).
- **PCA-justified weighting** — the first principal component explains **~72%** of indicator
  variance with near-equal loadings, so equal weighting is **data-justified**, not arbitrary.
- **Weight-sensitivity analysis** — RNI rankings are stable across weighting schemes
  (**Spearman 0.98–0.99**).

## 6. Machine learning — predicting surveillance blind spots

`src/predict.py` (+ `src/geo.py`). A gradient-boosted model predicts resistance from
**generalizable** features (pathogen, drug, drug class, Gram, continent, income tier, year) —
**never country identity** — so it can estimate resistance for countries/regions with **no
surveillance**. Validated by holding out **entire countries** (GroupKFold):

- **R² = 0.73** on unseen countries · weighted MAE **8.1%** vs **18.4%** for a global-mean baseline.

This turns the framework into a tool that can *estimate* resistance where data is missing.

## 6b. Causal layer — is R&D reactive?

`src/causal.py`. A pathogen × year resistance panel with a year × R&D-attention-tertile interaction
and **pathogen fixed effects** tests whether R&D attention relates to the *rate* of resistance growth.
Result: the **highest-R&D pathogens are those still rising** (+0.18 %R/yr, p<0.001) while
lower-attention pathogens are stable/declining — R&D **chases** resistance rather than pre-empting it
(observational, not causal proof). A transparent counterfactual quantifies avertable 2035 resistance.

## 7. Key findings

- **The clearest treatment gaps are Gram-positive:** *Enterococcus faecium* (VRE) and
  *Staphylococcus epidermidis* — high resistance need, little targeted R&D.
- **Broad-spectrum Gram-negative R&D already covers the Enterobacterales** (*Providencia*,
  *Proteus*, *Serratia*, *Citrobacter*) — attributing R&D only to *named* species would falsely
  flag these; our class-aware model corrects it.
- **WHO-critical *A. baumannii* / *K. pneumoniae* and TB are well-served** (high need, high attention).
- **Resistance is actively worsening** for critical Gram-negatives (significant positive time trends).
- **The biggest blind spot is geographic.** Sub-Saharan Africa is barely surveilled, yet SPIDAAR
  shows **82% (*E. coli*) and 90% (*K. pneumoniae*) ceftriaxone resistance** — high burden, low
  visibility, a double failure of surveillance *and* R&D reach.
- **R&D is reactive, not pre-emptive** (causal layer). The highest-R&D pathogens are those whose
  resistance is still *rising* (+0.18%/yr, p<0.001); lower-attention pathogens are stable or
  declining — investment chases resistance instead of getting ahead of it.

### Stewardship implications
1. Targeted R&D for the Gram-positive gaps (VRE).
2. Sustained investment in the well-served critical Gram-negatives, where resistance is still climbing.
3. **Expand LMIC surveillance**, where the model flags the highest unmeasured resistance.

## 8. Deliverables

| Deliverable | Path |
|---|---|
| 📄 Manuscript-style report | `treatment_gap_radar/REPORT.md` |
| 📓 Reproducible analysis notebook (+ HTML export) | `treatment_gap_radar/notebooks/Treatment_Gap_Radar.ipynb` / `.html` |
| 🛰️ Interactive Streamlit dashboard | `treatment_gap_radar/app/dashboard.py` (6 tabs) |
| 🧱 Documented schema findings | `treatment_gap_radar/PROFILING_NOTES.md` |

## 9. How to run

```bash
# 1. Environment
conda create -n tgr python=3.12 -y && conda activate tgr
pip install "numpy" "pandas==2.2.2" openpyxl xlrd matplotlib seaborn plotly streamlit \
            pycountry rapidfuzz jupyter pyarrow pyyaml scikit-learn statsmodels

# 2. Place the raw Vivli datasets in the parent challenge folder (see src/paths.py). NOT in this repo.

# 3. Build everything with ONE command (from the treatment_gap_radar/ folder)
python -m src.pipeline                 # harmonize -> indicators -> RNI/RAI/gap -> rigor + ML
#   python -m src.pipeline --no-harmonize   # reuse existing isolates_long.parquet
#   python -m src.pipeline --fast           # skip the slower blind-spot ML model

# 4. Explore
streamlit run app/dashboard.py
jupyter notebook notebooks/Treatment_Gap_Radar.ipynb
```

## 10. Repository layout

```
treatment_gap_radar/
  config/        canonical antibiotics, pathogens, indicator weights
  src/
    paths.py         data locations
    canon.py         pathogen/drug/country normalization (+ fuzzy match)
    mic.py           censored-MIC parsing (<=, >, unicode ≤/≥)
    breakpoints.py   curated CLSI M100 + TB critical concentrations
    loaders.py       per-dataset loaders -> common long schema
    harmonize.py     unify all sources -> isolates_long.parquet
    indicators.py    the six resistance indicators
    rni.py / rai.py / gap.py   the two indices + gap
    models.py        trends, PCA weights, sensitivity, bootstrap CIs
    forecast.py      early-warning: project trajectories to 50% threshold
    geo.py / predict.py   blind-spot prediction model
    causal.py        is R&D reactive? panel growth model + counterfactual
    validate.py      external validation vs WHO GLASS 2022
    util.py          shared helpers (minmax)
    viz.py           shared Plotly figures
    pipeline.py      single end-to-end build (harmonize -> analysis -> rigor + ML + causal)
  app/dashboard.py   Streamlit app
  notebooks/         analysis notebook + builder
  data_processed/    aggregate results (parquet) — raw/isolate-level data excluded
  REPORT.md, PROFILING_NOTES.md, README.md
```

## 11. Data governance

The raw Vivli surveillance data is **restricted under the Vivli Data Use Agreement** and is **not**
in this repository — an allowlist `.gitignore` and pre-commit guards ensure it is never committed.
Only **aggregate results** (minimum cell count **n ≥ 30**) are published, consistent with the
challenge's open-science requirement. The 135 MB isolate-level `isolates_long.parquet` stays local.

## 12. Limitations (stated plainly)

- Only ATLAS provides native S/I/R; other datasets use a curated, simplified CLSI/TB breakpoint
  table (documented; not for clinical use).
- R&D Hub records are tagged mostly at Gram-class level and by funder geography; R&D attention is
  attributed at Gram-class level (plus a species bonus). We do not over-claim species-level R&D.
- Continents / income tiers in the prediction model are coarse, documented approximations.
- Surveillance over-represents high-income hospital settings, biasing global estimates.

## Acknowledgement

> This publication or presentation is based on research using data from GSK, Innoviva Specialty
> Therapeutics, Johnson & Johnson, Paratek, Pfizer, Shionogi, Venatorx, Venus Remedies Limited,
> obtained through https://amr.vivli.org

*Built with Claude Code.*
