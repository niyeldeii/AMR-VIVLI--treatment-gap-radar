# Treatment Gap Radar — Implementation Plan

## Context

This is a **Vivli AMR Data Challenge** submission (AMR ID 00013367, "Treatment Gap Radar").
The approved proposal (`DrSnapshot-00013367-...pdf`) asks us to answer one question:

> **Is global antimicrobial R&D effort (funding + drug pipelines) going where antimicrobial
> resistance is actually getting worse?**

We answer it by building and overlaying two indices:

1. **Resistance Need Index (RNI)** — from the Vivli surveillance datasets (isolate-level MIC /
   S-I-R lab data).
2. **R&D Attention Index (RAI)** — from `Projects.xlsx`, which (confirmed by inspection) **is a
   Global AMR R&D Hub export**: funder, amount USD/EUR, product/pipeline name, research area
   (Basic / Therapeutics / Diagnostics / Operational), and infectious-agent tags.

Both halves are already present locally. The gap between them = the "treatment gaps" the radar
surfaces.

**Decisions locked in with the user:**
- Deliverable: **both** a reproducible Jupyter notebook *and* an interactive dashboard.
- Scope: **full framework** — all eligible datasets, all six resistance indicators.
- Pathogens: **all bacteria** (Gram-negative, Gram-positive, TB). **Exclude fungal**
  (`ATLAS_Antifungals`) from the core; keep as an optional appendix only.

## Data inventory (already on disk, NOT in git)

Root: the challenge folder containing this repo. Data referenced by relative path; never committed.

| Dataset | File | Format | Notes |
|---|---|---|---|
| ATLAS Antibiotics | `ATLAS_Antibiotics/atlas_vivli_2004_2024.csv` (~370 MB) | CSV | **Primary.** Wide; MIC + `_I` (S/I/R) per drug; Species, Country, Year, Age Group, Source |
| SOAR 201818 | `SOAR 201818/gsk_201818_published.csv` | CSV | Wide; MIC only (no interpretation col); ORGANISMNAME, COUNTRY, REGION, YEARCOLLECTED, AGE |
| SOAR 201910 / 207965 | `SOAR 201910/...xlsx`, `SOAR 207965/...xlsx` | XLSX | Need profiling |
| Innoviva Acinetobacter | `Surveillance .../IST-Entasis_Acinetobacter...xlsx` | XLSX | A. baumannii focus |
| DREAM (Bedaquiline/MDR-TB) | `Bedaquiline.../...xlsx` | XLSX | TB; different drug panel |
| KEYSTONE (Omadacycline) | `KEYSTONE/Omadacycline...xlsx` | XLSX | |
| SIDERO-WT (Cefiderocol) | `SIDERO-WT/...xlsx` | XLSX | Gram-negative, cefiderocol |
| GEARS | `GEARS/Venatorx...xlsx` | XLSX | |
| PLEA I / II, GASAR III | `PLEA.../...xlsx`, `GASAR.../...xlsx` | XLSX | Venus; protocol PDFs included |
| SPIDAAR RWE | `SPIDAAR.../spidaar_isolatedata.xls` + `_patientdata.xls` | **XLS (legacy)** | needs `xlrd`; separate isolate + patient tables |
| **R&D Hub** | `Projects.xlsx` (sheet `data`, ~58 MB) | XLSX | The RAI source |
| ATLAS Antifungals | `ATLAS_Antifungals/...xlsx` | XLSX | **Out of core scope** |

## Prerequisite: fix the Python environment (blocker)

Current `base` conda env is broken: **pandas 2.2.2 cannot import under NumPy 2.5.0** (binary
built for NumPy 1.x), and **`xlrd` is missing** (required for the legacy `.xls` SPIDAAR files).

Fix in a dedicated env so we don't disturb `base`:

```bash
conda create -n tgr python=3.12 -y
conda activate tgr
pip install "numpy==1.26.4" "pandas==2.2.2" openpyxl xlrd \
            matplotlib seaborn plotly streamlit pycountry rapidfuzz jupyter
```

(`numpy==1.26.4` is the safe pin against pandas 2.2.2; `pycountry` for ISO country
normalization; `rapidfuzz` for fuzzy pathogen/drug-name matching across datasets.)

## Proposed project structure

```
treatment_gap_radar/
  config/
    pathogens.yaml        # canonical pathogen names + Gram class + WHO priority tier
    antibiotics.yaml      # canonical drug names + drug class + priority-drug flag
    dataset_map.yaml      # per-dataset column → canonical-field mapping
    weights.yaml          # indicator weights for RNI/RAI
  src/
    loaders.py            # one loader per dataset → tidy DataFrame
    harmonize.py          # wide→long, name normalization, S/I/R derivation
    breakpoints.py        # CLSI/EUCAST MIC breakpoints where _I absent
    indicators.py         # the six resistance indicators
    rni.py                # assemble Resistance Need Index
    rai.py                # R&D Attention Index from Projects.xlsx
    gap.py                # join RNI vs RAI, quadrant classification
  data_processed/         # parquet outputs (gitignored)
  figures/                # exported static charts/maps
  notebooks/
    Treatment_Gap_Radar.ipynb   # end-to-end narrative deliverable
  app/
    dashboard.py          # Streamlit app
```

## Phase 1 — Profiling & harmonization (the real engineering work)

Every dataset has a different schema. Build a **common long-format isolate table**:

`source | isolate_id | pathogen | gram_class | country | region | year | age | age_group | specimen | antibiotic | drug_class | mic | interpretation(S/I/R)`

Steps:
1. **Profile every file**: dump columns, dtypes, unique pathogens/countries/years, which drug
   columns exist, whether interpretation columns exist. Record into `config/dataset_map.yaml`.
2. **Normalize entities**:
   - Pathogen names → canonical (`pathogens.yaml`), tag Gram class + WHO priority tier.
   - Antibiotic names → canonical (`antibiotics.yaml`), tag drug class + priority-drug flag.
     ATLAS uses `Drug` + `Drug_I`; SOAR uses bare MIC columns. Use `rapidfuzz` for near-matches.
   - Countries → ISO-3 via `pycountry`; keep region.
3. **Derive S/I/R**: prefer the dataset's own `_I` column (ATLAS). Where absent (SOAR etc.),
   parse censored MIC strings (`<=0.015`, `>8`) to numeric and apply breakpoints in
   `breakpoints.py`. **Caveat to document:** breakpoint coverage is partial; indicators are
   computed only on (pathogen, drug) pairs with a defined breakpoint or a provided `_I`.
4. Write `data_processed/isolates_long.parquet`.

## Phase 2 — Resistance Need Index (six indicators)

Computed per **pathogen × drug × country × year**, then aggregated (`src/indicators.py`):

1. **Resistance prevalence** — %R (and %NS = I+R) among tested isolates.
2. **MIC drift** — slope of MIC50/MIC90 (or geometric-mean MIC) over `year`; positive = worsening.
3. **MDR frequency** — % isolates resistant to ≥3 drug classes (class map from `antibiotics.yaml`).
4. **Geographic spread** — share of surveyed countries where %R exceeds a threshold.
5. **Therapeutic scarcity** — among priority drugs for that pathogen, count/fraction still
   effective (low %R); fewer effective options ⇒ higher scarcity.
6. **Pediatric involvement** — %R within pediatric age bands (from `Age Group`/`AGE`).

Each indicator min-max normalized to 0–1; **RNI = weighted sum** (`config/weights.yaml`,
default equal weights, sensitivity check in notebook). Output `data_processed/rni.parquet`,
aggregated to pathogen level and pathogen×country level.

## Phase 3 — R&D Attention Index

From `Projects.xlsx` sheet `data` (`src/rai.py`):
- Map each project to a pathogen/Gram group via `Infectious Agent` + `Individual Infectious Agent`
  (mostly `Gram negative` / `Gram positive`, some specific e.g. `Aspergillus`, TB via `Disease`).
- Aggregate per pathogen group: **project count, total `Amount USD`, count of pipeline products**
  (`Product Name` non-null), split by `Research Area` (Therapeutics / Diagnostics / Basic).
- Normalize to 0–1 ⇒ **RAI**. Optionally weight Therapeutics/Diagnostics over Basic.
- Output `data_processed/rai.parquet`.

**Granularity caveat (document prominently):** R&D Hub records are tagged by *funder/institution
country* and a coarse infectious-agent category — **not** by the country where resistance occurs,
and rarely below Gram-class. Therefore:
- The **RNI-vs-RAI gap is computed at pathogen / Gram-class level** (the level both sides share).
- **Geographic "blind spots"** are derived from *surveillance coverage* (which countries are
  under-sampled in the Vivli data), not from R&D geography.

## Phase 4 — Gap analysis & outputs

`src/gap.py`: join RNI and RAI at pathogen level → quadrant classification:
- **High need / Low attention** → priority treatment gaps (headline result).
- High need / High attention → well-served.
- Low need / High attention → possible over-investment.
- Low need / Low attention → low priority.

Gap score = `RNI − RAI` (both normalized). Output `data_processed/gap.parquet`.

**Deliverable A — Notebook** (`notebooks/Treatment_Gap_Radar.ipynb`): narrative from raw data →
harmonization → six indicators → RNI → RAI → gap quadrants, with the headline figures:
- Gap **quadrant scatter** (RNI x-axis, RAI y-axis, bubble = isolate count).
- Resistance **choropleth maps** (Plotly) per key pathogen.
- **MIC-drift** trend lines for priority pathogen-drug pairs.
- Surveillance-coverage map (blind spots).
- The required acknowledgement text (GSK, Innoviva, J&J, Paratek, Pfizer, Shionogi, Venatorx,
  Venus, via amr.vivli.org).

**Deliverable B — Dashboard** (`app/dashboard.py`, Streamlit): reads the parquet outputs;
filters by pathogen / drug / country / year; renders the gap radar (quadrant + the resistance
choropleth + indicator drill-downs). Reuses functions from `src/` — no recomputation logic
duplicated in the app.

## Verification

1. **Env**: `python -c "import pandas, numpy, xlrd, streamlit, plotly; print('ok')"` succeeds.
2. **ETL sanity**: after Phase 1, assert every loader returns rows; print per-dataset row counts,
   pathogen/country/year coverage. Spot-check censored-MIC parsing on known rows.
3. **Indicator sanity**: all %R in [0,1]; known epidemiology holds (e.g. carbapenem resistance
   markedly higher in *A. baumannii* than *E. coli*).
4. **Index sanity**: RNI/RAI in [0,1]; re-run with alternate `weights.yaml` to confirm quadrant
   conclusions are stable (sensitivity analysis recorded in notebook).
5. **Notebook**: `jupyter nbconvert --to notebook --execute notebooks/Treatment_Gap_Radar.ipynb`
   runs top-to-bottom with no errors.
6. **Dashboard**: `streamlit run app/dashboard.py`, load it, exercise filters, screenshot the gap
   radar to confirm it renders against the processed data.

## Open items to confirm during build (not blockers)
- Exact column layouts of the un-profiled XLSX/XLS datasets (resolved in Phase 1 profiling).
- Breakpoint source: default to **CLSI**; note EUCAST alternative. Coverage gaps documented.
- Final indicator weights (start equal; tune with sensitivity analysis).
