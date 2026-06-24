# Dataset Profiling Notes

Generated from a first-pass schema scan of every dataset (header + sample rows). These notes
drive the harmonization layer (`src/loaders.py`, `src/harmonize.py`). Raw data is **not** in
this repo — paths are relative to the challenge folder one level up.

## Cross-cutting findings

| Concern | Detail | Action |
|---|---|---|
| S/I/R availability | **Only ATLAS** has interpretation (`*_I`) columns. All others are raw MIC. | Derive S/I/R from MIC via CLSI breakpoints (`breakpoints.py`) for non-ATLAS. |
| Drug naming | 3 conventions: full names; abbrev codes (SOAR201910 `AMC/AMP`, GEARS `CAZ_MIC`, DREAM `INH/RMP`); newline-broken names (KEYSTONE). | Canonical drug map in `config/antibiotics.yaml` keyed on normalized (strip newlines/whitespace, lowercase) names + an alias table for codes. |
| Species naming | `E. coli` vs `Escherichia coli`; `M. tuberculosis`; GEARS appends phenotype (`Staphylococcus aureus, MSSA`). | Canonical pathogen map in `config/pathogens.yaml`; strip phenotype suffix; expand abbreviations; `rapidfuzz` fallback. |
| Year/date formats | plain `Year`; epoch-ns datetime (SOAR201910 `Collection Date`); `YYYYMMDD` int (SIDERO `Date Collected`); `YearCollected` int. | Per-dataset year extractor in loader; normalize to integer `year`. |
| Header row | PLEA I & II have data starting below row 0 (all `Unnamed:` cols). | Detect header row (first row with >50% non-null string cells) before parse. |
| MIC censoring | Values like `<=0.015`, `>8`, `>2`. | Parse to numeric + censor flag (`mic_op` in {<=,>=,<,>,=}); keep raw string too. |

## Per-dataset

- **ATLAS_Antibiotics** (CSV, ~370 MB, 127 cols) — PRIMARY. Cols: `Isolate Id, Study, Species,
  Country, State, Gender, Age Group, Speciality, Source, Year`, then per-drug `Drug` + `Drug_I`
  (47 interp cols), then β-lactamase gene columns (`NDM, KPC, OXA, CTX-M-*, VIM, IMP, SHV, TEM`…).
  Read with `usecols` / chunks given size. Year range 2004–2024 (sample showed 2004 only).
- **SOAR 201818** (CSV, 24 cols) — MIC only. `ORGANISMNAME, COUNTRY, REGION, YEARCOLLECTED, AGE,
  DEID_CAT_AGE, BETALACTAMASE`, full-name drug cols. S. pneumoniae + H. influenzae. 9 countries.
- **SOAR 201910** (xlsx, sheet `3550 valid MIC data (2)`, 26 cols) — MIC only, **abbrev drug
  codes** (`AMC,AMP,AMX,AXO,AZM,CDN,CEC,CLA,CXM,DIN,ERY,FIX,LEV,MXF,PEN,POD,SXT`). `Collection
  Date` is epoch-ns datetime. `Organism, Country, Age`.
- **SOAR 207965** (xlsx, sheet `Sheet2`, 37 cols) — MIC only, full-name drugs. Has `GramType`,
  `OriginalOrganismName/FinalOrganismName/OrganismFamilyName`, `YearCollected` (2018–2021).
- **Innoviva Acinetobacter** (xlsx `Sheet1`, 19 cols) — MIC only. Acinetobacter spp. (baumannii
  dominant). `OrganismName, Country (33), YearCollected 2016–2021`, drugs incl. Cefiderocol-era
  panel + `Sulbactam/514`.
- **DREAM TB** (xlsx `DREAM Dataset`, 27 cols) — *M. tuberculosis*. Drug codes
  `INH,RMP,LZD,CFZ,LVX,OFX,MXF,CAP,KAN,AMI,EMB` + `BDQ Broth/BDQ MGIT` + resistance-gene mutation
  cols (`Rv0678, atpE, pepQ, Rv1979c`). TB needs its own breakpoint logic / treat separately.
- **KEYSTONE** (xlsx `Line List`, 46 cols) — Gram-positive heavy (S. aureus, Enterococcus). MIC
  only, **column names contain `\n`**. US-only (`US Census Division`), rich clinical context
  (ICU, VAP, infection source). Omadacycline focus.
- **SIDERO-WT** (xlsx `Five year Surveillance data` [+ Japanese `作業用` sheet], 20 cols) —
  Gram-negative, cefiderocol focus. MIC only. `Date Collected` = YYYYMMDD int. 34 countries.
- **GEARS** (xlsx `Data` [+ `Define` dictionary sheet], 23 cols) — Gram-negative. MIC cols with
  `_MIC` suffix + codes (`CAZ_MIC,FEP_MIC,IPM_MIC,MEM_MIC,TZP_MIC`…). 11 countries.
- **PLEA Study I** (xlsx `Sheet 1`) & **PLEA Study II** (xlsx `MIC of Ploymyxin`) — Venus polymyxin
  studies. **Header not in row 0** — needs header detection. Decode against the protocol PDFs.
- **GASAR Study III** (xlsx `Sheet1` [+ `SpecimenSources`], 8 cols) — `Species` (abbrev: E. coli,
  K. pneumoniae, A. baumannii, P. aeruginosa), `Polymyxin B MIC`, `Gene Combination`,
  `Phenotypic Combination`, Year 2022–2023, 1 country (India).
- **SPIDAAR RWE** (xls, legacy) — **patient-level real-world-evidence**, not MIC surveillance.
  `spidaar_isolatedata` (`data` sheet): coded flags `c3r, mdr, mrsa, amrtx, stype, ctry`.
  `spidaar_patientdata`: 60 cols of clinical/demographic (`agegr, ctry, mdr-ish, comorbidities`).
  Has `codebook`/`definitions`/`variable labels` sheets — **decode via those before use**. Best
  used for MDR-frequency and pediatric indicators, joined isolate↔patient on id.

## Implications for the plan
- Breakpoint coverage is the critical path: ATLAS gives free S/I/R; everything else depends on
  CLSI tables per (pathogen, drug). TB (DREAM) and polymyxin-only (Venus) need special handling.
- SPIDAAR is structurally different (RWE) — fold in only for MDR% / pediatric, don't force it into
  the per-drug MIC schema.
- Build a single `config/dataset_map.yaml` capturing, per dataset: sheet, header row, id col,
  pathogen col, country col, year extractor, drug-column list, naming convention.
