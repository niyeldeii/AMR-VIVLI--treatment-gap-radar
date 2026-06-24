# Treatment Gap Radar — code

Integrates 13 Vivli AMR surveillance datasets with the Global AMR R&D Hub export to compare a
**Resistance Need Index (RNI)** against an **R&D Attention Index (RAI)** and surface pathogen-level
"treatment gaps" (high resistance, low R&D).

## Setup
```bash
conda create -n tgr python=3.12 -y && conda activate tgr
pip install "numpy==1.26.4" "pandas==2.2.2" openpyxl xlrd matplotlib seaborn \
            plotly streamlit pycountry rapidfuzz jupyter pyarrow pyyaml
```
Raw data is **not** in this repo (Vivli DUA). It must sit in the parent challenge folder; paths are
in `src/paths.py`.

## Run (from this `treatment_gap_radar/` folder)
```bash
python -m src.pipeline                 # harmonize -> indicators -> RNI/RAI -> gap (writes data_processed/*.parquet)
python -m src.pipeline --no-harmonize  # reuse existing isolates_long.parquet
streamlit run app/dashboard.py         # interactive dashboard
python notebooks/build_notebook.py     # (re)generate the analysis notebook
jupyter nbconvert --to notebook --execute --inplace notebooks/Treatment_Gap_Radar.ipynb
```

## Layout
| Path | Purpose |
|---|---|
| `config/` | canonical drugs, pathogens, indicator weights |
| `src/canon.py`, `mic.py`, `breakpoints.py` | normalization, MIC parsing, CLSI/TB breakpoints |
| `src/loaders.py`, `harmonize.py` | per-dataset loaders -> unified long table |
| `src/indicators.py` | six resistance indicators |
| `src/rni.py`, `rai.py`, `gap.py` | the two indices + gap quadrants |
| `src/viz.py` | Plotly figures (shared by notebook + dashboard) |
| `app/dashboard.py` | Streamlit app |
| `notebooks/` | analysis notebook (built from `build_notebook.py`) |
| `PROFILING_NOTES.md` | per-dataset schema findings |

## Headline finding
Priority treatment gaps are **neglected Enterobacterales** (*Providencia*, *Klebsiella aerogenes*,
*Proteus*, *Serratia*, *Citrobacter*): high resistance, minimal R&D. WHO-critical headline pathogens
(*A. baumannii*, *K. pneumoniae*) and TB are well-served. Sub-Saharan Africa is a surveillance blind
spot despite >80% 3GC resistance in the SPIDAAR RWE cohort.

## Method caveats
- Only ATLAS carries native S/I/R; others use a curated CLSI/TB breakpoint table (`src/breakpoints.py`).
- RAI is a coarse, text-derived attribution (R&D Hub tags are mostly Gram-class + funder geography),
  so the gap is interpreted at pathogen level.
