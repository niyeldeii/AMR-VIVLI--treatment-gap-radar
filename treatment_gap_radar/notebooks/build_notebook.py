"""Assemble Treatment_Gap_Radar.ipynb from source cells (reproducible).

Run from the treatment_gap_radar folder:  python notebooks/build_notebook.py
Then execute:  jupyter nbconvert --to notebook --execute --inplace notebooks/Treatment_Gap_Radar.ipynb
"""
import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip()))
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip()))

md("""
# 🛰️ Treatment Gap Radar
### Does global antimicrobial R&D go where resistance is actually getting worse?
*Vivli AMR Data Challenge — AMR ID 00013367*

This notebook integrates **13 antimicrobial-resistance surveillance datasets** (Vivli) with the
**Global AMR R&D Hub** investment/pipeline data to build two indices and compare them:

- **Resistance Need Index (RNI)** — six indicators of resistance pressure per pathogen.
- **R&D Attention Index (RAI)** — investment, project count, and pipeline products per pathogen.

The **gap** (RNI − RAI) surfaces pathogens where resistance is high but R&D attention is low.
""")

md("""
> *This publication or presentation is based on research using data from GSK, Innoviva Specialty
> Therapeutics, Johnson & Johnson, Paratek, Pfizer, Shionogi, Venatorx, Venus Remedies Limited,
> obtained through https://amr.vivli.org*
""")

code("""
import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()))

import pandas as pd
import plotly.io as pio
pio.renderers.default = "notebook_connected"   # embed figures as HTML (CDN plotly.js)
from src import viz
from src.paths import PROCESSED_DIR

# Build outputs if missing (reuses isolates_long.parquet if already harmonized).
if not (PROCESSED_DIR / "gap.parquet").exists():
    from src import pipeline
    pipeline.main(skip_harmonize=(PROCESSED_DIR / "isolates_long.parquet").exists())

long_cols = ["source", "pathogen", "antibiotic", "sir", "year", "country_iso3"]
print("Outputs ready in", PROCESSED_DIR)
""")

md("## 1 · Data harmonization\nAll datasets are normalized into one long table (isolate × antibiotic × S/I/R).")

code("""
import pyarrow.parquet as pq
pf = pq.ParquetFile(PROCESSED_DIR / "isolates_long.parquet")
print(f"Harmonized rows: {pf.metadata.num_rows:,}")
src_counts = pd.read_parquet(PROCESSED_DIR / "isolates_long.parquet", columns=["source"])["source"].value_counts()
src_counts
""")

md("## 2 · Resistance indicators → Resistance Need Index\nSix indicators per pathogen, normalized 0–1 and combined.")

code("""
rni = pd.read_parquet(PROCESSED_DIR / "rni.parquet")
rni[["who", "n_isolates", "prevalence", "mic_drift", "mdr", "geo_spread", "scarcity", "pediatric", "RNI"]].round(3)
""")

code("viz.indicator_heatmap(rni)")

md("## 3 · R&D Attention Index\nText-attributed investment, project count, and pipeline products per pathogen (Global AMR R&D Hub).")

code("""
rai = pd.read_parquet(PROCESSED_DIR / "rai.parquet")
rai[["projects", "investment", "pipeline", "therapeutics", "RAI"]].round(2)
""")

md("""## 4 · The Treatment Gap
RNI (x) vs RAI (y). The **red quadrant** = high resistance need, low R&D attention — the priority gaps.""")

code("viz.gap_quadrant()")

code("""
gap = pd.read_parquet(PROCESSED_DIR / "gap.parquet")
print("PRIORITY TREATMENT GAPS (high need / low attention):")
gap[gap["quadrant"].str.startswith("PRIORITY")].sort_values("gap_score", ascending=False)[
    ["who", "n_isolates", "RNI", "RAI", "gap_score"]].round(3)
""")

md("## 5 · Pathogen deep-dives\nGeographic distribution and time trend for a critical pathogen–drug pair.")

code("viz.resistance_choropleth('Acinetobacter baumannii', 'Meropenem')")
code("viz.resistance_trend('Klebsiella pneumoniae', 'Meropenem')")

md("""## 6 · Surveillance blind spots & LMICs
Data volume by country. Blind spots (white) mark low-visibility regions — often LMICs where burden
may be high. The SPIDAAR real-world dataset (Kenya/Ghana/Uganda/Malawi) shows the point sharply.""")

code("viz.coverage_map()")

code("""
long = pd.read_parquet(PROCESSED_DIR / "isolates_long.parquet",
                       columns=["source", "pathogen", "antibiotic", "sir", "country"])
sp = long[(long["source"] == "SPIDAAR") & (long["antibiotic"] == "Ceftriaxone")]
print("SPIDAAR (Sub-Saharan Africa) — % ceftriaxone-resistant (3GC-R):")
(sp.groupby("pathogen")["sir"].apply(lambda s: round((s == "R").mean()*100, 1))
   .rename("%R").reset_index())
""")

md("""## 7 · Key findings

- **Priority treatment gaps** are dominated by **neglected Enterobacterales** — *Providencia*,
  *Klebsiella aerogenes*, *Proteus mirabilis*, *Serratia marcescens*, *Citrobacter freundii* —
  which carry substantial resistance but attract minimal R&D investment.
- The WHO-critical headline pathogens (*A. baumannii*, *K. pneumoniae*, *E. cloacae*) and **TB**
  are **well-served**: high need *and* high R&D attention.
- *S. aureus*, *P. aeruginosa* and *E. coli* receive R&D attention well above their measured need
  in this surveillance footprint.
- **Blind spots:** Sub-Saharan Africa is barely represented in the surveillance datasets, yet the
  SPIDAAR RWE cohort shows **>80% 3GC resistance** in *E. coli*/*K. pneumoniae* — high burden,
  low visibility.

### Caveats
- Only ATLAS ships native S/I/R; other datasets use a curated CLSI/TB breakpoint table (approximate,
  documented in `src/breakpoints.py`).
- R&D Hub records are tagged coarsely (mostly Gram class) and by funder geography, so RAI is an
  approximate, text-derived attribution; the gap is interpreted at pathogen level.
- Indicator weights are equal by default (`config/weights.yaml`); see sensitivity below.
""")

code("""
# Sensitivity: re-rank with prevalence-weighted RNI to confirm the priority gaps are stable.
from src.rni import compute_rni
w = {"prevalence": 3, "mic_drift": 1, "mdr": 1, "geo_spread": 1, "scarcity": 1, "pediatric": 1}
alt = compute_rni(weights=w, save=False).sort_values("RNI", ascending=False)
alt[["who", "RNI"]].head(8).round(3)
""")

nb["cells"] = cells
out = Path(__file__).resolve().parent / "Treatment_Gap_Radar.ipynb"
nbf.write(nb, out)
print("wrote", out)
