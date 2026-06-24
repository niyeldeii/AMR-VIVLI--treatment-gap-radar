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
rai[["gram", "projects", "investment", "pipeline", "named_projects", "RAI"]].round(2)
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

md("""## 7 · Statistical rigor
We stress-test the framework so the conclusions are defensible, not arbitrary.""")

md("""**(a) Resistance is rising — quantified.** Logistic regression of resistance on year gives an
odds ratio per year with a 95% confidence interval for each key pathogen–drug pair.""")
code("viz.trend_forest()")
code("pd.read_parquet(PROCESSED_DIR / 'trend_models.parquet').round(4)")

md("""**(b) The index weighting is justified, not arbitrary.** Principal component analysis shows the
six indicators load on one dominant axis (PC1 explains ~72% of variance), and the data-driven
weights are close to equal — so equal weighting is defensible.""")
code("""
from src.models import pca_weights
pca = pca_weights()
print("PC1 explains %.0f%% of variance" % (pca['explained_variance'][0]*100))
pd.Series(pca['weights'], name='PCA-derived weight').round(3)
""")

md("**(c) Conclusions are robust to the weighting choice** (Spearman rank correlation vs equal weights).")
code("pd.read_parquet(PROCESSED_DIR / 'sensitivity.parquet')")

md("""**(d) Uncertainty is reported.** Bootstrap 95% confidence intervals for resistance prevalence
of key pathogen–drug pairs.""")
code("pd.read_parquet(PROCESSED_DIR / 'ci_keypairs.parquet').round(2)")

md("""## 8 · Predicting surveillance blind spots (machine learning)
A gradient-boosted model predicts resistance from *generalizable* features (pathogen, drug, Gram,
continent, income tier, year) — never country identity — so it can estimate resistance for
countries/regions with **no surveillance**. We validate by holding out whole countries.""")
code("""
rs = pd.read_parquet(PROCESSED_DIR / 'rigor_summary.parquet').iloc[0]
print(f"Held-out-country cross-validation R2 : {rs['cv_R2_unseen_countries']:.2f}")
print(f"Weighted mean absolute error          : {rs['cv_weighted_MAE']*100:.1f}%")
print(f"Naive baseline (global mean) MAE      : {rs['baseline_MAE_global_mean']*100:.1f}%")
""")
code("viz.blindspot_continent('Klebsiella pneumoniae', 'Meropenem')")

md("""## 9 · Key findings

- **The clearest treatment gaps are Gram-positive:** ***Enterococcus faecium* (VRE)** — high
  resistance need, little targeted R&D — and ***Staphylococcus epidermidis***. These are genuine
  gaps, not artifacts of how R&D is tagged.
- **Broad-spectrum Gram-negative R&D covers the Enterobacterales**, including the obscure ones
  (*Providencia*, *Proteus*, *Serratia*, *Citrobacter*). Attributing R&D only to *named* species
  would falsely flag these — our class-aware model corrects that.
- **WHO-critical headliners (*A. baumannii*, *K. pneumoniae*) and TB are well-served** — high need
  *and* high R&D attention. The system works where it focuses.
- **Resistance is actively worsening** for the critical Gram-negatives (e.g. *A. baumannii*
  meropenem, OR ≈ 1.07 per year, p < 0.001).
- **The biggest blind spot is geographic.** Sub-Saharan Africa is barely represented in surveillance,
  yet the SPIDAAR cohort shows **>80% third-generation-cephalosporin resistance** in *E. coli* and
  *K. pneumoniae*, and our model predicts high resistance across under-surveilled regions. High
  burden, low visibility — a double failure of surveillance *and* R&D reach.

### Stewardship implications
Prioritize (1) targeted R&D for VRE / Gram-positive gaps, (2) sustained investment in the
well-served critical Gram-negatives where resistance is still climbing, and (3) **expanding
surveillance in LMICs**, where the model flags the highest unmeasured resistance.

### Caveats (stated plainly)
- Only ATLAS ships native S/I/R; other datasets use a curated CLSI/TB breakpoint table
  (approximate, documented in `src/breakpoints.py`).
- R&D Hub records are tagged mostly at Gram-class level and by funder geography, so R&D attention
  is attributed at the Gram-class level (with a species-specific bonus); we do **not** over-claim
  species-level R&D differences.
- Income tiers/continents used in the prediction model are coarse, documented approximations.
""")

nb["cells"] = cells
out = Path(__file__).resolve().parent / "Treatment_Gap_Radar.ipynb"
nbf.write(nb, out)
print("wrote", out)
