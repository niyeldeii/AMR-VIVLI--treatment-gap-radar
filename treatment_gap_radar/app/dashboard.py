"""Treatment Gap Radar — interactive dashboard.

Run:  streamlit run treatment_gap_radar/app/dashboard.py
Reads the precomputed parquet outputs (build them first with src.pipeline).
"""
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # treatment_gap_radar/
from src import viz                                               # noqa: E402
from src.paths import PROCESSED_DIR                               # noqa: E402

st.set_page_config(page_title="Treatment Gap Radar", layout="wide")


@st.cache_data
def load(name):
    return pd.read_parquet(PROCESSED_DIR / name)


if not (PROCESSED_DIR / "gap.parquet").exists():
    st.error("Outputs not found. Run:  python -m src.pipeline   (from the treatment_gap_radar folder)")
    st.stop()

gap = load("gap.parquet")
rni = load("rni.parquet")
combo = load("combo_resistance.parquet")

st.title("🛰️ Treatment Gap Radar")
st.caption("Antimicrobial Resistance Need vs Global R&D Attention — Vivli AMR Data Challenge (AMR ID 00013367)")

tab_gap, tab_path, tab_trend, tab_blind, tab_cov, tab_about = st.tabs(
    ["Gap Radar", "Pathogen explorer", "Trends & rigor", "Blind-spot prediction",
     "Surveillance coverage", "About / data"])

with tab_gap:
    c1, c2, c3 = st.columns(3)
    n_gap = (gap["quadrant"].str.startswith("PRIORITY")).sum()
    c1.metric("Pathogens analyzed", len(gap))
    c2.metric("Priority treatment gaps", int(n_gap))
    c3.metric("Total isolates", f"{int(gap['n_isolates'].sum()):,}")
    st.plotly_chart(viz.gap_quadrant(gap), use_container_width=True)
    st.subheader("Priority gaps (high resistance need, low R&D attention)")
    pg = gap[gap["quadrant"].str.startswith("PRIORITY")].sort_values("gap_score", ascending=False)
    st.dataframe(pg[["who", "n_isolates", "RNI", "RAI", "gap_score"]].round(3),
                 use_container_width=True)
    with st.expander("Indicator heatmap (all pathogens)"):
        st.plotly_chart(viz.indicator_heatmap(rni), use_container_width=True)

with tab_path:
    pathogens = sorted(combo["pathogen"].dropna().unique())
    default = pathogens.index("Acinetobacter baumannii") if "Acinetobacter baumannii" in pathogens else 0
    patho = st.selectbox("Pathogen", pathogens, index=default)
    drugs = sorted(combo[combo["pathogen"] == patho]["antibiotic"].dropna().unique())
    drug = st.selectbox("Antibiotic", drugs,
                        index=(drugs.index("Meropenem") if "Meropenem" in drugs else 0))
    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(viz.resistance_choropleth(patho, drug), use_container_width=True)
    with col2:
        st.plotly_chart(viz.resistance_trend(patho, drug), use_container_width=True)
    if patho in rni.index:
        r = rni.loc[patho]
        st.write(f"**{patho}** — RNI {r['RNI']:.2f}  ·  WHO tier: {r.get('who') or '—'}  "
                 f"·  isolates {int(r['n_isolates']):,}")

with tab_trend:
    st.subheader("Resistance trends over time")
    st.caption("Logistic regression of resistance on year — odds ratio per year with 95% CIs. "
               "OR > 1 means resistance is rising.")
    if (PROCESSED_DIR / "trend_models.parquet").exists():
        st.plotly_chart(viz.trend_forest(), use_container_width=True)
        st.dataframe(load("trend_models.parquet").round(4), use_container_width=True)
    st.divider()
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Methodology is robust")
        if (PROCESSED_DIR / "rigor_summary.parquet").exists():
            rs = load("rigor_summary.parquet").iloc[0]
            st.metric("PCA — variance on 1st component", f"{rs['pc1_var']*100:.0f}%",
                      help="Indicators load on one dominant 'need' axis, so equal weighting is justified.")
            st.write("**Data-driven (PCA) weights vs equal weights:**")
            wts = {k.replace("w_", ""): rs[k] for k in rs.index if k.startswith("w_")}
            st.dataframe(pd.Series(wts, name="PCA weight").round(3))
    with c2:
        st.subheader("Rankings are weight-insensitive")
        if (PROCESSED_DIR / "sensitivity.parquet").exists():
            st.caption("Spearman correlation of the pathogen ranking under different weightings "
                       "vs the equal-weight default. ~1.0 = conclusions don't depend on weights.")
            st.dataframe(load("sensitivity.parquet"), use_container_width=True)

with tab_blind:
    st.subheader("Predicting resistance where there is no surveillance")
    if (PROCESSED_DIR / "rigor_summary.parquet").exists():
        rs = load("rigor_summary.parquet").iloc[0]
        m1, m2, m3 = st.columns(3)
        m1.metric("Model R² on unseen countries", f"{rs['cv_R2_unseen_countries']:.2f}",
                  help="Held-out-country cross-validation: how well resistance is predicted "
                       "for countries the model never saw.")
        m2.metric("Mean abs. error", f"{rs['cv_weighted_MAE']*100:.1f}%")
        m3.metric("vs naive baseline", f"{rs['baseline_MAE_global_mean']*100:.1f}%")
    if (PROCESSED_DIR / "blindspot_predictions.parquet").exists():
        preds = load("blindspot_predictions.parquet")
        paths = sorted(preds["pathogen"].unique())
        dflt = paths.index("Klebsiella pneumoniae") if "Klebsiella pneumoniae" in paths else 0
        bp = st.selectbox("Pathogen", paths, index=dflt, key="bp_path")
        bdrugs = sorted(preds[preds["pathogen"] == bp]["antibiotic"].unique())
        bd = st.selectbox("Antibiotic", bdrugs, key="bp_drug")
        st.plotly_chart(viz.blindspot_continent(bp, bd), use_container_width=True)
        st.caption("Red bars = regions with thin/absent surveillance (≤2 countries reporting) — "
                   "the model's best estimate of resistance there.")

with tab_cov:
    st.plotly_chart(viz.coverage_map(), use_container_width=True)
    st.info("White / pale countries are surveillance blind spots — little or no isolate data, "
            "regardless of likely clinical burden (esp. low- and middle-income countries).")

with tab_about:
    st.markdown("""
**Resistance Need Index (RNI)** combines six indicators (prevalence, MIC drift, MDR frequency,
geographic spread, therapeutic scarcity, pediatric involvement) from the Vivli surveillance data.
**R&D Attention Index (RAI)** is derived from the Global AMR R&D Hub export (investment, projects,
pipeline products). The gap = RNI − RAI.

*Caveat:* R&D records are tagged coarsely (mostly Gram class) and by funder geography, so the
need-vs-attention comparison is at pathogen level and is an approximate attribution.

> Data from GSK, Innoviva Specialty Therapeutics, Johnson & Johnson, Paratek, Pfizer, Shionogi,
> Venatorx, Venus Remedies Limited, obtained through https://amr.vivli.org
""")
