"""Treatment Gap Radar — interactive dashboard.

Run:  streamlit run treatment_gap_radar/app/dashboard.py
Reads precomputed parquet outputs (build them with `python -m src.pipeline`).
"""
import importlib
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # treatment_gap_radar/
from src import viz                                               # noqa: E402
from src.paths import PROCESSED_DIR                               # noqa: E402

importlib.reload(viz)   # pick up the latest viz on Streamlit Cloud hot-reload (avoid stale module)

st.set_page_config(page_title="Treatment Gap Radar", layout="wide", page_icon="🛰️")
STRETCH = "stretch"


@st.cache_data
def load(name):
    return pd.read_parquet(PROCESSED_DIR / name)


def have(name):
    return (PROCESSED_DIR / name).exists()


def chart(fig, key):
    st.plotly_chart(fig, width=STRETCH, key=key)


def table(df, **kw):
    st.dataframe(df, width=STRETCH, **kw)


def dl(df, label, fname):
    st.download_button(label, df.to_csv().encode(), file_name=fname, mime="text/csv", key="dl_" + fname)


if not have("gap.parquet"):
    st.error("Outputs not found. Run `python -m src.pipeline` from the treatment_gap_radar folder.")
    st.stop()

gap = load("gap.parquet")
rni = load("rni.parquet")
combo = load("combo_resistance.parquet")
rs = load("rigor_summary.parquet").iloc[0] if have("rigor_summary.parquet") else None

# ---------------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("🛰️ Treatment Gap Radar")
    st.caption("Vivli AMR Surveillance Open Data Re-use Data Challenge · AMR ID 00013367")
    st.metric("Pathogens scored", len(gap))
    st.metric("Priority treatment gaps", int(gap["quadrant"].str.startswith("PRIORITY").sum()))
    st.metric("Isolates analyzed", f"{int(gap['n_isolates'].sum()):,}")
    if rs is not None and "cv_R2_unseen_countries" in rs:
        st.metric("Blind-spot model R² (unseen countries)", f"{rs['cv_R2_unseen_countries']:.2f}")
    st.markdown("---")
    st.markdown("**Data:** 13 Vivli AMR datasets + Global AMR R&D Hub ($18.9 B, 19,023 projects).")
    st.markdown("[GitHub repository](https://github.com/niyeldeii/AMR-VIVLI--treatment-gap-radar)")

st.title("🛰️ Treatment Gap Radar")
st.markdown("#### Does global antimicrobial R&D go where resistance is actually getting worse?")

tabs = st.tabs(["Overview", "Gap Radar", "Pathogen explorer", "Trends · rigor · causal",
                "Blind-spot prediction", "Surveillance coverage", "Methodology",
                "Findings & implications", "About / data"])

# ============================================================= 0. OVERVIEW
with tabs[0]:
    st.markdown("""
The **Treatment Gap Radar** scores every pathogen on two axes and overlays them:
a **Resistance Need Index** (how bad resistance is getting, from six surveillance indicators) and
an **R&D Attention Index** (investment + pipeline from the Global AMR R&D Hub). The gap between them
surfaces where the world is *under-developing* therapies relative to need.
""")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Surveillance datasets", "13")
    c2.metric("Isolate–antibiotic results", "17.9 M")
    c3.metric("R&D tracked", "$18.9 B")
    if rs is not None:
        c4.metric("Indicators on 1 PCA axis", f"{rs['pc1_var']*100:.0f}%")
    st.caption("**How to read the radar:** each bubble is a pathogen — x = resistance need, "
               "y = R&D attention, size = isolates tested. Dashed lines are the medians; the "
               "**bottom-right** quadrant (high need, low attention) holds the priority treatment gaps.")
    chart(viz.gap_quadrant(gap), "ov_gap")
    st.success("**Headline:** the clearest gaps are Gram-positive (VRE / *E. faecium*, "
               "*S. epidermidis*); broad-spectrum Gram-negative R&D already covers the "
               "Enterobacterales; the biggest blind spot is geographic (LMIC surveillance).")

# ============================================================= 1. GAP RADAR
with tabs[1]:
    st.subheader("Resistance Need vs R&D Attention")
    chart(viz.gap_quadrant(gap), "gr_gap")
    colA, colB = st.columns([2, 1])
    with colA:
        st.subheader("Priority gaps (high need, low attention)")
        pg = gap[gap["quadrant"].str.startswith("PRIORITY")].sort_values("gap_score", ascending=False)
        table(pg[["who", "n_isolates", "RNI", "RAI", "gap_score"]].round(3))
        dl(gap.round(4), "⬇ Download full gap table (CSV)", "treatment_gap.csv")
    with colB:
        st.subheader("Quadrant counts")
        table(gap["quadrant"].str.split(" \\(").str[0].value_counts().rename("pathogens"))
    with st.expander("Indicator heatmap — all pathogens"):
        chart(viz.indicator_heatmap(rni), "gr_heat")

# ============================================================= 2. PATHOGEN EXPLORER
with tabs[2]:
    pathogens = sorted(combo["pathogen"].dropna().unique())
    dflt = pathogens.index("Acinetobacter baumannii") if "Acinetobacter baumannii" in pathogens else 0
    patho = st.selectbox("Pathogen", pathogens, index=dflt)
    drugs = sorted(combo[combo["pathogen"] == patho]["antibiotic"].dropna().unique())
    drug = st.selectbox("Antibiotic", drugs,
                        index=(drugs.index("Meropenem") if "Meropenem" in drugs else 0))

    if patho in rni.index:
        r = rni.loc[patho]
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Resistance Need Index", f"{r['RNI']:.2f}")
        if patho in gap.index:
            m2.metric("R&D Attention Index", f"{gap.loc[patho, 'RAI']:.2f}")
            m3.metric("Gap score", f"{gap.loc[patho, 'gap_score']:+.2f}")
        m4.metric("Isolates", f"{int(r['n_isolates']):,}")
        st.caption(f"WHO priority tier: **{r.get('who') or '—'}**")

    c1, c2 = st.columns(2)
    with c1:
        if patho in rni.index:
            chart(viz.indicator_radar(patho, rni), "pe_radar")
        chart(viz.resistance_choropleth(patho, drug), "pe_chor")
    with c2:
        chart(viz.resistance_trend(patho, drug), "pe_trend")
        if have("ci_keypairs.parquet"):
            ci = load("ci_keypairs.parquet")
            row = ci[(ci.pathogen == patho) & (ci.drug == drug)]
            if len(row):
                rr = row.iloc[0]
                st.info(f"**% resistant: {rr['pctR']:.1f}%** (95% CI {rr['lo']:.1f}–{rr['hi']:.1f}, "
                        f"n={int(rr['n']):,})")
        if have("trend_models.parquet"):
            tm = load("trend_models.parquet")
            row = tm[(tm.pathogen == patho) & (tm.drug == drug)]
            if len(row):
                t = row.iloc[0]
                arrow = "rising 📈" if t["OR_per_year"] > 1 else "falling 📉"
                st.warning(f"**Trend: {arrow}** — OR {t['OR_per_year']:.3f}/yr "
                           f"(95% CI {t['ci_lo']:.3f}–{t['ci_hi']:.3f}, p={t['p_value']:.1e})")
    if have("forecast_curves.parquet"):
        fc = load("forecast_curves.parquet")
        if bool(((fc["pathogen"] == patho) & (fc["drug"] == drug)).any()):
            st.markdown("**Early-warning forecast** — observed trajectory + logistic projection to 2035.")
            chart(viz.forecast_plot(patho, drug), "pe_forecast")

# ============================================================= 3. TRENDS & RIGOR
with tabs[3]:
    st.subheader("Resistance trends over time")
    st.caption("Logistic regression of resistance on year — odds ratio per year (95% CI). >1 = rising.")
    if have("trend_models.parquet"):
        chart(viz.trend_forest(), "tr_forest")
        table(load("trend_models.parquet").round(4))
    if have("forecasts.parquet"):
        st.subheader("Early-warning forecast")
        st.caption("Projected %R and the year each *rising* pathogen–drug pair is expected to cross "
                   "50% resistance (logistic projection, 95% CI on the crossing year).")
        ft = load("forecasts.parquet")
        cols = [c for c in ["pathogen", "drug", "pctR_last", "OR_per_year", "pctR_2030",
                            "pctR_2035", "cross50_year", "cross50_lo", "cross50_hi"] if c in ft.columns]
        table(ft[cols])
    st.divider()
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Weighting is data-justified")
        if rs is not None:
            st.metric("Variance on 1st principal component", f"{rs['pc1_var']*100:.0f}%",
                      help="High value = the six indicators measure one underlying 'need' axis, "
                           "so equal weighting is justified.")
            wts = {k.replace("w_", ""): rs[k] for k in rs.index if k.startswith("w_")}
            table(pd.Series(wts, name="PCA-derived weight").round(3))
    with c2:
        st.subheader("Rankings are weight-insensitive")
        if have("sensitivity.parquet"):
            st.caption("Spearman correlation of the pathogen ranking vs equal weights (≈1 = robust).")
            table(load("sensitivity.parquet"))
    if have("ci_keypairs.parquet"):
        with st.expander("Bootstrap 95% confidence intervals (key pairs)"):
            table(load("ci_keypairs.parquet").round(2))
    if have("rai_attribution.parquet"):
        with st.expander("R&D attribution robustness — the gap is stable across schemes"):
            st.caption("Priority gaps under three R&D-attribution schemes. The class-aware schemes "
                       "agree (Spearman ≈ 1); only naive *name-only* attribution differs — and we "
                       "reject it because broad-spectrum R&D demonstrably covers unnamed species.")
            table(load("rai_attribution.parquet"))
    if have("glass_validation.parquet"):
        with st.expander("External validation vs WHO GLASS 2022"):
            table(load("glass_validation.parquet"))
            st.caption("MRSA matches GLASS almost exactly. Our 3GC estimates sit below GLASS medians "
                       "(ATLAS high-income sampling bias); the Sub-Saharan-Africa subset brackets the "
                       "high end — consistent with the surveillance blind-spot argument.")
    if have("causal_growth.parquet"):
        st.divider()
        st.subheader("Causal layer — is R&D reactive?")
        st.caption("Resistance growth rate by R&D-attention tertile, controlling for each pathogen's "
                   "baseline (fixed effects). High-R&D pathogens growing *fastest* indicates "
                   "investment is **chasing** resistance, not pre-empting it (observational).")
        chart(viz.causal_growth_bar(), "causal_bar")
        if have("causal_counterfactual.parquet"):
            with st.expander("Counterfactual — 2035 resistance averted under best-case growth (illustrative)"):
                cf = load("causal_counterfactual.parquet")
                table(cf[["pathogen", "drug", "pctR_2035", "pctR_2035_bestcase", "averted_2035_pts"]])

# ============================================================= 4. BLIND-SPOT PREDICTION
with tabs[4]:
    st.subheader("Predicting resistance where there is no surveillance")
    st.markdown("A gradient-boosted model predicts resistance from *generalizable* features "
                "(pathogen, drug, Gram, continent, income tier, year) — **not** country identity — "
                "validated by holding out whole countries.")
    if rs is not None and "cv_R2_unseen_countries" in rs:
        m1, m2, m3 = st.columns(3)
        m1.metric("R² on unseen countries", f"{rs['cv_R2_unseen_countries']:.2f}")
        m2.metric("Mean absolute error", f"{rs['cv_weighted_MAE']*100:.1f}%")
        m3.metric("vs naive baseline", f"{rs['baseline_MAE_global_mean']*100:.1f}%")
    if have("blindspot_predictions.parquet"):
        preds = load("blindspot_predictions.parquet")
        paths = sorted(preds["pathogen"].unique())
        d0 = paths.index("Klebsiella pneumoniae") if "Klebsiella pneumoniae" in paths else 0
        bp = st.selectbox("Pathogen", paths, index=d0, key="bp_path")
        bdrugs = sorted(preds[preds.pathogen == bp]["antibiotic"].unique())
        bd = st.selectbox("Antibiotic", bdrugs, key="bp_drug")
        chart(viz.blindspot_continent(bp, bd), "bs_chart")
        st.caption("Red bars = thin/absent surveillance (≤2 countries) — the model's best estimate.")
        with st.expander("All predicted blind spots (thin surveillance, highest predicted resistance)"):
            thin = preds[preds["thin_surveillance"]].sort_values("pred_pctR", ascending=False)
            table(thin[["pathogen", "antibiotic", "continent", "n_countries", "pred_pctR"]].round(3))
            dl(preds.round(4), "⬇ Download all predictions (CSV)", "blindspot_predictions.csv")

# ============================================================= 5. COVERAGE
with tabs[5]:
    st.subheader("Surveillance coverage — where are the blind spots?")
    chart(viz.coverage_map(), "cov_map")
    st.info("Pale/white countries are surveillance blind spots — little or no isolate data, "
            "regardless of likely clinical burden (especially low- and middle-income countries). "
            "The SPIDAAR cohort shows >80% 3GC resistance in Sub-Saharan Africa, yet that region "
            "is barely represented in the global surveillance datasets.")

# ============================================================= 6. METHODOLOGY
with tabs[6]:
    st.markdown("""
### How it works

**1. Harmonization.** 13 datasets with different schemas, drug names, date formats and headers are
normalized into one long table (`isolate × antibiotic × susceptibility`, 17.9 M rows) with canonical
pathogen, antibiotic + drug class, and ISO-3 country.

**2. Susceptibility (S/I/R).** ATLAS ships native interpretation; the 12 MIC-only datasets are
interpreted with a curated **CLSI M100** breakpoint table (8 organism groups) plus **tuberculosis
critical concentrations**. Validated by cross-dataset agreement (e.g. *A. baumannii* meropenem 60.9 %
in ATLAS vs 62.3 % in Innoviva).

**3. Resistance Need Index (RNI).** Six indicators per pathogen — prevalence, time trend (MIC drift),
multidrug-resistance frequency, geographic spread, therapeutic scarcity, pediatric involvement —
each normalized 0–1 and combined. Equal weights are **validated by PCA** (PC1 ≈ 72 % variance).

**4. R&D Attention Index (RAI).** R&D Hub records resolve mostly to Gram class, so each pathogen
inherits its Gram-class broad-spectrum pool (GN $4.9 B · GP $3.0 B · TB $4.6 B) **plus** a bonus for
work that names it. This avoids falsely flagging broadly-covered species (e.g. *Providencia*).

**5. Gap.** RNI − RAI → four quadrants; the priority quadrant is high need / low attention.

**6. Rigor & ML.** Bootstrap 95 % CIs; logistic resistance-vs-year trend models; weight-sensitivity
analysis; and a gradient-boosted blind-spot model validated on held-out countries (R² ≈ 0.73).
""")

# ============================================================= 7. FINDINGS
with tabs[7]:
    st.markdown("""
### Key findings

- 🔴 **The clearest treatment gaps are Gram-positive:** *Enterococcus faecium* (VRE) and
  *Staphylococcus epidermidis* — high resistance need, little targeted R&D.
- 🟢 **Broad-spectrum Gram-negative R&D already covers the Enterobacterales** (*Providencia*,
  *Proteus*, *Serratia*, *Citrobacter*). A name-only attribution would falsely flag these.
- 🟢 **WHO-critical *A. baumannii* / *K. pneumoniae* and TB are well-served** (high need, high attention).
- 📈 **Resistance is actively worsening** for critical Gram-negatives (significant positive trends).
- 🌍 **The biggest blind spot is geographic.** Sub-Saharan Africa is barely surveilled, yet SPIDAAR
  shows **82 % (*E. coli*) and 90 % (*K. pneumoniae*) ceftriaxone resistance** — high burden, low
  visibility.
- 🔁 **R&D is reactive, not pre-emptive.** The pathogens receiving the most R&D are precisely those
  whose resistance is still *rising* (+0.18 %/yr, p<0.001); lower-attention pathogens are stable or
  declining. Investment chases resistance rather than getting ahead of it.

### Stewardship & policy implications
1. **Targeted R&D** for the Gram-positive gaps (VRE).
2. **Sustained investment** in the well-served critical Gram-negatives, where resistance is still climbing.
3. **Expand LMIC surveillance**, where the model flags the highest *unmeasured* resistance.
""")

# ============================================================= 8. ABOUT
with tabs[8]:
    st.markdown("""
### What's in this dashboard
- **Overview** — the pitch plus the headline gap quadrant.
- **Gap Radar** — Resistance Need vs R&D Attention scatter, the priority-gap table, the all-pathogen
  indicator heatmap, and a CSV download of the full results.
- **Pathogen explorer** — pick any pathogen + antibiotic to see its six-indicator *radar profile*, a
  world resistance map, the time trend, and the bootstrap 95% CI + logistic trend statistics.
- **Trends & rigor** — the resistance-trend forest plot, PCA-justified indicator weighting, and the
  weight-sensitivity analysis (showing rankings are robust).
- **Blind-spot prediction** — the ML model's accuracy on unseen countries and its predicted
  resistance for under-surveilled regions, with a CSV download.
- **Surveillance coverage** — a world map of data volume highlighting blind spots.
- **Methodology** & **Findings & implications** — how the framework works and what it means.

### Data governance
Raw Vivli data is **restricted (Data Use Agreement)** and is **not** published — only aggregate
results (minimum cell count **n ≥ 30**). Isolate-level data stays local.

### Limitations
- Only ATLAS provides native S/I/R; others use a curated, simplified CLSI/TB breakpoint table.
- R&D Hub tags are mostly Gram-class + funder geography → attention attributed at Gram-class level.
- Continents / income tiers in the prediction model are coarse, documented approximations.
- Surveillance over-represents high-income hospital settings.

### Acknowledgement
> This publication or presentation is based on research using data from GSK, Innoviva Specialty
> Therapeutics, Johnson & Johnson, Paratek, Pfizer, Shionogi, Venatorx, Venus Remedies Limited,
> obtained through https://amr.vivli.org
""")

st.divider()
st.caption("🛰️ Treatment Gap Radar · Vivli AMR Surveillance Open Data Re-use Data Challenge "
           "(AMR ID 00013367) · data via https://amr.vivli.org · aggregate results only (n ≥ 30)")
