"""Plotly figures for the Treatment Gap Radar, shared by the notebook and dashboard."""
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .paths import PROCESSED_DIR

QUAD_COLORS = {
    "PRIORITY GAP (high need / low attention)": "#d62728",
    "Well-served (high need / high attention)": "#2ca02c",
    "Possible over-investment (low need / high attention)": "#1f77b4",
    "Low priority (low need / low attention)": "#7f7f7f",
}


def load(name):
    return pd.read_parquet(PROCESSED_DIR / name)


# --------------------------------------------------------------- 1. gap quadrant
def gap_quadrant(gap=None):
    gap = gap if gap is not None else load("gap.parquet")
    g = gap.reset_index()
    rni_mid, rai_mid = g["RNI"].median(), g["RAI"].median()
    fig = px.scatter(
        g, x="RNI", y="RAI", size="n_isolates", color="quadrant",
        color_discrete_map=QUAD_COLORS, text="pathogen", size_max=55,
        hover_data={"who": True, "gap_score": ":.2f", "n_isolates": ":,"},
        labels={"RNI": "Resistance Need Index", "RAI": "R&D Attention Index"})
    fig.update_traces(textposition="top center", textfont_size=9)
    fig.add_vline(x=rni_mid, line_dash="dot", line_color="grey")
    fig.add_hline(y=rai_mid, line_dash="dot", line_color="grey")
    fig.add_annotation(x=g["RNI"].max(), y=rai_mid, text="← PRIORITY GAPS (high need, low R&D)",
                       showarrow=False, yshift=-12, font=dict(color="#d62728", size=11))
    fig.update_layout(title="Treatment Gap Radar — Resistance Need vs R&D Attention",
                      height=650, legend=dict(orientation="h", y=-0.18))
    return fig


# --------------------------------------------------------- 2. indicator heatmap
def indicator_heatmap(rni=None):
    rni = rni if rni is not None else load("rni.parquet")
    cols = ["prevalence_n", "mic_drift_n", "mdr_n", "geo_spread_n", "scarcity_n", "pediatric_n"]
    nice = ["Prevalence", "MIC drift", "MDR", "Geo spread", "Scarcity", "Pediatric"]
    d = rni.sort_values("RNI", ascending=True)
    fig = go.Figure(go.Heatmap(
        z=d[cols].values, x=nice, y=d.index.tolist(),
        colorscale="Reds", zmin=0, zmax=1, colorbar=dict(title="norm.")))
    fig.update_layout(title="Resistance indicators by pathogen (normalized 0–1)",
                      height=max(400, 22 * len(d)))
    return fig


# --------------------------------------------------------- 3. resistance choropleth
def resistance_choropleth(pathogen, drug=None, min_n=20):
    pc = load("combo_resistance.parquet")
    sub = pc[pc["pathogen"] == pathogen]
    if drug:
        sub = sub[sub["antibiotic"] == drug]
    agg = sub.groupby("country_iso3", as_index=False).agg(nR=("nR", "sum"), n=("n", "sum"))
    agg["pctR"] = 100 * agg["nR"] / agg["n"]
    agg = agg[agg["n"] >= min_n]
    title = f"{pathogen}" + (f" — {drug}" if drug else "") + " : % resistant by country"
    fig = px.choropleth(agg, locations="country_iso3", color="pctR",
                        color_continuous_scale="Reds", range_color=(0, 100),
                        hover_data={"n": ":,"}, labels={"pctR": "% R"})
    fig.update_layout(title=title, height=500, geo=dict(showframe=False))
    return fig


# --------------------------------------------------------- 4. MIC-drift trend
def resistance_trend(pathogen, drug):
    pc = load("combo_resistance.parquet")
    sub = pc[(pc["pathogen"] == pathogen) & (pc["antibiotic"] == drug)]
    ts = sub.groupby("year", as_index=False).agg(nR=("nR", "sum"), n=("n", "sum"))
    ts["pctR"] = 100 * ts["nR"] / ts["n"]
    fig = px.line(ts, x="year", y="pctR", markers=True,
                  labels={"pctR": "% resistant", "year": "Year"})
    fig.update_layout(title=f"{pathogen} — {drug}: resistance trend", height=420,
                      yaxis_range=[0, max(100, ts["pctR"].max() * 1.1 if len(ts) else 100)])
    return fig


# --------------------------------------------------- 5. surveillance coverage (blind spots)
def coverage_map():
    pc = load("combo_resistance.parquet")
    cov = (pc.groupby("country_iso3").agg(records=("n", "sum"),
           pathogens=("pathogen", "nunique")).reset_index())
    cov["log_records"] = np.log10(cov["records"] + 1)
    fig = px.choropleth(cov, locations="country_iso3", color="log_records",
                        color_continuous_scale="Blues",
                        hover_data={"records": ":,", "pathogens": True, "log_records": False},
                        labels={"log_records": "log10(records)"})
    fig.update_layout(title="Surveillance coverage (data volume by country — blind spots in white)",
                      height=500, geo=dict(showframe=False))
    return fig
