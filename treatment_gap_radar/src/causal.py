"""Causal layer: is R&D attention associated with how fast a pathogen's resistance grows?

Design (observational, causal-inference-style): a pathogen x year panel of resistance, with a
year x R&D-level interaction and pathogen fixed effects. Because a pathogen's R&D-attention level
is time-invariant it is absorbed by the pathogen fixed effects; the *interaction* (year x R&D level)
identifies whether resistance grows at different rates for high- vs low-R&D pathogens, net of each
pathogen's baseline. We also give a transparent counterfactual scenario.

This is association, not proof of causation (R&D is largely reactive, allocated after resistance
emerges, and confounded). We state that plainly; the value is a defensible, policy-relevant estimate
of whether current R&D attention is keeping pace with — or chasing — resistance.
"""
import numpy as np
import pandas as pd

from .paths import PROCESSED_DIR

MIN_YEARS = 6
MIN_N_PER_CELL = 50


def _panel():
    """pathogen x year weighted %R, joined to pathogen R&D-attention tertile."""
    combo = pd.read_parquet(PROCESSED_DIR / "combo_resistance.parquet")
    pan = (combo.groupby(["pathogen", "year"], observed=True)
           .agg(nR=("nR", "sum"), n=("n", "sum")).reset_index())
    pan = pan[pan["n"] >= MIN_N_PER_CELL]
    pan["pctR"] = 100 * pan["nR"] / pan["n"]
    keep = pan.groupby("pathogen")["year"].nunique()
    pan = pan[pan["pathogen"].isin(keep[keep >= MIN_YEARS].index)]

    rai = pd.read_parquet(PROCESSED_DIR / "rai.parquet")[["RAI"]]
    rai["rd_level"] = pd.qcut(rai["RAI"], 3, labels=["Low", "Mid", "High"])
    pan = pan.merge(rai, left_on="pathogen", right_index=True, how="inner")
    pan["year_c"] = pan["year"] - pan["year"].mean()
    return pan


def growth_by_rd_level(save=True):
    """Resistance growth rate (%R per year) by R&D-attention tertile, pathogen FE, weighted."""
    import statsmodels.formula.api as smf
    pan = _panel()
    # year_c main = growth for reference level (High); year_c:rd_level = differential growth.
    m = smf.wls("pctR ~ year_c + year_c:C(rd_level) + C(pathogen)",
                data=pan, weights=pan["n"]).fit()
    base = m.params.get("year_c", np.nan)          # High (reference)
    rows = []
    for lvl in ["Low", "Mid", "High"]:
        key = f"year_c:C(rd_level)[T.{lvl}]"
        slope = base + (m.params[key] if key in m.params else 0.0)
        p = (m.pvalues[key] if key in m.pvalues else np.nan)
        rows.append({"rd_level": lvl, "growth_pctR_per_year": round(float(slope), 3),
                     "diff_vs_High_p": (round(float(p), 4) if not np.isnan(p) else np.nan),
                     "n_pathogens": int(pan[pan.rd_level == lvl]["pathogen"].nunique())})
    out = pd.DataFrame(rows)
    if save:
        out.to_parquet(PROCESSED_DIR / "causal_growth.parquet")
    return out, m


def counterfactual(save=True):
    """Transparent what-if: if every pathogen's resistance grew at the slowest observed
    (best-case) rate, how much 2035 resistance would be averted vs the status-quo forecast.
    Illustrative scenario, not a causal guarantee."""
    growth, _ = growth_by_rd_level(save=False)
    best = growth["growth_pctR_per_year"].min()
    fc = pd.read_parquet(PROCESSED_DIR / "forecasts.parquet")
    cf = fc[["pathogen", "drug", "pctR_last", "last_year", "OR_per_year", "pctR_2035"]].copy()
    yrs = 2035 - cf["last_year"]
    cf["pctR_2035_bestcase"] = (cf["pctR_last"] + best * yrs).clip(0, 100).round(1)
    cf["averted_2035_pts"] = (cf["pctR_2035"] - cf["pctR_2035_bestcase"]).round(1)
    cf = cf.sort_values("averted_2035_pts", ascending=False)
    if save:
        cf.to_parquet(PROCESSED_DIR / "causal_counterfactual.parquet")
    return cf


if __name__ == "__main__":
    g, m = growth_by_rd_level()
    print("=== Resistance growth (%R/yr) by R&D-attention tertile (pathogen FE, weighted) ===")
    print(g.to_string(index=False))
    print("\nInterpretation: if HIGH-R&D pathogens grow *faster*, R&D is reactive — chasing"
          " resistance rather than pre-empting it.\n")
    print("=== Counterfactual: 2035 resistance averted under best-case growth ===")
    print(counterfactual()[["pathogen", "drug", "pctR_2035", "pctR_2035_bestcase",
                            "averted_2035_pts"]].head(8).to_string(index=False))
