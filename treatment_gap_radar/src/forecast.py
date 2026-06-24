"""Forecasting / early-warning: project resistance trajectories forward and estimate
when a pathogen-drug crosses a danger threshold (e.g. 50% resistant).

Method: logistic regression of resistance on calendar year (isolate level). The fitted
log-odds are linear in year, so we project P(resistant) forward and solve for the
threshold-crossing year, propagating the slope's 95% CI into a CI on that year.
"""
import numpy as np
import pandas as pd

from .paths import PROCESSED_DIR
from .indicators import load_long

THRESHOLD = 0.50
HORIZON = 2035


def _logit(p):
    return np.log(p / (1 - p))


def forecast_pair(pathogen, drug, threshold=THRESHOLD, horizon=HORIZON, min_n=300):
    import statsmodels.api as sm
    df = load_long(cols=["pathogen", "antibiotic", "sir", "year"])
    s = df[(df.pathogen == pathogen) & (df.antibiotic == drug)].dropna(subset=["sir", "year"])
    if len(s) < min_n or s["year"].nunique() < 5:
        return None
    y = (s["sir"] == "R").astype(int).to_numpy()
    if y.mean() in (0.0, 1.0):
        return None
    yr0 = int(s["year"].min())
    X = sm.add_constant((s["year"] - yr0).to_numpy().astype(float))
    res = sm.Logit(y, X).fit(disp=0)
    a, b = res.params
    (b_lo, b_hi) = res.conf_int()[1]
    last = int(s["year"].max())

    def p_at(year, slope=b, intercept=a):
        return float(1 / (1 + np.exp(-(intercept + slope * (year - yr0)))))

    def cross_year(slope):
        if slope <= 0:
            return np.nan
        t = (_logit(threshold) - a) / slope
        return yr0 + t

    pct_last = 100 * y[s["year"] == last].mean() if (s["year"] == last).any() else p_at(last) * 100
    already = p_at(last) >= threshold
    cross = np.nan if already else cross_year(b)
    cross_ci = ([np.nan, np.nan] if already else sorted([cross_year(b_lo), cross_year(b_hi)]))
    # only report future crossings within a sensible horizon
    if np.isfinite(cross) and (cross < last or cross > 2100):
        cross, cross_ci = np.nan, [np.nan, np.nan]
    return {
        "pathogen": pathogen, "drug": drug, "n": int(len(s)),
        "last_year": last, "pctR_last": round(pct_last, 1),
        "OR_per_year": round(float(np.exp(b)), 3),
        "pctR_2030": round(p_at(2030) * 100, 1),
        "pctR_2035": round(p_at(2035) * 100, 1),
        "already_over_50": bool(already),
        "cross50_year": (round(cross, 1) if np.isfinite(cross) else np.nan),
        "cross50_lo": (round(cross_ci[0], 1) if np.isfinite(cross_ci[0]) else np.nan),
        "cross50_hi": (round(cross_ci[1], 1) if np.isfinite(cross_ci[1]) else np.nan),
        "rising": bool(b > 0),
    }


def forecast_table(pairs, save=True):
    rows = [r for r in (forecast_pair(p, d) for p, d in pairs) if r]
    out = pd.DataFrame(rows).sort_values("pctR_2035", ascending=False)
    if save and len(out):
        out.to_parquet(PROCESSED_DIR / "forecasts.parquet")
        build_curves(pairs)            # precompute observed + projected curves for plotting
    return out


def build_curves(pairs, save=True):
    """Precompute observed + fitted/projected %R curves (so the dashboard needs no model)."""
    frames = []
    for p, d in pairs:
        try:
            obs, fit = fitted_curve(p, d)
        except Exception:
            continue
        obs = obs.assign(pathogen=p, drug=d, kind="observed")[["pathogen", "drug", "year", "pctR", "kind"]]
        fit = fit.assign(pathogen=p, drug=d, kind="projected")[["pathogen", "drug", "year", "pctR", "kind"]]
        frames += [obs, fit]
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if save and len(out):
        out.to_parquet(PROCESSED_DIR / "forecast_curves.parquet")
    return out


def fitted_curve(pathogen, drug, horizon=HORIZON):
    """Observed yearly %R + fitted logistic projection (for plotting)."""
    import statsmodels.api as sm
    df = load_long(cols=["pathogen", "antibiotic", "sir", "year"])
    s = df[(df.pathogen == pathogen) & (df.antibiotic == drug)].dropna(subset=["sir", "year"])
    obs = (s.assign(R=(s.sir == "R").astype(int)).groupby("year")["R"]
           .agg(["mean", "size"]).reset_index().rename(columns={"mean": "pctR", "size": "n"}))
    yr0 = int(s["year"].min())
    y = (s["sir"] == "R").astype(int).to_numpy()
    X = sm.add_constant((s["year"] - yr0).to_numpy().astype(float))
    res = sm.Logit(y, X).fit(disp=0)
    fut = np.arange(yr0, horizon + 1)
    pred = 1 / (1 + np.exp(-(res.params[0] + res.params[1] * (fut - yr0))))
    return obs, pd.DataFrame({"year": fut, "pctR": pred})


if __name__ == "__main__":
    from .pipeline import KEY_PAIRS
    t = forecast_table(KEY_PAIRS)
    cols = ["pathogen", "drug", "pctR_last", "OR_per_year", "pctR_2030", "pctR_2035",
            "cross50_year", "cross50_lo", "cross50_hi"]
    print(t[cols].to_string(index=False))
