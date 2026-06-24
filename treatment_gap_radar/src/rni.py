"""Resistance Need Index: normalize the six indicators and combine via weights."""
import numpy as np
import pandas as pd
import yaml

from .paths import PROCESSED_DIR, CONFIG_DIR

IND_COLS = ["prevalence", "mic_drift", "mdr", "geo_spread", "scarcity", "pediatric"]
MIN_ISOLATES = 200   # pathogens below this are excluded from the ranked index (too sparse)


def _weights():
    with open(CONFIG_DIR / "weights.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)["rni"]


def minmax(s):
    s = s.astype(float)
    lo, hi = s.min(), s.max()
    if not np.isfinite(lo) or hi == lo:
        return pd.Series(0.5, index=s.index)
    return (s - lo) / (hi - lo)


def compute_rni(weights=None, save=True):
    ind = pd.read_parquet(PROCESSED_DIR / "indicators.parquet")
    ind = ind[ind["n_isolates"] >= MIN_ISOLATES].copy()
    w = weights or _weights()

    norm = pd.DataFrame(index=ind.index)
    for c in IND_COLS:
        col = ind[c].copy()
        # missing => no evidence of need => fill with column min before scaling
        col = col.fillna(col.min())
        norm[c + "_n"] = minmax(col)

    wsum = sum(w[c] for c in IND_COLS)
    ind["RNI"] = sum(norm[c + "_n"] * w[c] for c in IND_COLS) / wsum
    out = pd.concat([ind, norm], axis=1).sort_values("RNI", ascending=False)
    if save:
        out.to_parquet(PROCESSED_DIR / "rni.parquet")
    return out


if __name__ == "__main__":
    r = compute_rni()
    print(r[["who", "n_isolates", "RNI"] + IND_COLS].round(3).to_string())
