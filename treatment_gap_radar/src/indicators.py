"""Compute the six resistance indicators per pathogen, plus supporting tables."""
import numpy as np
import pandas as pd
import yaml

from .paths import PROCESSED_DIR, CONFIG_DIR

LONG = PROCESSED_DIR / "isolates_long.parquet"


def _cfg():
    with open(CONFIG_DIR / "weights.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _slope(x, y):
    """OLS slope of y on x; nan if <3 points or no variance."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3 or x.std() == 0:
        return np.nan
    return np.polyfit(x, y, 1)[0]


def load_long(cols=None):
    return pd.read_parquet(LONG, columns=cols)


def combo_resistance(df, min_n):
    """%R per (pathogen, antibiotic, country, year) with >= min_n tested."""
    r = df.dropna(subset=["sir", "pathogen"]).copy()
    r["is_R"] = (r["sir"] == "R").astype(int)
    g = (r.groupby(["pathogen", "antibiotic", "drug_class", "country_iso3", "year"], observed=True)
           .agg(n=("is_R", "size"), nR=("is_R", "sum")).reset_index())
    g = g[g["n"] >= min_n]
    g["pctR"] = g["nR"] / g["n"]
    return g


def compute_indicators():
    cfg = _cfg()
    min_n = cfg["min_tested_combo"]
    geo_thr = cfg["country_resistant_threshold"]
    comp_thr = cfg["class_compromised_threshold"]

    df = load_long(cols=["isolate_id", "pathogen", "gram", "who", "country_iso3",
                         "year", "pediatric", "antibiotic", "drug_class", "sir"])
    df = df.dropna(subset=["pathogen"])
    df["is_R"] = (df["sir"] == "R").map({True: 1, False: 0})  # NA stays NA

    # ---- 1. prevalence: overall %R per pathogen
    prev = (df.dropna(subset=["sir"]).assign(R=lambda d: (d.sir == "R").astype(int))
              .groupby("pathogen", observed=True)["R"].mean().rename("prevalence"))

    # ---- 2. mic_drift: trend of yearly %R per pathogen (avg slope across drugs)
    yr = combo_resistance(df, min_n)  # reuse, but we want pathogen-year over drugs
    py = (df.dropna(subset=["sir"]).assign(R=lambda d: (d.sir == "R").astype(int))
            .groupby(["pathogen", "drug_class", "year"], observed=True)["R"]
            .agg(["mean", "size"]).reset_index())
    py = py[py["size"] >= min_n]
    drift = (py.groupby(["pathogen", "drug_class"], observed=True)
               .apply(lambda d: _slope(d["year"], d["mean"]), include_groups=False)
               .groupby("pathogen").mean().rename("mic_drift"))

    # ---- 3. MDR: % isolates resistant to >=3 drug classes (>=3 classes tested)
    iso = df.dropna(subset=["sir"]).copy()
    grp = iso.groupby("isolate_id", observed=True)
    classes_tested = grp["drug_class"].nunique()
    classes_R = (iso[iso.sir == "R"].groupby("isolate_id", observed=True)["drug_class"]
                 .nunique().reindex(classes_tested.index).fillna(0))
    patho_of = grp["pathogen"].first()
    iso_tbl = pd.DataFrame({"pathogen": patho_of, "ntested": classes_tested, "nR": classes_R})
    iso_tbl = iso_tbl[iso_tbl["ntested"] >= 3]
    iso_tbl["mdr"] = (iso_tbl["nR"] >= 3).astype(int)
    mdr = iso_tbl.groupby("pathogen", observed=True)["mdr"].mean().rename("mdr")

    # ---- 4. geo_spread: share of surveyed countries with pathogen %R > threshold
    cc = (df.dropna(subset=["sir"]).assign(R=lambda d: (d.sir == "R").astype(int))
            .groupby(["pathogen", "country_iso3"], observed=True)["R"]
            .agg(["mean", "size"]).reset_index())
    cc = cc[cc["size"] >= min_n]
    geo = (cc.assign(hot=lambda d: (d["mean"] > geo_thr).astype(int))
             .groupby("pathogen", observed=True)["hot"].mean().rename("geo_spread"))

    # ---- 5. scarcity: share of a pathogen's drug classes that are compromised (%R > threshold)
    dcls = (df.dropna(subset=["sir"]).assign(R=lambda d: (d.sir == "R").astype(int))
              .groupby(["pathogen", "drug_class"], observed=True)["R"]
              .agg(["mean", "size"]).reset_index())
    dcls = dcls[dcls["size"] >= min_n]
    scarcity = (dcls.assign(bad=lambda d: (d["mean"] > comp_thr).astype(int))
                  .groupby("pathogen", observed=True)["bad"].mean().rename("scarcity"))

    # ---- 6. pediatric: %R among pediatric isolates
    ped = (df[(df.pediatric == True)].dropna(subset=["sir"])
             .assign(R=lambda d: (d.sir == "R").astype(int))
             .groupby("pathogen", observed=True)["R"].mean().rename("pediatric"))

    meta = (df.groupby("pathogen", observed=True)[["gram", "who"]].first())
    n_iso = df.groupby("pathogen", observed=True)["isolate_id"].nunique().rename("n_isolates")

    ind = pd.concat([meta, n_iso, prev, drift, mdr, geo, scarcity, ped], axis=1)
    ind.to_parquet(PROCESSED_DIR / "indicators.parquet")
    yr.to_parquet(PROCESSED_DIR / "combo_resistance.parquet")
    cc.to_parquet(PROCESSED_DIR / "pathogen_country_resistance.parquet")
    return ind


if __name__ == "__main__":
    ind = compute_indicators()
    cols = ["who", "n_isolates", "prevalence", "mic_drift", "mdr", "geo_spread", "scarcity", "pediatric"]
    print(ind.sort_values("prevalence", ascending=False)[cols].round(3).to_string())
