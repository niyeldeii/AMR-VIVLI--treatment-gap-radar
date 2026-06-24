"""ML-C: predict resistance for surveillance blind spots.

Trains a gradient-boosted model to predict %R from GENERALIZABLE features
(pathogen, drug, drug class, Gram, continent, income tier, year) -- deliberately
NOT country identity -- so it can estimate resistance in countries/regions with no
surveillance. Validated by holding out whole countries (GroupKFold): performance on
unseen countries is the blind-spot-estimation capability.
"""
import numpy as np
import pandas as pd

from .paths import PROCESSED_DIR
from .canon import _load as _load_yaml
from .geo import CONTINENT, INCOME

FEATURES = ["pathogen", "antibiotic", "drug_class", "gram", "continent", "income", "year"]


def _gram_map():
    return {p: m["gram"] for p, m in _load_yaml("pathogens.yaml").items()}


def build_frame():
    df = pd.read_parquet(PROCESSED_DIR / "combo_resistance.parquet")
    df = df.dropna(subset=["country_iso3"]).copy()
    df["continent"] = df["country_iso3"].map(CONTINENT)
    df["income"] = df["country_iso3"].map(INCOME)
    df["gram"] = df["pathogen"].map(_gram_map())
    df = df.dropna(subset=["continent", "income", "gram"])
    df["pctR"] = df["nR"] / df["n"]
    return df


def _design(df):
    X = pd.get_dummies(df[["pathogen", "antibiotic", "drug_class", "gram", "continent", "income"]],
                       dtype=float)
    X["year"] = df["year"].astype(float)
    return X


def evaluate():
    """Hold out whole countries; report weighted MAE & R2 on unseen countries."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import mean_absolute_error, r2_score
    df = build_frame().reset_index(drop=True)
    X, y, wt, grp = _design(df), df["pctR"].to_numpy(), df["n"].to_numpy(), df["country_iso3"]
    Xv = X.to_numpy()
    maes, r2s = [], []
    gkf = GroupKFold(n_splits=5)
    for tr, te in gkf.split(Xv, y, groups=grp):
        m = HistGradientBoostingRegressor(max_depth=4, learning_rate=0.08,
                                          max_iter=400, l2_regularization=1.0)
        m.fit(Xv[tr], y[tr], sample_weight=wt[tr])
        p = np.clip(m.predict(Xv[te]), 0, 1)
        maes.append(mean_absolute_error(y[te], p, sample_weight=wt[te]))
        r2s.append(r2_score(y[te], p, sample_weight=wt[te]))
    # baseline: predict global mean %R
    base = mean_absolute_error(y, np.full_like(y, np.average(y, weights=wt)), sample_weight=wt)
    return {"cv_weighted_MAE": round(float(np.mean(maes)), 4),
            "cv_R2_unseen_countries": round(float(np.mean(r2s)), 3),
            "baseline_MAE_global_mean": round(float(base), 4),
            "n_country_combos": len(df), "n_countries": int(grp.nunique())}


def predict_blindspots(min_obs=3):
    """Fit on all data; estimate %R for (continent x pathogen x drug) cells that are
    poorly observed (few contributing countries). Flags high predicted resistance where
    surveillance is thin."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    df = build_frame().reset_index(drop=True)
    X, y, wt = _design(df), df["pctR"].to_numpy(), df["n"].to_numpy()
    model = HistGradientBoostingRegressor(max_depth=4, learning_rate=0.08,
                                          max_iter=400, l2_regularization=1.0)
    model.fit(X.to_numpy(), y, sample_weight=wt)

    gram = _gram_map()
    cont_income = {c: INCOME[k] for k in INCOME for c in [CONTINENT[k]]}  # representative tier/continent
    # representative income per continent (most common)
    rep_income = (pd.DataFrame({"c": [CONTINENT[k] for k in INCOME],
                                "i": [INCOME[k] for k in INCOME]})
                  .groupby("c")["i"].agg(lambda s: s.mode().iloc[0]).to_dict())

    combos = df.groupby(["pathogen", "antibiotic", "drug_class", "continent"], observed=True).agg(
        n_countries=("country_iso3", "nunique"), obs_pctR=("pctR", "mean"),
        obs_n=("n", "sum")).reset_index()
    combos["gram"] = combos["pathogen"].map(gram)
    combos["income"] = combos["continent"].map(rep_income)
    combos["year"] = int(df["year"].max())
    Xc = pd.get_dummies(combos[["pathogen", "antibiotic", "drug_class", "gram", "continent", "income"]],
                        dtype=float)
    Xc["year"] = combos["year"].astype(float)
    Xc = Xc.reindex(columns=X.columns, fill_value=0.0)
    combos["pred_pctR"] = np.clip(model.predict(Xc.to_numpy()), 0, 1)
    combos["thin_surveillance"] = combos["n_countries"] < min_obs
    combos = combos.sort_values("pred_pctR", ascending=False)
    combos.to_parquet(PROCESSED_DIR / "blindspot_predictions.parquet")
    return combos


if __name__ == "__main__":
    print("Cross-validated (held-out countries):")
    for k, v in evaluate().items():
        print(f"  {k}: {v}")
    bs = predict_blindspots()
    thin = bs[bs["thin_surveillance"]]
    print(f"\nThin-surveillance cells: {len(thin)} of {len(bs)}")
    print("\nTop predicted-resistance blind spots (few/no surveyed countries):")
    print(thin.head(12)[["pathogen", "antibiotic", "continent", "n_countries", "pred_pctR"]]
          .round(3).to_string(index=False))
