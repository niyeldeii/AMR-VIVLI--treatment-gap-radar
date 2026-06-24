"""Statistical rigor + light ML for defensibility.

A. Trend models     - logistic regression of resistance ~ year (OR/year, p-value, 95% CI).
B. PCA weighting     - data-driven indicator weights vs the equal-weight default.
   + bootstrap CIs for %R, and a weight sensitivity analysis (rank stability).
"""
import numpy as np
import pandas as pd

from .paths import PROCESSED_DIR
from .indicators import load_long
from .rni import compute_rni, IND_COLS


# ----------------------------------------------------------- bootstrap CI for %R
def bootstrap_pctR(pathogen, drug, n_boot=2000, seed=0):
    df = load_long(cols=["pathogen", "antibiotic", "sir"])
    s = df[(df.pathogen == pathogen) & (df.antibiotic == drug)]["sir"].dropna()
    y = (s == "R").to_numpy().astype(float)
    if len(y) < 10:
        return {"pctR": np.nan, "lo": np.nan, "hi": np.nan, "n": len(y)}
    rng = np.random.default_rng(seed)
    boots = [rng.choice(y, size=len(y), replace=True).mean() for _ in range(n_boot)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"pctR": 100 * y.mean(), "lo": 100 * lo, "hi": 100 * hi, "n": len(y)}


# ------------------------------------------------------- A. logistic trend model
def trend_model(pathogen, drug, min_n=200):
    """Logistic regression: P(resistant) ~ year. Returns OR/year, p-value, 95% CI."""
    import statsmodels.api as sm
    df = load_long(cols=["pathogen", "antibiotic", "sir", "year"])
    s = df[(df.pathogen == pathogen) & (df.antibiotic == drug)].dropna(subset=["sir", "year"])
    if len(s) < min_n or s["year"].nunique() < 4:
        return None
    y = (s["sir"] == "R").astype(int).to_numpy()
    if y.sum() == 0 or y.sum() == len(y):
        return None
    X = sm.add_constant((s["year"] - s["year"].min()).to_numpy().astype(float))
    res = sm.Logit(y, X).fit(disp=0)
    coef = res.params[1]
    ci = res.conf_int()[1]
    return {"pathogen": pathogen, "drug": drug, "n": int(len(s)),
            "OR_per_year": float(np.exp(coef)),
            "ci_lo": float(np.exp(ci[0])), "ci_hi": float(np.exp(ci[1])),
            "p_value": float(res.pvalues[1])}


def trend_table(pairs):
    rows = [r for r in (trend_model(p, d) for p, d in pairs) if r]
    out = pd.DataFrame(rows)
    if len(out):
        out.to_parquet(PROCESSED_DIR / "trend_models.parquet")
    return out


# ------------------------------------------------ B. PCA-derived indicator weights
def pca_weights():
    """Data-driven weights from PC1 loadings on the six normalized indicators."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    rni = compute_rni(save=False)
    cols = [c + "_n" for c in IND_COLS]
    X = StandardScaler().fit_transform(rni[cols].fillna(rni[cols].mean()).to_numpy())
    pca = PCA(n_components=3).fit(X)
    pc1 = np.abs(pca.components_[0])
    weights = dict(zip(IND_COLS, (pc1 / pc1.sum()).round(3)))
    return {"weights": weights,
            "explained_variance": pca.explained_variance_ratio_.round(3).tolist()}


# ------------------------------------------------------- weight sensitivity check
def sensitivity():
    """Rank stability of RNI under several weighting schemes (Spearman vs equal weights)."""
    from scipy.stats import spearmanr
    schemes = {
        "equal": {c: 1 for c in IND_COLS},
        "prevalence_heavy": {**{c: 1 for c in IND_COLS}, "prevalence": 3},
        "clinical": {"prevalence": 2, "mic_drift": 1, "mdr": 2, "geo_spread": 1,
                     "scarcity": 2, "pediatric": 1},
        "pca": pca_weights()["weights"],
    }
    base = compute_rni(weights=schemes["equal"], save=False)["RNI"]
    rows = []
    for name, w in schemes.items():
        r = compute_rni(weights=w, save=False)["RNI"].reindex(base.index)
        rho = spearmanr(base.values, r.values).statistic
        rows.append({"scheme": name, "spearman_vs_equal": round(float(rho), 3),
                     "top_pathogen": r.sort_values(ascending=False).index[0]})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print("PCA weights:", pca_weights())
    print("\nSensitivity (rank stability):")
    print(sensitivity().to_string(index=False))
    print("\nExample trend (A. baumannii meropenem):", trend_model("Acinetobacter baumannii", "Meropenem"))
    print("Bootstrap CI (A. baumannii meropenem):", bootstrap_pctR("Acinetobacter baumannii", "Meropenem"))
