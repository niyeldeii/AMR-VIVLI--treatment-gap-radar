"""Precompute rigor + ML artifacts to parquet so the dashboard/notebook just read them
(no scikit-learn / statsmodels needed at serve time).

Outputs (data_processed/):
  trend_models.parquet      logistic resistance~year OR/yr + 95% CI + p
  ci_keypairs.parquet       bootstrap 95% CI of %R for key pathogen-drug pairs
  sensitivity.parquet       RNI rank stability across weight schemes
  rigor_summary.parquet     PCA weights + explained variance + blind-spot model metrics
  blindspot_predictions.parquet   predicted %R by continent x pathogen x drug
"""
import json
import pandas as pd

from .paths import PROCESSED_DIR
from . import models, predict

KEY_PAIRS = [
    ("Acinetobacter baumannii", "Meropenem"),
    ("Klebsiella pneumoniae", "Meropenem"),
    ("Klebsiella pneumoniae", "Ceftriaxone"),
    ("Escherichia coli", "Ceftriaxone"),
    ("Escherichia coli", "Ciprofloxacin"),
    ("Pseudomonas aeruginosa", "Meropenem"),
    ("Enterococcus faecium", "Vancomycin"),
    ("Staphylococcus aureus", "Oxacillin"),
    ("Streptococcus pneumoniae", "Penicillin"),
    ("Escherichia coli", "Trimethoprim-sulfamethoxazole"),
]


def main():
    print("trend models ..."); trends = models.trend_table(KEY_PAIRS)
    print("bootstrap CIs ...")
    ci = pd.DataFrame([{"pathogen": p, "drug": d, **models.bootstrap_pctR(p, d)} for p, d in KEY_PAIRS])
    ci.to_parquet(PROCESSED_DIR / "ci_keypairs.parquet")
    print("sensitivity ..."); sens = models.sensitivity(); sens.to_parquet(PROCESSED_DIR / "sensitivity.parquet")
    print("pca ..."); pca = models.pca_weights()
    print("blind-spot model CV ..."); metrics = predict.evaluate()
    print("blind-spot predictions ..."); predict.predict_blindspots()

    summary = {**{f"w_{k}": v for k, v in pca["weights"].items()},
               "pc1_var": pca["explained_variance"][0], **metrics}
    pd.DataFrame([summary]).to_parquet(PROCESSED_DIR / "rigor_summary.parquet")
    (PROCESSED_DIR / "rigor_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print("done. metrics:", metrics)
    return summary


if __name__ == "__main__":
    main()
