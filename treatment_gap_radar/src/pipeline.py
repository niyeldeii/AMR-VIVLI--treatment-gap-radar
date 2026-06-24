"""Single end-to-end pipeline: raw data -> everything the project needs.

    python -m src.pipeline                 # full rebuild (harmonize + analysis + rigor + ML)
    python -m src.pipeline --no-harmonize  # reuse existing isolates_long.parquet
    python -m src.pipeline --fast          # skip the slower ML blind-spot model

Stages: harmonize -> indicators -> RNI/RAI/gap -> rigor (CIs, trends, PCA, sensitivity)
-> ML blind-spot prediction.  All outputs land in data_processed/ as parquet.
"""
import sys

import pandas as pd

from . import harmonize, indicators, gap, models, predict, rai, forecast, validate, causal
from .paths import PROCESSED_DIR

# pathogen-drug pairs used for trend models + bootstrap CIs
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


def _rigor(with_ml=True):
    models.trend_table(KEY_PAIRS)
    pd.DataFrame([{"pathogen": p, "drug": d, **models.bootstrap_pctR(p, d)}
                 for p, d in KEY_PAIRS]).to_parquet(PROCESSED_DIR / "ci_keypairs.parquet")
    models.sensitivity().to_parquet(PROCESSED_DIR / "sensitivity.parquet")
    pca = models.pca_weights()
    rai.attribution_robustness()                 # gap stable across R&D attribution schemes
    forecast.forecast_table(KEY_PAIRS)           # early-warning: project trajectories to threshold
    validate.validate()                          # external benchmark vs WHO GLASS
    causal.growth_by_rd_level()                  # is R&D reactive? growth by attention level
    causal.counterfactual()                      # illustrative averted-resistance scenario
    metrics = {}
    if with_ml:
        metrics = predict.evaluate()
        predict.predict_blindspots()
    summary = {**{f"w_{k}": v for k, v in pca["weights"].items()},
               "pc1_var": pca["explained_variance"][0], **metrics}
    pd.DataFrame([summary]).to_parquet(PROCESSED_DIR / "rigor_summary.parquet")
    return summary


def main(skip_harmonize=False, with_ml=True):
    if not skip_harmonize:
        print("[1/4] Harmonizing all datasets ..."); harmonize.build()
    print("[2/4] Computing resistance indicators ..."); indicators.compute_indicators()
    print("[3/4] Computing RNI / RAI / gap ...")
    g = gap.compute_gap()
    print("[4/4] Rigor + ML (CIs, trends, PCA, sensitivity, blind-spot model) ...")
    summary = _rigor(with_ml=with_ml)
    n_gap = int((g["quadrant"].str.startswith("PRIORITY")).sum())
    print(f"\nDone. {len(g)} pathogens scored; {n_gap} priority gaps. "
          f"PCA PC1 {summary.get('pc1_var', float('nan')):.2f}; "
          f"blind-spot R2 {summary.get('cv_R2_unseen_countries', float('nan'))}.")
    return g


if __name__ == "__main__":
    main(skip_harmonize="--no-harmonize" in sys.argv, with_ml="--fast" not in sys.argv)
