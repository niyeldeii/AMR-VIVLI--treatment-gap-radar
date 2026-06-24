"""End-to-end build: raw data -> long table -> indicators -> RNI/RAI -> gap.

Usage (from the treatment_gap_radar folder):
    python -m src.pipeline              # full rebuild
    python -m src.pipeline --no-harmonize   # reuse existing isolates_long.parquet
"""
import sys

from . import harmonize, indicators, gap


def main(skip_harmonize=False):
    if not skip_harmonize:
        print("[1/3] Harmonizing all datasets ...")
        harmonize.build()
    print("[2/3] Computing resistance indicators ...")
    indicators.compute_indicators()
    print("[3/3] Computing RNI / RAI / gap ...")
    g = gap.compute_gap()
    n_gap = (g["quadrant"].str.startswith("PRIORITY")).sum()
    print(f"Done. {len(g)} pathogens scored; {int(n_gap)} priority treatment gaps.")
    return g


if __name__ == "__main__":
    main(skip_harmonize="--no-harmonize" in sys.argv)
