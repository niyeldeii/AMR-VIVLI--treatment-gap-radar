"""Treatment Gap: join RNI vs RAI, classify quadrants, score the gap."""
import pandas as pd

from .paths import PROCESSED_DIR
from .rni import compute_rni
from .rai import compute_rai


def compute_gap():
    rni = compute_rni()[["who", "n_isolates", "RNI"]]
    rai = compute_rai()[["projects", "investment", "pipeline", "RAI"]]
    g = rni.join(rai, how="left")
    g["RAI"] = g["RAI"].fillna(g["RAI"].min())

    rni_mid = g["RNI"].median()
    rai_mid = g["RAI"].median()

    def quad(r):
        hi_need = r["RNI"] >= rni_mid
        hi_att = r["RAI"] >= rai_mid
        if hi_need and not hi_att:
            return "PRIORITY GAP (high need / low attention)"
        if hi_need and hi_att:
            return "Well-served (high need / high attention)"
        if not hi_need and hi_att:
            return "Possible over-investment (low need / high attention)"
        return "Low priority (low need / low attention)"

    g["quadrant"] = g.apply(quad, axis=1)
    g["gap_score"] = g["RNI"] - g["RAI"]      # >0 => under-served relative to need
    g = g.sort_values("gap_score", ascending=False)
    g.to_parquet(PROCESSED_DIR / "gap.parquet")
    return g


if __name__ == "__main__":
    g = compute_gap()
    pd.set_option("display.width", 160)
    print(g[["who", "n_isolates", "RNI", "RAI", "gap_score", "quadrant"]].round(3).to_string())
