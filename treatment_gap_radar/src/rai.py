"""R&D Attention Index from the Global AMR R&D Hub export (Projects.xlsx).

Per canonical pathogen we keyword-match the project text and aggregate:
  investment (total USD), projects (count), pipeline (# distinct named products).
Then normalize and combine. Granularity caveat: R&D records are tagged coarsely
(mostly Gram class) and by funder geography, not by where resistance occurs -- so this
is an approximate, text-derived attribution. Documented in the notebook.
"""
import re
import numpy as np
import pandas as pd
import yaml

from .paths import raw_path, PROCESSED_DIR, CONFIG_DIR
from . import canon

# keyword(s) per canonical pathogen (lowercased, regex-OR). Genus-level terms included.
KEYWORDS = {
    "Escherichia coli":           [r"escherichia coli", r"\be\.? coli\b"],
    "Klebsiella pneumoniae":      [r"klebsiella"],
    "Klebsiella oxytoca":         [r"klebsiella oxytoca"],
    "Klebsiella aerogenes":       [r"aerogenes"],
    "Acinetobacter baumannii":    [r"acinetobacter"],
    "Pseudomonas aeruginosa":     [r"pseudomonas"],
    "Enterobacter cloacae":       [r"enterobacter"],
    "Serratia marcescens":        [r"serratia"],
    "Citrobacter freundii":       [r"citrobacter"],
    "Proteus mirabilis":          [r"proteus"],
    "Providencia rettgeri":       [r"providencia"],
    "Providencia stuartii":       [r"providencia"],
    "Haemophilus influenzae":     [r"haemophilus"],
    "Neisseria gonorrhoeae":      [r"gonorrh", r"neisseria gonorrh"],
    "Salmonella enterica":        [r"salmonella"],
    "Staphylococcus aureus":      [r"staphylococcus aureus", r"\bmrsa\b", r"\bmssa\b", r"s\.? aureus"],
    "Staphylococcus epidermidis": [r"epidermidis"],
    "Streptococcus pneumoniae":   [r"pneumococc", r"streptococcus pneumoniae", r"s\.? pneumoniae"],
    "Streptococcus agalactiae":   [r"agalactiae", r"group b strep"],
    "Streptococcus pyogenes":     [r"pyogenes", r"group a strep"],
    "Enterococcus faecalis":      [r"enterococc", r"faecalis"],
    "Enterococcus faecium":       [r"enterococc", r"faecium", r"\bvre\b"],
    "Mycobacterium tuberculosis": [r"tuberculosis", r"\bm\.? tb\b", r"\bmdr-?tb\b"],
}


def _weights():
    with open(CONFIG_DIR / "weights.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)["rai"]


def minmax(s):
    s = s.astype(float)
    lo, hi = s.min(), s.max()
    if hi == lo:
        return pd.Series(0.5, index=s.index)
    return (s - lo) / (hi - lo)


def load_projects():
    df = pd.read_excel(raw_path("RND_HUB"), sheet_name="data")
    df["text"] = (df["Title"].astype(str) + " " + df["Abstract"].astype(str) + " " +
                  df["Individual Infectious Agent"].astype(str) + " " +
                  df["Infectious Agent"].astype(str) + " " + df["Disease"].astype(str) + " " +
                  df["Product Name"].astype(str)).str.lower()
    df["usd"] = pd.to_numeric(df["Amount USD"], errors="coerce")
    return df


def compute_rai(weights=None):
    df = load_projects()
    w = weights or _weights()
    rows = []
    for patho, pats in KEYWORDS.items():
        rgx = re.compile("|".join(pats))
        hit = df[df["text"].str.contains(rgx, regex=True, na=False)]
        rows.append({
            "pathogen": patho,
            "projects": len(hit),
            "investment": hit["usd"].sum(min_count=1),
            "pipeline": hit.loc[hit["Product Name"].notna(), "Product Name"].nunique(),
            "therapeutics": (hit["Research Area"].astype(str).str.contains("Therapeutics")).sum(),
            "diagnostics": (hit["Research Area"].astype(str).str.contains("Diagnostics")).sum(),
        })
    rai = pd.DataFrame(rows).set_index("pathogen")
    rai["investment"] = rai["investment"].fillna(0)

    # log-scale heavy-tailed money/counts before normalizing
    comp = {}
    for c in ["investment", "projects", "pipeline"]:
        comp[c + "_n"] = minmax(np.log1p(rai[c]))
    comp = pd.DataFrame(comp, index=rai.index)
    wsum = sum(w[c] for c in ["investment", "projects", "pipeline"])
    rai["RAI"] = sum(comp[c + "_n"] * w[c] for c in ["investment", "projects", "pipeline"]) / wsum
    out = pd.concat([rai, comp], axis=1).sort_values("RAI", ascending=False)
    out.to_parquet(PROCESSED_DIR / "rai.parquet")
    return out


if __name__ == "__main__":
    r = compute_rai()
    print(r[["projects", "investment", "pipeline", "therapeutics", "RAI"]].round(2).to_string())
