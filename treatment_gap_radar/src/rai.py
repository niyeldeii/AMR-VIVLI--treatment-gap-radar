"""R&D Attention Index from the Global AMR R&D Hub export (Projects.xlsx).

Hardened attribution. For each pathogen we count R&D that could *plausibly address it*:
  (broad-spectrum work tagged for its Gram class)  UNION  (work naming the species).
This is defensible because a broad Gram-negative antibiotic treats e.g. Providencia even when
that genus is never named -- so species are not falsely scored as "zero R&D".

Per pathogen we aggregate over that project set:
  investment (total USD), projects (count), pipeline (# distinct named products),
then optionally weight clinical work (therapeutics/diagnostics) above basic research,
restrict to the modern resistance era, log-scale, normalize, and combine.

Granularity caveat (documented): R&D Hub records are tagged mostly at Gram-class level and by
funder geography; species-level signal is text-derived. The gap is therefore interpreted at
pathogen level, with a robust Gram-class baseline carried for every species.
"""
import re
import numpy as np
import pandas as pd
import yaml

from .paths import raw_path, PROCESSED_DIR, CONFIG_DIR
from .canon import _load as _load_yaml

MIN_END_YEAR = 2010           # modern resistance era (projects ending >= this)
CLINICAL_WEIGHT = 2.0         # weight for therapeutics/diagnostics projects vs basic research

# species-specific keyword patterns (word-boundary, lowercased)
KEYWORDS = {
    "Escherichia coli":           [r"\bescherichia coli\b", r"\be\.?\s?coli\b"],
    "Klebsiella pneumoniae":      [r"\bklebsiella pneumoniae\b", r"\bk\.?\s?pneumoniae\b"],
    "Klebsiella oxytoca":         [r"\bklebsiella oxytoca\b"],
    "Klebsiella aerogenes":       [r"\baerogenes\b"],
    "Acinetobacter baumannii":    [r"\bacinetobacter\b"],
    "Pseudomonas aeruginosa":     [r"\bpseudomonas\b"],
    "Enterobacter cloacae":       [r"\benterobacter\b"],
    "Serratia marcescens":        [r"\bserratia\b"],
    "Citrobacter freundii":       [r"\bcitrobacter\b"],
    "Proteus mirabilis":          [r"\bproteus\b"],
    "Providencia rettgeri":       [r"\bprovidencia\b"],
    "Providencia stuartii":       [r"\bprovidencia\b"],
    "Haemophilus influenzae":     [r"\bhaemophilus\b"],
    "Neisseria gonorrhoeae":      [r"\bgonorrh", r"\bneisseria gonorrh"],
    "Salmonella enterica":        [r"\bsalmonella\b"],
    "Staphylococcus aureus":      [r"\bstaphylococcus aureus\b", r"\bmrsa\b", r"\bmssa\b", r"\bs\.?\s?aureus\b"],
    "Staphylococcus epidermidis": [r"\bepidermidis\b"],
    "Streptococcus pneumoniae":   [r"\bpneumococc", r"\bstreptococcus pneumoniae\b", r"\bs\.?\s?pneumoniae\b"],
    "Streptococcus agalactiae":   [r"\bagalactiae\b", r"\bgroup b strep"],
    "Streptococcus pyogenes":     [r"\bpyogenes\b", r"\bgroup a strep"],
    "Enterococcus faecalis":      [r"\bfaecalis\b", r"\benterococc"],
    "Enterococcus faecium":       [r"\bfaecium\b", r"\bvre\b", r"\benterococc"],
    "Mycobacterium tuberculosis": [r"\btuberculosis\b", r"\bmdr-?tb\b", r"\bm\.?\s?tb\b"],
}


def _weights():
    with open(CONFIG_DIR / "weights.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)["rai"]


def _gram_of():
    out = {}
    for p, meta in _load_yaml("pathogens.yaml").items():
        out[p] = meta["gram"]
    return out


def minmax(s):
    s = s.astype(float)
    lo, hi = s.min(), s.max()
    return pd.Series(0.5, index=s.index) if hi == lo else (s - lo) / (hi - lo)


def load_projects():
    df = pd.read_excel(raw_path("RND_HUB"), sheet_name="data")
    _tc = ["Title", "Abstract", "Individual Infectious Agent", "Infectious Agent",
           "Disease", "Product Name"]
    df["text"] = df[_tc[0]].astype(str)
    for c in _tc[1:]:
        df["text"] = df["text"] + " " + df[c].astype(str)
    df["text"] = df["text"].str.lower()
    df["usd"] = pd.to_numeric(df["Amount USD"], errors="coerce").fillna(0)
    df["end_year"] = pd.to_numeric(df["End Year"], errors="coerce")
    ia = df["Infectious Agent"].astype(str).str.lower()
    dis = df["Disease"].astype(str).str.lower()
    df["is_gn"] = ia.str.contains("gram negative")
    df["is_gp"] = ia.str.contains("gram positive")
    df["is_tb"] = dis.str.contains("tuberculosis") | df["text"].str.contains("tuberculosis")
    ra = df["Research Area"].astype(str)
    df["clinical"] = ra.str.contains("Therapeutics") | ra.str.contains("Diagnostics")
    df["w"] = np.where(df["clinical"], CLINICAL_WEIGHT, 1.0)
    return df


def compute_rai(weights=None, save=True):
    df = load_projects()
    df = df[df["end_year"].fillna(9999) >= MIN_END_YEAR]      # modern era
    w = weights or _weights()
    gram = _gram_of()
    class_mask = {"negative": df["is_gn"], "positive": df["is_gp"], "afb": df["is_tb"]}

    # group-level (Gram-class / TB) pools -- the level the R&D data actually resolves to
    def pool(mask):
        h = df[mask]
        return {"investment": float((h["usd"] * h["w"]).sum()),
                "projects": float(h["w"].sum()),
                "pipeline": int(h.loc[h["Product Name"].notna(), "Product Name"].nunique())}
    group_tot = {g: pool(m) for g, m in class_mask.items()}

    rows = []
    for patho, pats in KEYWORDS.items():
        g = gram.get(patho)
        gt = group_tot.get(g, {"investment": 0.0, "projects": 0.0, "pipeline": 0})
        rgx = re.compile("|".join(pats))
        specific = df["text"].str.contains(rgx, regex=True, na=False)
        sp = df[specific]
        # species attention = its Gram-class broad-spectrum pool + a bonus for work naming it
        rows.append({
            "pathogen": patho, "gram": g,
            "projects": gt["projects"] + float(sp["w"].sum()),
            "investment": gt["investment"] + float((sp["usd"] * sp["w"]).sum()),
            "pipeline": gt["pipeline"] + int(sp.loc[sp["Product Name"].notna(), "Product Name"].nunique()),
            "named_projects": int(specific.sum()),                    # transparency: species-named only
        })
    rai = pd.DataFrame(rows).set_index("pathogen")

    comp = pd.DataFrame(index=rai.index)
    for c in ["investment", "projects", "pipeline"]:
        comp[c + "_n"] = minmax(np.log1p(rai[c]))
    wsum = sum(w[c] for c in ["investment", "projects", "pipeline"])
    rai["RAI"] = sum(comp[c + "_n"] * w[c] for c in ["investment", "projects", "pipeline"]) / wsum
    out = pd.concat([rai, comp], axis=1).sort_values("RAI", ascending=False)
    if save:
        out.to_parquet(PROCESSED_DIR / "rai.parquet")
    return out


if __name__ == "__main__":
    r = compute_rai()
    print(r[["gram", "projects", "investment", "pipeline", "named_projects", "RAI"]].round(2).to_string())
