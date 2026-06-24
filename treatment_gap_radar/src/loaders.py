"""Per-dataset loaders -> common long-format isolate table.

Common schema (one row per isolate x antibiotic result):
    source, isolate_id, pathogen, gram, who, country_iso3, country, region,
    year, age, age_group, pediatric, specimen, antibiotic, drug_class, mic, mic_op, sir

`sir` is 'S'/'I'/'R' where known (ATLAS native; others via breakpoints later), else <NA>.
"""
import numpy as np
import pandas as pd

from . import canon
from .mic import parse_mic_series
from .paths import raw_path

LONG_COLS = ["source", "isolate_id", "pathogen", "gram", "who", "country_iso3",
             "country", "region", "year", "age", "age_group", "pediatric",
             "specimen", "antibiotic", "drug_class", "mic", "mic_op", "sir"]

SIR_MAP = {"susceptible": "S", "intermediate": "I", "resistant": "R",
           "s": "S", "i": "I", "r": "R", "nonsusceptible": "R", "non-susceptible": "R"}


def _map_unique(series, fn):
    """Apply fn over unique values only, then broadcast back."""
    u = series.dropna().unique()
    table = {v: fn(v) for v in u}
    return series.map(table)


def _add_pathogen_cols(df, species_series):
    pj = _map_unique(species_series.astype(str), canon.canon_pathogen)
    df["pathogen"] = pj.map(lambda t: t[0] if isinstance(t, tuple) else None)
    df["gram"] = pj.map(lambda t: t[1] if isinstance(t, tuple) else None)
    df["who"] = pj.map(lambda t: t[2] if isinstance(t, tuple) else None)


def _pediatric(age_group):
    """ATLAS-style band string or numeric age -> pediatric bool (<18)."""
    if age_group is None:
        return pd.NA
    s = str(age_group).strip().lower()
    if s in ("0 - 17", "0-17"):
        return True
    if s in ("18 - 30", "31 - 60", "61+", "18-30", "31-60"):
        return False
    return pd.NA


def _finalize(df):
    """Ensure all LONG_COLS exist, correct order/dtypes; drop empty result rows."""
    for c in LONG_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    df = df[LONG_COLS]
    # keep rows that have at least a parsed MIC or an S/I/R
    df = df[df["mic"].notna() | df["sir"].notna()].copy()
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["mic"] = pd.to_numeric(df["mic"], errors="coerce").astype(float)
    return df


# --------------------------------------------------------------------------- ATLAS
def iter_atlas_long(chunksize=200_000):
    """Stream ATLAS as long-format chunks (native S/I/R interpretation)."""
    f = raw_path("ATLAS")
    head = pd.read_csv(f, nrows=5, low_memory=False)
    cols = list(head.columns)
    drugs = [c for c in cols if (c + "_I") in cols]            # MIC col with paired _I
    meta = ["Isolate Id", "Species", "Country", "Age Group", "Source", "Year"]

    for chunk in pd.read_csv(f, usecols=meta + drugs + [d + "_I" for d in drugs],
                             dtype=str, chunksize=chunksize, low_memory=False):
        n = len(chunk)
        base = chunk[meta].reset_index(drop=True)
        # var-major melt keeps MIC and _I aligned row-for-row
        mic_long = chunk[drugs].reset_index(drop=True).melt(
            var_name="antibiotic_raw", value_name="mic_raw")
        sir_long = chunk[[d + "_I" for d in drugs]].reset_index(drop=True).melt(
            value_name="sir_raw")["sir_raw"]
        meta_rep = pd.concat([base] * len(drugs), ignore_index=True)

        out = meta_rep
        out["antibiotic_raw"] = mic_long["antibiotic_raw"].values
        out["mic_raw"] = mic_long["mic_raw"].values
        out["sir_raw"] = sir_long.values
        out = out[out["mic_raw"].notna() | out["sir_raw"].notna()].copy()
        if out.empty:
            continue

        dj = _map_unique(out["antibiotic_raw"], lambda d: canon.canon_drug(d))
        out["antibiotic"] = dj.map(lambda t: t[0])
        out["drug_class"] = dj.map(lambda t: t[1])
        out = out[out["antibiotic"].notna()]

        out["mic"], out["mic_op"] = parse_mic_series(out["mic_raw"])
        out["sir"] = out["sir_raw"].astype(str).str.strip().str.lower().map(SIR_MAP).astype("string")

        _add_pathogen_cols(out, out["Species"])
        out["country"] = out["Country"]
        out["country_iso3"] = _map_unique(out["Country"].astype(str), canon.canon_country)
        out["region"] = pd.NA
        out["year"] = out["Year"]
        out["age_group"] = out["Age Group"]
        out["age"] = pd.NA
        out["pediatric"] = _map_unique(out["Age Group"].astype(str), _pediatric)
        out["specimen"] = out["Source"]
        out["isolate_id"] = "ATLAS:" + out["Isolate Id"].astype(str)
        out["source"] = "ATLAS"
        yield _finalize(out)
