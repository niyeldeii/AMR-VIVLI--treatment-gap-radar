"""Per-dataset loaders -> common long-format isolate table.

Common schema (one row per isolate x antibiotic result):
    source, isolate_id, pathogen, gram, who, country_iso3, country, region,
    year, age, age_group, pediatric, specimen, antibiotic, drug_class, mic, mic_op, sir

`sir` is 'S'/'I'/'R' where known (ATLAS native; others via breakpoints later), else <NA>.
"""
import numpy as np
import pandas as pd

import re as _re

from . import canon
from . import breakpoints as bp
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


# ----------------------------------------------------- generic MIC-only datasets
def _clean_drug_header(col):
    """Strip MIC/result noise from a column header before canonical lookup."""
    s = str(col).replace("\n", " ")
    s = _re.sub(r"[/_]", " ", s)          # separators -> space (so \b sees 'mic' in CAZ_MIC)
    s = _re.sub(r"\(.*?\)", " ", s)
    s = _re.sub(r"\b(mic|mic50|mic90|broth|mgit|fixed at \d+|at \d+\s*\w*)\b", " ", s, flags=_re.I)
    return _re.sub(r"\s+", " ", s).strip()


def detect_drug_cols(columns, meta_cols, cutoff=94):
    """Return {raw_col: (canonical, class)} for columns that resolve to a known drug."""
    out = {}
    for c in columns:
        if c in meta_cols:
            continue
        canon_name, cls = canon.canon_drug(_clean_drug_header(c), fuzzy=True, cutoff=cutoff)
        if canon_name:
            out[c] = (canon_name, cls)
    return out


def melt_mic_dataset(df, source, *, id_col=None, species_col, country_col=None,
                     region_col=None, year=None, age_col=None, age_group_col=None,
                     specimen_col=None, drug_cols=None, meta_cols=None):
    """Generic wide-MIC -> long loader. Derives S/I/R via CLSI breakpoints.

    `year` may be a column name or a precomputed Series aligned to df.index.
    """
    df = df.copy()
    meta_cols = set(meta_cols or [])
    meta_cols |= {c for c in [id_col, species_col, country_col, region_col, age_col,
                              age_group_col, specimen_col] if c}
    if isinstance(year, str):
        meta_cols.add(year)
    drugmap = drug_cols or detect_drug_cols(df.columns, meta_cols)
    if not drugmap:
        raise ValueError(f"{source}: no drug columns detected")

    base = pd.DataFrame(index=df.index)
    base["species_raw"] = df[species_col].astype(str)
    base["country"] = df[country_col] if country_col else pd.NA
    base["region"] = df[region_col] if region_col else pd.NA
    base["age"] = pd.to_numeric(df[age_col], errors="coerce") if age_col else pd.NA
    base["age_group"] = df[age_group_col] if age_group_col else pd.NA
    base["specimen"] = df[specimen_col] if specimen_col else pd.NA
    if isinstance(year, str):
        base["year"] = pd.to_numeric(df[year], errors="coerce")
    elif year is not None:
        base["year"] = pd.to_numeric(pd.Series(year, index=df.index), errors="coerce")
    else:
        base["year"] = pd.NA
    base["isolate_id"] = (source + ":" +
                          (df[id_col].astype(str) if id_col else pd.Series(df.index.astype(str), index=df.index)))

    frames = []
    for raw, (cdrug, cls) in drugmap.items():
        sub = base.copy()
        sub["antibiotic"] = cdrug
        sub["drug_class"] = cls
        sub["mic"], sub["mic_op"] = parse_mic_series(df[raw])
        frames.append(sub[sub["mic"].notna()])
    out = pd.concat(frames, ignore_index=True)

    _add_pathogen_cols(out, out["species_raw"])
    out["country_iso3"] = _map_unique(out["country"].astype(str), canon.canon_country)
    # pediatric from numeric age if present
    if age_col:
        out["pediatric"] = out["age"].map(lambda a: (a < 18) if pd.notna(a) else pd.NA)
    elif age_group_col:
        out["pediatric"] = _map_unique(out["age_group"].astype(str), _pediatric)
    else:
        out["pediatric"] = pd.NA
    out["source"] = source

    # S/I/R via breakpoints (per unique pathogen/gram/drug/mic to limit calls)
    key = out[["pathogen", "gram", "antibiotic", "mic", "mic_op"]].drop_duplicates()
    key["sir"] = key.apply(lambda r: bp.interpret(r["pathogen"], r["gram"], r["antibiotic"],
                                                  r["mic"], r["mic_op"]), axis=1)
    out = out.merge(key, on=["pathogen", "gram", "antibiotic", "mic", "mic_op"], how="left")
    out["sir"] = out["sir"].astype("string")
    return _finalize(out)


# --------------- per-dataset thin wrappers (handle each file's quirks) ----------
def load_soar_201818():
    df = pd.read_csv(raw_path("SOAR_201818"), dtype=str, low_memory=False)
    return melt_mic_dataset(df, "SOAR_201818", id_col="IHMANUMBER", species_col="ORGANISMNAME",
                            country_col="COUNTRY", region_col="REGION", year="YEARCOLLECTED",
                            age_col="AGE", specimen_col="BODYLOCATION",
                            meta_cols={"GENDER", "BETALACTAMASE", "DEID_CAT_AGE", "INVESTIGATORNAME"})


def load_soar_201910():
    df = pd.read_excel(raw_path("SOAR_201910"))
    yr = pd.to_datetime(df["Collection Date"], errors="coerce").dt.year
    return melt_mic_dataset(df, "SOAR_201910", id_col="Isolate Number", species_col="Organism",
                            country_col="Country", year=yr, age_col="Age",
                            specimen_col="BodyLocation",
                            meta_cols={"Centre", "Gender", "Betalactamase", "Collection Date"})


def load_soar_207965():
    df = pd.read_excel(raw_path("SOAR_207965"), sheet_name="Sheet2")
    return melt_mic_dataset(df, "SOAR_207965", id_col="IHMA #", species_col="FinalOrganismName",
                            country_col="Country", region_col="Region", year="YearCollected",
                            age_col="Age", specimen_col="BodyLocation",
                            meta_cols={"Investigator", "InvestigatorName", "OriginalOrganismName",
                                       "OrganismFamilyName", "GramType", "Gender", "FacilityName",
                                       "Evaluable", "Beta Lactamase"})


def load_innoviva():
    df = pd.read_excel(raw_path("INNOVIVA_ACINETO"))
    return melt_mic_dataset(df, "INNOVIVA", id_col="Vivli No.", species_col="OrganismName",
                            country_col="Country", region_col="Region", year="YearCollected",
                            age_col="Age", specimen_col="BodyLocation",
                            meta_cols={"Gender", "FacilityName"})


def load_sidero():
    df = pd.read_excel(raw_path("SIDERO_WT"), sheet_name="Five year Surveillance data")
    return melt_mic_dataset(df, "SIDERO_WT", species_col="Organism Name", country_col="Country",
                            region_col="Region", year="Year Collected", specimen_col="Body Location",
                            meta_cols={"Date Collected"})


def load_gears():
    df = pd.read_excel(raw_path("GEARS"), sheet_name="Data")
    return melt_mic_dataset(df, "GEARS", id_col="Isolate", species_col="Organism",
                            country_col="Country", region_col="Region", year="Year",
                            age_col="Age", specimen_col="BodySite",
                            meta_cols={"Family", "Gender", "Facility"})


def load_keystone():
    df = pd.read_excel(raw_path("KEYSTONE"), sheet_name="Line List")
    df.columns = [str(c).replace("\n", " ").strip() for c in df.columns]
    return melt_mic_dataset(df, "KEYSTONE", species_col="Organism", country_col="Country",
                            year="Study Year", age_col="Age", specimen_col="Specimen Type",
                            meta_cols={"Collection Number", "Continent", "US Census Division",
                                       "Nosocomial", "Gender", "Medical Service", "Infection Source",
                                       "Infection Type", "Source of Bloodstream infection",
                                       "Ventilator-Associated Pneumonia", "Intensive Care Unit (ICU)",
                                       "Cystic Fibrosis (CF) Patient"})


def load_gasar():
    df = pd.read_excel(raw_path("GASAR_III"), sheet_name="Sheet1")  # header row 0; leading NaN rows drop out
    return melt_mic_dataset(df, "GASAR_III", id_col="Isolate ID", species_col="Species",
                            country_col="Country", year="Year", specimen_col="Source",
                            meta_cols={"Gene Combination", "Phenotypic Combination"})


def load_plea_i():
    df = pd.read_excel(raw_path("PLEA_I"), sheet_name="Sheet 1", header=1)
    return melt_mic_dataset(df, "PLEA_I", id_col="Isolate ID", species_col="Species",
                            country_col="Country", year="Year", specimen_col="Source",
                            meta_cols={"Gene Combination", "Phenotypic Combination"})


def load_plea_ii():
    df = pd.read_excel(raw_path("PLEA_II"), sheet_name="MIC of Ploymyxin ", header=2)
    return melt_mic_dataset(df, "PLEA_II", id_col="Isolate ID", species_col="Species",
                            country_col="Country", year="Year", specimen_col="Source",
                            meta_cols={"Phenotypic Combination"})


def load_dream():
    """M. tuberculosis (DREAM). MIC columns -> TB critical-concentration breakpoints."""
    df = pd.read_excel(raw_path("DREAM_TB"), sheet_name="DREAM Dataset")
    gene_cols = [c for c in df.columns if any(k in str(c) for k in
                 ("Rv0678", "atpE", "pepQ", "Rv1979", "_NT", "_AA"))]
    return melt_mic_dataset(df, "DREAM_TB", species_col="Organism", country_col="Country",
                            region_col="Continent", year="Year Collected", specimen_col="Specimen",
                            meta_cols=set(gene_cols) | {"SubType"})


def load_spidaar():
    """SPIDAAR patient-level RWE (Kenya/Ghana/Uganda/Malawi). Phenotype flags -> long rows.

    c3r (3GC resistance) -> Ceftriaxone S/R for Gram-negatives (excl. intrinsic-AmpC
    Enterobacter); mrsa -> Oxacillin S/R for S. aureus. mdr flag saved separately.
    """
    from .paths import PROCESSED_DIR
    df = pd.read_excel(raw_path("SPIDAAR_ISOLATE"), sheet_name="data")
    df["year"] = df["sampdat"].astype(str).str.extract(r"(\d{4})")[0]
    pj = _map_unique(df["isolate"].astype(str), canon.canon_pathogen)
    df["pathogen"] = pj.map(lambda t: t[0]); df["gram"] = pj.map(lambda t: t[1]); df["who"] = pj.map(lambda t: t[2])
    df = df[df["pathogen"].notna()].copy()
    df["isolate_id"] = "SPIDAAR:" + df["iid"].astype(str)

    rows = []
    # c3r -> Ceftriaxone (Gram-negative, exclude Enterobacter intrinsic AmpC)
    gn = df[(df["gram"] == "negative") & (~df["pathogen"].str.startswith("Enterobacter"))]
    c = gn[gn["c3r"].isin([0, 1])].copy()
    c["antibiotic"] = "Ceftriaxone"; c["drug_class"] = "Cephalosporin"
    c["sir"] = c["c3r"].map({1: "R", 0: "S"})
    rows.append(c)
    # mrsa -> Oxacillin (S. aureus)
    sa = df[(df["pathogen"] == "Staphylococcus aureus") & (df["mrsa"].isin([0, 1]))].copy()
    sa["antibiotic"] = "Oxacillin"; sa["drug_class"] = "Penicillin"
    sa["sir"] = sa["mrsa"].map({1: "R", 0: "S"})
    rows.append(sa)

    out = pd.concat(rows, ignore_index=True)
    out["country"] = out["ctry"]
    out["country_iso3"] = _map_unique(out["ctry"].astype(str), canon.canon_country)
    out["specimen"] = out["stype"]
    out["region"] = "Africa"
    out["age"] = pd.NA; out["age_group"] = pd.NA; out["pediatric"] = pd.NA
    out["mic"] = pd.NA; out["mic_op"] = pd.NA
    out["source"] = "SPIDAAR"
    out["sir"] = out["sir"].astype("string")

    # supplementary: isolate-level MDR flag (0/1; 99=unknown dropped)
    mdr = df[df["mdr"].isin([0, 1])][["isolate_id", "pathogen", "gram", "who", "ctry", "mdr"]].copy()
    mdr.rename(columns={"ctry": "country"}).to_parquet(PROCESSED_DIR / "spidaar_mdr.parquet")
    return _finalize(out)


MIC_LOADERS = [load_soar_201818, load_soar_201910, load_soar_207965,
               load_innoviva, load_sidero, load_gears, load_keystone,
               load_gasar, load_plea_i, load_plea_ii, load_dream, load_spidaar]
