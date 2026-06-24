"""Canonical normalization for drugs, pathogens, and countries."""
import re
import functools
import yaml
import pycountry
from rapidfuzz import process, fuzz
from .paths import CONFIG_DIR

_WS = re.compile(r"[\s/_\-]+")
_PAREN = re.compile(r"\(.*?\)")


def norm(s) -> str:
    """Lowercase, drop parentheticals, collapse separators/newlines to single spaces."""
    if s is None:
        return ""
    s = str(s).replace("\n", " ").strip().lower()
    s = _PAREN.sub(" ", s)
    s = _WS.sub(" ", s)
    return s.strip()


@functools.lru_cache(maxsize=1)
def _load(name):
    with open(CONFIG_DIR / name, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


@functools.lru_cache(maxsize=1)
def _drug_index():
    """alias(normalized) -> (canonical, class)."""
    idx = {}
    for canon, meta in _load("antibiotics.yaml").items():
        idx[norm(canon)] = (canon, meta["class"])
        for a in meta.get("aliases", []) or []:
            idx[norm(a)] = (canon, meta["class"])
    return idx


@functools.lru_cache(maxsize=1)
def _pathogen_index():
    """alias(normalized) -> (canonical, gram, who)."""
    idx = {}
    for canon, meta in _load("pathogens.yaml").items():
        rec = (canon, meta["gram"], meta.get("who", "") or "")
        idx[norm(canon)] = rec
        for a in meta.get("aliases", []) or []:
            idx[norm(a)] = rec
    return idx


# phenotype tags some datasets append to species (e.g. "Staphylococcus aureus, MSSA")
_PHENO = re.compile(r",?\s*(mssa|mrsa|vre|esbl|cre|mdr|xdr)\b", re.I)


def canon_drug(name, fuzzy=True, cutoff=90):
    """Return (canonical, class) or (None, None)."""
    n = norm(name)
    idx = _drug_index()
    if n in idx:
        return idx[n]
    if fuzzy and n:
        m = process.extractOne(n, idx.keys(), scorer=fuzz.token_sort_ratio)
        if m and m[1] >= cutoff:
            return idx[m[0]]
    return (None, None)


def canon_pathogen(name, fuzzy=True, cutoff=90):
    """Return (canonical, gram, who) or (None, None, None)."""
    raw = _PHENO.sub("", str(name)) if name is not None else ""
    n = norm(raw)
    idx = _pathogen_index()
    if n in idx:
        return idx[n]
    if fuzzy and n:
        m = process.extractOne(n, idx.keys(), scorer=fuzz.token_sort_ratio)
        if m and m[1] >= cutoff:
            return idx[m[0]]
    return (None, None, None)


@functools.lru_cache(maxsize=4096)
def canon_country(name):
    """Return ISO-3166 alpha-3 code or None."""
    if not name or str(name).strip().lower() in ("", "nan", "unknown"):
        return None
    raw = str(name).strip()
    fixes = {
        "Russia": "RUS", "Russian Federation": "RUS", "South Korea": "KOR",
        "Korea, South": "KOR", "Korea": "KOR", "Czech Republic": "CZE",
        "Slovak Republic": "SVK", "USA": "USA", "United States": "USA",
        "UK": "GBR", "Vietnam": "VNM", "Iran": "IRN", "Venezuela": "VEN",
        "Taiwan": "TWN", "Hong Kong": "HKG", "Ivory Coast": "CIV",
        "Cote D'Ivoire": "CIV",
    }
    if raw in fixes:
        return fixes[raw]
    try:
        c = pycountry.countries.lookup(raw)
        return c.alpha_3
    except LookupError:
        try:
            res = pycountry.countries.search_fuzzy(raw)
            if res:
                return res[0].alpha_3
        except LookupError:
            return None
    return None
