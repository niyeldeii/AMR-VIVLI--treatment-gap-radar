"""Parse censored MIC strings into (numeric value, operator)."""
import re
import numpy as np
import pandas as pd

_MIC = re.compile(r"^\s*(<=|>=|<|>|=)?\s*([0-9]*\.?[0-9]+)\s*$")
# non-result tokens that appear in MIC cells
_NULLISH = {"", "nan", "na", "n/a", "nt", "ntd", "-", ".", "none", "null"}


def parse_mic(value):
    """Return (mic_value: float|nan, op: str). op in {<=,>=,<,>,=}.

    '<=0.015' -> (0.015,'<='); '>8' -> (8.0,'>'); '2' -> (2.0,'='); junk -> (nan,'').
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return (np.nan, "")
    s = str(value).strip()
    if s.lower() in _NULLISH:
        return (np.nan, "")
    m = _MIC.match(s)
    if not m:
        return (np.nan, "")
    op = m.group(1) or "="
    try:
        return (float(m.group(2)), op)
    except ValueError:
        return (np.nan, "")


def parse_mic_series(s: pd.Series):
    """Vectorized parse of a MIC column. Returns (values: Series[float], ops: Series[str])."""
    ss = s.astype("string").str.strip()
    ss = ss.mask(ss.str.lower().isin(_NULLISH))
    ext = ss.str.extract(r"^\s*(<=|>=|<|>|=)?\s*([0-9]*\.?[0-9]+)\s*$")
    val = pd.to_numeric(ext[1], errors="coerce").astype(float)
    op = ext[0].fillna("=")
    op = op.where(val.notna(), "")
    return val, op
