"""Shared small helpers."""
import pandas as pd


def minmax(s: pd.Series) -> pd.Series:
    """Min-max scale a series to 0-1; constant/empty -> 0.5."""
    s = s.astype(float)
    lo, hi = s.min(), s.max()
    if not pd.notna(lo) or hi == lo:
        return pd.Series(0.5, index=s.index)
    return (s - lo) / (hi - lo)
