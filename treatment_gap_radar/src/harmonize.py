"""Build the unified long-format isolate table to data_processed/isolates_long.parquet."""
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .loaders import iter_atlas_long, MIC_LOADERS, LONG_COLS
from .paths import PROCESSED_DIR

OUT = PROCESSED_DIR / "isolates_long.parquet"


def _chunks():
    """Yield (label, long-DataFrame) for every source."""
    for chunk in iter_atlas_long():
        yield "ATLAS", chunk
    for fn in MIC_LOADERS:
        try:
            yield fn.__name__, fn()
        except Exception as e:
            print(f"  !! {fn.__name__} FAILED: {type(e).__name__}: {e}", flush=True)


def build(verbose=True):
    writer = None
    total = 0
    by_src = {}
    try:
        for label, chunk in _chunks():
            if chunk is None or chunk.empty:
                continue
            # enforce consistent dtypes so the parquet schema is stable across sources
            chunk = chunk[LONG_COLS].astype({c: "string" for c in
                        ["source", "isolate_id", "pathogen", "gram", "who", "country_iso3",
                         "country", "region", "age_group", "specimen", "antibiotic",
                         "drug_class", "mic_op", "sir"]})
            chunk["year"] = chunk["year"].astype("Int64")
            chunk["mic"] = pd.to_numeric(chunk["mic"], errors="coerce").astype("float64")
            chunk["age"] = pd.to_numeric(chunk["age"], errors="coerce").astype("float64")
            chunk["pediatric"] = chunk["pediatric"].astype("boolean")
            table = pa.Table.from_pandas(chunk, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(OUT, table.schema)
            else:
                table = table.cast(writer.schema)
            writer.write_table(table)
            total += len(chunk)
            by_src[label] = by_src.get(label, 0) + len(chunk)
            if verbose:
                print(f"  +{len(chunk):>9,} rows  [{label}]  (running {total:,})", flush=True)
    finally:
        if writer is not None:
            writer.close()
    if verbose:
        print(f"\nWrote {total:,} long rows -> {OUT}")
        for k, v in by_src.items():
            print(f"   {k:16} {v:>10,}")
    return total


if __name__ == "__main__":
    build()
