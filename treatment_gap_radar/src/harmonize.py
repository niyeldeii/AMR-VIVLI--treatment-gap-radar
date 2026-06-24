"""Build the unified long-format isolate table to data_processed/isolates_long.parquet."""
import pyarrow as pa
import pyarrow.parquet as pq

from .loaders import iter_atlas_long, LONG_COLS
from .paths import PROCESSED_DIR

OUT = PROCESSED_DIR / "isolates_long.parquet"


def build(verbose=True):
    writer = None
    total = 0
    try:
        sources = list(iter_atlas_long())  # generator of chunks
        for chunk in sources:
            if chunk.empty:
                continue
            table = pa.Table.from_pandas(chunk[LONG_COLS], preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(OUT, table.schema)
            writer.write_table(table)
            total += len(chunk)
            if verbose:
                print(f"  +{len(chunk):>9,} rows (running {total:,})", flush=True)
    finally:
        if writer is not None:
            writer.close()
    if verbose:
        print(f"Wrote {total:,} long rows -> {OUT}")
    return total


if __name__ == "__main__":
    build()
