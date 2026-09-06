"""Content-addressed cache for the respiratory-support waterfall.

The waterfall is the dominant cost of a Phase 0 run and its inputs change rarely.
It is verified per-encounter independent (validation/waterfall_span_equivalence.py),
so rows are cached PER HOSPITALIZATION: a narrowing cohort change is a full hit, a
widening one computes only the new patients, and the cohort never enters the key.

The key pins everything that can change the output. A cache cannot go stale
silently, because any change to an input produces a different key and a miss.
"""
from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
CACHE_DIR = REPO / "output" / "intermediate_phi" / "_cache"


def _source_identity(config: dict) -> str:
    d = Path(config["data_directory"])
    hits = sorted(d.glob("*respiratory_support*"))
    if not hits:
        raise SystemExit(f"no respiratory_support file under {d}")
    st = hits[0].stat()
    return f"{hits[0].name}:{st.st_size}:{int(st.st_mtime)}"


def cache_key(config: dict, transform_fns: list) -> str:
    """Everything that can change the waterfalled rows, hashed.

    transform_fns are the functions we apply around the vendor call; hashing their
    SOURCE means editing them invalidates the cache automatically, rather than
    relying on a version constant someone remembers to bump.
    """
    import clifpy

    parts = [
        _source_identity(config),
        f"clifpy:{getattr(clifpy, '__version__', 'unknown')}",
    ]
    for fn in transform_fns:
        try:
            parts.append(hashlib.sha256(
                inspect.getsource(fn).encode()).hexdigest()[:12])
        except (OSError, TypeError):
            parts.append("nosource")
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def _path(key: str) -> Path:
    return CACHE_DIR / f"waterfall_{key}.parquet"


def load(key: str, hosp_ids: list[str]) -> tuple[pd.DataFrame | None, list[str]]:
    """Return (cached rows for the requested ids, ids still to compute)."""
    f = _path(key)
    if not f.exists():
        return None, list(hosp_ids)
    df = pd.read_parquet(f)
    have = set(df["hospitalization_id"].astype(str))
    want = set(map(str, hosp_ids))
    missing = sorted(want - have)
    hit = df[df["hospitalization_id"].astype(str).isin(want)]
    return hit, missing


def store(key: str, new_rows: pd.DataFrame) -> int:
    """Append newly computed hospitalizations. Returns the cache's total row count."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    f = _path(key)
    if f.exists():
        old = pd.read_parquet(f)
        fresh = set(new_rows["hospitalization_id"].astype(str))
        old = old[~old["hospitalization_id"].astype(str).isin(fresh)]
        out = pd.concat([old, new_rows], ignore_index=True)
    else:
        out = new_rows
    out.to_parquet(f, index=False)
    (CACHE_DIR / f"waterfall_{key}.json").write_text(json.dumps({
        "key": key, "rows": len(out),
        "hospitalizations": int(out["hospitalization_id"].nunique()),
    }, indent=2))
    return len(out)


def describe(key: str) -> str:
    f = _path(key)
    if not f.exists():
        return f"waterfall cache {key}: empty"
    meta = CACHE_DIR / f"waterfall_{key}.json"
    if meta.exists():
        m = json.loads(meta.read_text())
        return (f"waterfall cache {key}: {m['rows']:,} rows, "
                f"{m['hospitalizations']:,} hospitalizations")
    return f"waterfall cache {key}: present"
