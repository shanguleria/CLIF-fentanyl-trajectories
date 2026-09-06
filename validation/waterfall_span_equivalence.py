"""Does trimming respiratory_support before the waterfall change the result?

The waterfall is the whole cost of a Phase 0 run: about an hour on 4.0M rows,
against a few minutes for everything else, because it processes each
hospitalization's entire stay when the analysis window is [anchor, anchor+72h].

Trimming first is only safe if it does not move the values inside that window.
This compares, on a sample:
    A  waterfall the full history, then trim to the window   (current)
    B  trim to [anchor - LEAD, anchor + EXTENT + TRAIL], then waterfall
and reports the disagreement rate on device_category and fio2_set.

Aggregates only. No row-level output.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from clifpy import RespiratorySupport  # noqa: E402

CFG = json.loads((REPO / "config" / "config.json").read_text())
COV = json.loads((REPO / "config" / "covariates.json").read_text())
EXTENT_H = COV["windows"]["granular"]["extent_hours"]
LEAD_H = int(sys.argv[2]) if len(sys.argv) > 2 else 24
TRAIL_H = 24
N_SAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 400

kw = dict(data_directory=CFG["data_directory"], filetype=CFG["filetype"],
          timezone=CFG["timezone"], output_directory=str(REPO / "logs"))


def canonical(df):
    import importlib.util
    sp = importlib.util.spec_from_file_location("b", REPO / "code" / "01_build_cohort.py")
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m._canonicalise_devices(df)


def main() -> None:
    screen = RespiratorySupport.from_file(
        **kw, columns=["hospitalization_id", "device_category", "recorded_dttm"]).df
    imv = screen[screen["device_category"] == "IMV"]
    anchors = imv.groupby("hospitalization_id")["recorded_dttm"].min().rename("anchor")
    ids = sorted(anchors.index.astype(str))[:N_SAMPLE]
    anchors = anchors.loc[ids]
    print(f"sample: {len(ids):,} hospitalizations with an IMV anchor")

    rs = RespiratorySupport.from_file(**kw, filters={"hospitalization_id": ids})
    raw = rs.df.merge(anchors, on="hospitalization_id", how="inner")
    print(f"raw rows: {len(raw):,}")

    import time
    t0 = time.time()
    full = canonical(RespiratorySupport.from_file(
        **kw, filters={"hospitalization_id": ids}).waterfall(
            verbose=False, return_dataframe=True))
    t_full = time.time() - t0
    full = full.merge(anchors, on="hospitalization_id", how="inner")

    keep = ((raw["recorded_dttm"] >= raw["anchor"] - pd.Timedelta(hours=LEAD_H))
            & (raw["recorded_dttm"] <= raw["anchor"] + pd.Timedelta(hours=EXTENT_H + TRAIL_H)))
    trimmed_raw = raw[keep].drop(columns=["anchor"])
    print(f"trimmed rows: {len(trimmed_raw):,} "
          f"({100*len(trimmed_raw)/len(raw):.1f}% of raw)")

    rs2 = RespiratorySupport.from_file(**kw, filters={"hospitalization_id": ids})
    rs2.df = trimmed_raw
    t0 = time.time()
    part = canonical(rs2.waterfall(verbose=False, return_dataframe=True))
    t_part = time.time() - t0
    part = part.merge(anchors, on="hospitalization_id", how="inner")

    def in_window(d):
        return d[(d["recorded_dttm"] >= d["anchor"])
                 & (d["recorded_dttm"] < d["anchor"] + pd.Timedelta(hours=EXTENT_H))]

    a, b = in_window(full), in_window(part)
    key = ["hospitalization_id", "recorded_dttm"]
    a = a.drop_duplicates(key).set_index(key)
    b = b.drop_duplicates(key).set_index(key)
    shared = a.index.intersection(b.index)

    print(f"\nwaterfall time: full {t_full:.1f}s   trimmed {t_part:.1f}s "
          f"({t_full/max(t_part, .01):.1f}x faster)")
    print(f"rows in the analysis window: full {len(a):,}  trimmed {len(b):,}  "
          f"shared {len(shared):,}")

    for col in ("device_category", "fio2_set"):
        av, bv = a.loc[shared, col], b.loc[shared, col]
        both_null = av.isna() & bv.isna()
        if col == "fio2_set":
            diff = (~both_null) & ~np.isclose(
                av.astype("float64").fillna(-1), bv.astype("float64").fillna(-1))
        else:
            diff = (~both_null) & (av.astype("string") != bv.astype("string"))
        n = int(diff.sum())
        print(f"  {col}: {n:,} of {len(shared):,} differ ({100*n/max(len(shared),1):.3f}%)")
        print(f"    non-null full {av.notna().sum():,} | trimmed {bv.notna().sum():,}")
        recovered = int((av.isna() & bv.notna()).sum())
        lost = int((av.notna() & bv.isna()).sum())
        both = av.notna() & bv.notna()
        if col == "fio2_set":
            disagree = int((both & ~np.isclose(av.astype("float64").fillna(0),
                                               bv.astype("float64").fillna(0))).sum())
        else:
            disagree = int((both & (av.astype("string") != bv.astype("string"))).sum())
        print(f"    recovered by trimming (full null -> trimmed value): {recovered:,}")
        print(f"    lost by trimming     (full value -> trimmed null):  {lost:,}")
        print(f"    both present but DISAGREE:                          {disagree:,}")


if __name__ == "__main__":
    main()
