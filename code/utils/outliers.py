"""
Outlier bounds: the single place bounds are applied.

Every bound comes from config/outlier_config.json. No bound is ever hardcoded in
a pipeline script -- change one there and it applies everywhere.

THREE LOUD-FAILURE BEHAVIOURS, all deliberate, all tested. Each replaces a
silent failure that a sibling CLIF repo actually shipped:

  1. A missing config file RAISES. Skipping outlier handling silently changes
     results, so "the file wasn't there" must never be a quiet no-op.
  2. A category present in the data with no bound in the config is PRINTED.
     A category with no entry looks exactly like a category that was checked and
     found clean, so the gap is reported rather than assumed away.
  3. Raw and converted medication doses go through DIFFERENT functions. A bound
     is only comparable to a value already in its unit; applying a converted
     bound to a raw charted dose nulls correct data in the wrong unit.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO / "config" / "outlier_config.json"


class OutlierConfigError(RuntimeError):
    """Raised when the bounds cannot be loaded or are asked for in the wrong unit."""


def load_config(path: Path | None = None) -> dict:
    p = Path(path) if path else CONFIG_PATH
    if not p.exists():
        raise OutlierConfigError(
            f"Refusing to continue: outlier config not found at {p}. "
            f"Skipping outlier handling silently changes results, so this is fatal "
            f"rather than a warning."
        )
    return json.loads(p.read_text())


@dataclass
class OutlierReport:
    scope: str
    n_checked: int = 0
    n_nulled: int = 0
    per_key: dict = field(default_factory=dict)     # key -> (n_checked, n_nulled)
    unbounded: list = field(default_factory=list)   # keys present in data, absent from config

    def __str__(self) -> str:
        lines = [f"outliers [{self.scope}]: {self.n_nulled:,} of {self.n_checked:,} nulled"]
        for k, (chk, nul) in sorted(self.per_key.items()):
            if nul:
                lines.append(f"    {k}: {nul:,} of {chk:,}")
        for k in sorted(self.unbounded):
            lines.append(f"    NO BOUND for {k}  <- present in data, absent from config")
        return "\n".join(lines)


def _clip(values: pd.Series, lo: float, hi: float) -> tuple[pd.Series, int]:
    v = pd.to_numeric(values, errors="coerce")
    bad = v.notna() & ((v < lo) | (v > hi))
    return v.where(~bad), int(bad.sum())


def apply_long(
    df: pd.DataFrame, table: str, category_col: str, value_col: str,
    config: dict | None = None,
) -> tuple[pd.DataFrame, OutlierReport]:
    """Bound a long CLIF table (labs, vitals, patient_assessments) per category."""
    cfg = config or load_config()
    bounds = cfg["analysis_unit_bounds"].get(table)
    if bounds is None:
        raise OutlierConfigError(
            f"no analysis_unit_bounds entry for table {table!r}. "
            f"Known: {sorted(k for k in cfg['analysis_unit_bounds'] if not k.startswith('_'))}. "
            f"An unbounded frame must never look bounded, so this raises."
        )
    bounds = {k: v for k, v in bounds.items() if not k.startswith("_")}

    out = df.copy()
    rep = OutlierReport(scope=f"{table} (long)")
    for cat, idx in out.groupby(category_col, observed=True).groups.items():
        n = len(idx)
        rep.n_checked += n
        if cat not in bounds:
            rep.unbounded.append(str(cat))
            rep.per_key[str(cat)] = (n, 0)
            continue
        lo, hi = bounds[cat]
        cleaned, n_bad = _clip(out.loc[idx, value_col], lo, hi)
        out.loc[idx, value_col] = cleaned
        rep.n_nulled += n_bad
        rep.per_key[str(cat)] = (n, n_bad)
    return out, rep


def apply_wide(
    df: pd.DataFrame, table: str, config: dict | None = None,
) -> tuple[pd.DataFrame, OutlierReport]:
    """Bound a wide frame whose column names are the category names."""
    cfg = config or load_config()
    bounds = cfg["analysis_unit_bounds"].get(table)
    if bounds is None:
        raise OutlierConfigError(f"no analysis_unit_bounds entry for table {table!r}")
    bounds = {k: v for k, v in bounds.items() if not k.startswith("_")}

    out = df.copy()
    rep = OutlierReport(scope=f"{table} (wide)")
    for col in out.columns:
        if col not in bounds:
            continue
        lo, hi = bounds[col]
        n = int(out[col].notna().sum())
        cleaned, n_bad = _clip(out[col], lo, hi)
        out[col] = cleaned
        rep.n_checked += n
        rep.n_nulled += n_bad
        rep.per_key[col] = (n, n_bad)
    return out, rep


def apply_med_raw(
    df: pd.DataFrame, drug_col: str, dose_col: str, unit_col: str,
    config: dict | None = None,
) -> tuple[pd.DataFrame, OutlierReport]:
    """Bound RAW charted medication doses, per (drug, charted unit).

    Must run BEFORE unit conversion. A (drug, unit) pair with no bound is
    reported, because that is the half of the problem the unit guard cannot see:
    those values are already in their charted unit, so nothing about them is a
    conversion failure -- they are simply wrong at source.
    """
    cfg = config or load_config()
    raw = {k: v for k, v in cfg["med_dose_raw"].items() if not k.startswith("_")}

    out = df.copy()
    rep = OutlierReport(scope="med_dose RAW (pre-conversion)")
    keys = out.groupby([drug_col, unit_col], observed=True).groups
    for (drug, unit), idx in keys.items():
        n = len(idx)
        rep.n_checked += n
        label = f"{drug} [{unit}]"
        per_unit = raw.get(drug)
        if not per_unit or unit not in per_unit:
            rep.unbounded.append(label)
            rep.per_key[label] = (n, 0)
            continue
        lo, hi = per_unit[unit]
        cleaned, n_bad = _clip(out.loc[idx, dose_col], lo, hi)
        out.loc[idx, dose_col] = cleaned
        rep.n_nulled += n_bad
        rep.per_key[label] = (n, n_bad)
    return out, rep


def apply_med_converted(
    df: pd.DataFrame, drug_col: str, dose_col: str, config: dict | None = None,
) -> tuple[pd.DataFrame, OutlierReport]:
    """Bound CONVERTED medication doses, per drug, in the analysis unit.

    Must run AFTER unit conversion and after the unit guard has confirmed the
    returned unit string.
    """
    cfg = config or load_config()
    conv = {k: v for k, v in cfg["med_dose_converted"].items() if not k.startswith("_")}

    out = df.copy()
    rep = OutlierReport(scope="med_dose CONVERTED (post-conversion)")
    for drug, idx in out.groupby(drug_col, observed=True).groups.items():
        n = len(idx)
        rep.n_checked += n
        if drug not in conv:
            rep.unbounded.append(str(drug))
            rep.per_key[str(drug)] = (n, 0)
            continue
        lo, hi = conv[drug]
        cleaned, n_bad = _clip(out.loc[idx, dose_col], lo, hi)
        out.loc[idx, dose_col] = cleaned
        rep.n_nulled += n_bad
        rep.per_key[str(drug)] = (n, n_bad)
    return out, rep


def fentanyl_sanity_ceiling(config: dict | None = None) -> float:
    """The value past which a total_dose is PROOF the bounds did not run.

    Derived, not written down, so it cannot drift when a bound changes.
    Not a clinical limit -- a real value is 0.5 to 5 mcg/kg/hr.
    """
    cfg = config or load_config()
    max_mcg_hr = cfg["med_dose_raw"]["fentanyl"]["mcg/hr"][1]
    min_weight = cfg["analysis_unit_bounds"]["vitals"]["weight_kg"][0]
    return max_mcg_hr / min_weight


def nee_sanity_ceiling(coefficients: dict, config: dict | None = None) -> float:
    """Arithmetic ceiling for NEE, given the coefficient table.

    CRRT-dose-lmtp measured nee = 8,001 against a ceiling of 17.45 before its
    unit guard existed. Asserting against this surfaces that in one line.
    """
    cfg = config or load_config()
    conv = cfg["med_dose_converted"]
    return sum(conv[d][1] * c for d, c in coefficients.items() if d in conv)
