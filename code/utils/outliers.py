"""Applies config/outlier_config.json. The only place bounds are applied.

Medication bounds are two-layered: apply_med_raw() runs before clifpy's unit
conversion and keys on the charted unit; apply_med_converted() runs after. A
missing config raises; a category with no bound is reported, never skipped
silently.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO / "config" / "outlier_config.json"


class OutlierConfigError(RuntimeError):
    """Bounds cannot be loaded, or were asked for in the wrong unit."""


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
    """Bound raw charted doses per (drug, charted unit). Must run before conversion."""
    cfg = config or load_config()
    raw = {k: v for k, v in cfg["med_dose_raw"].items() if not k.startswith("_")}

    out = df.copy()
    rep = OutlierReport(scope="med_dose RAW (pre-conversion)")
    keys = out.groupby([drug_col, unit_col], observed=True).groups
    for (drug, unit), idx in keys.items():
        n = len(idx)
        rep.n_checked += n
        label = f"{drug} [{unit}]"
        # Units are compared case- and space-insensitively: the schema says
        # "units/min" and sites chart "Units/min".
        per_unit = {k.lower().strip(): v for k, v in (raw.get(drug) or {}).items()}
        unit = str(unit).lower().strip()
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
    """Bound converted doses per drug. Must run after conversion and the unit guard."""
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
    """Ceiling past which a dose proves the raw bounds did not run. Derived, not fixed.

    The analysis unit is mcg/hr, and three charted units feed it. A bound only
    constrains the arm it is written in, so the ceiling is the LARGEST value any
    arm can still produce after conversion -- weight-based arms at the maximum
    plausible weight. Taking the mcg/hr bound alone understates it, and the
    assertion then fires on data the bounds did accept.
    """
    cfg = config or load_config()
    raw = cfg["med_dose_raw"]["fentanyl"]
    w_max = cfg["analysis_unit_bounds"]["vitals"]["weight_kg"][1]
    per_arm = {
        "mcg/hr": lambda hi: hi,
        "mcg/kg/hr": lambda hi: hi * w_max,
        "mg/hr": lambda hi: hi * 1000.0,
        "mcg/kg/min": lambda hi: hi * 60.0 * w_max,
        "mcg/min": lambda hi: hi * 60.0,
    }
    return max(per_arm[u](hi) for u, (_, hi) in raw.items() if u in per_arm)


def fentanyl_charted_max_mcg_hr(config: dict | None = None) -> float:
    """The mcg/hr bound itself: high but bound-consistent above this, not a fault."""
    cfg = config or load_config()
    return cfg["med_dose_raw"]["fentanyl"]["mcg/hr"][1]


def nee_sanity_ceiling(coefficients: dict, config: dict | None = None) -> float:
    """Arithmetic ceiling for NEE given the coefficient table."""
    cfg = config or load_config()
    conv = cfg["med_dose_converted"]
    return sum(conv[d][1] * c for d, c in coefficients.items() if d in conv)
