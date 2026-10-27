"""Medication dose unit conversion, driven by config/protocol.json dose_units.

Done here rather than by clifpy, which leaves the raw value in place when it
cannot convert and reports the failure only in the returned unit string. A
(drug, charted unit) pair with no entry raises; a row that needs a weight and
has none is nulled and counted.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
_COV = json.loads((REPO / "config" / "protocol.json").read_text())
TABLES = {k: v for k, v in _COV["dose_units"].items() if not k.startswith("_")}

TARGET_OF: dict[str, str] = {}
for _target, _spec in TABLES.items():
    for _drug in _spec.get("_target_for", []):
        TARGET_OF[_drug] = _target


class DoseUnitError(ValueError):
    """A (drug, charted unit) pair the conversion table does not cover."""


@dataclass
class DoseReport:
    n_rows: int = 0
    n_converted: int = 0
    n_no_weight: int = 0
    accepted: dict = field(default_factory=dict)
    per_pair: dict = field(default_factory=dict)

    def __str__(self) -> str:
        lines = [f"dose conversion: {self.n_converted:,} of {self.n_rows:,} converted"]
        for k, v in sorted(self.per_pair.items()):
            lines.append(f"    {k}: {v:,}")
        for k, v in sorted(self.accepted.items()):
            lines.append(f"    accepted as unconverted (left NaN): {k}: {v:,}")
        if self.n_no_weight:
            lines.append(f"    nulled for want of a weight: {self.n_no_weight:,}")
        return "\n".join(lines)


def target_unit(drug: str) -> str:
    if drug not in TARGET_OF:
        raise DoseUnitError(
            f"no target unit declared for {drug!r}; known: {sorted(TARGET_OF)}"
        )
    return TARGET_OF[drug]


def convert(df: pd.DataFrame, drug_col: str, dose_col: str, unit_col: str,
            weight_col: str = "weight_kg",
            accepted: dict[str, list[str]] | None = None) -> tuple[pd.Series, DoseReport]:
    """Convert charted doses to each drug's target unit.

    Returns (converted series, report). Raises on a (drug, unit) pair the table
    does not cover -- dropping one would remove a drug from the score for exactly
    the patients who received it. The one exception is `accepted` ({drug: [units]},
    from the site's config.json): a pair the site has chosen to leave unconverted
    stays NaN and is counted in the report. Anything not listed still raises.
    """
    accepted = {d: {u.lower().strip() for u in us} for d, us in (accepted or {}).items()}
    drug = df[drug_col].astype("string")
    unit = df[unit_col].astype("string").str.lower().str.strip()
    dose = pd.to_numeric(df[dose_col], errors="coerce")
    weight = pd.to_numeric(df[weight_col], errors="coerce") if weight_col in df else None

    out = pd.Series(np.nan, index=df.index, dtype="float64")
    rep = DoseReport(n_rows=int(dose.notna().sum()))
    unknown: list[str] = []

    for (d, u), idx in df.groupby([drug, unit], observed=True, dropna=False).groups.items():
        if pd.isna(d) or pd.isna(u):
            continue
        table = TABLES.get(TARGET_OF.get(d, ""), {})
        spec = table.get(u)
        if spec is None:
            n_dose = int(dose.loc[idx].notna().sum())
            if n_dose and u in accepted.get(d, ()):
                rep.accepted[f"{d} [{u}]"] = n_dose
            elif n_dose:
                unknown.append(f"{d} [{u}]")
            continue
        vals = dose.loc[idx] * spec["factor"]
        if spec["weight"] != "none":
            if weight is None:
                unknown.append(f"{d} [{u}] (needs a weight column)")
                continue
            w = weight.loc[idx]
            vals = vals / w if spec["weight"] == "divide" else vals * w
            rep.n_no_weight += int((dose.loc[idx].notna() & w.isna()).sum())
        out.loc[idx] = vals
        n = int(vals.notna().sum())
        rep.n_converted += n
        rep.per_pair[f"{d} [{u}] -> {TARGET_OF[d]}"] = n

    if unknown:
        raise DoseUnitError(
            "no conversion for: " + ", ".join(sorted(set(unknown)))
            + ". Add them to config/protocol.json dose_units rather than dropping "
              "them -- a dropped unit removes a drug for the patients charted that way. If a "
              "unit cannot be converted (e.g. a bare amount) and you accept losing those "
              "rows, list it in config.json accepted_unhandled_dose_units and confirm the "
              "share of rows with your data owner first."
        )
    return out, rep
