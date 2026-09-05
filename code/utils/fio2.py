"""FiO2 scale detection and normalisation to a fraction in [0.21, 1.0].

Scale is decided per column; bounds are applied per value. A column is rescaled
only if its distribution says it is on the percent scale; an out-of-range value
in an otherwise-fractional column is nulled, not rescaled. Ambiguous columns
raise. Rationale: design_notes.md §11, config/covariates.json oxygenation.fio2_scale.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

_CFG = json.loads((Path(__file__).resolve().parents[2] / "config" / "covariates.json").read_text())
_F = _CFG["time_varying"]["oxygenation"]["fio2_scale"]

FRACTION_BAND = tuple(_F["fraction_band"])
PERCENT_BAND = tuple(_F["percent_band"])
MIN_SHARE = _F["column_scale_min_share"]


class Fio2ScaleError(ValueError):
    """A fio2 column is neither clearly fraction nor clearly percent."""


@dataclass
class Fio2Report:
    scale: str
    n_total: int
    n_nonnull: int
    n_in_fraction_band: int
    n_in_percent_band: int
    n_neither: int
    rescaled: bool
    n_nulled: int
    examples_nulled: list = field(default_factory=list)

    def __str__(self) -> str:
        if not self.n_nonnull:
            return "fio2: no non-null values"
        pct = lambda n: 100.0 * n / self.n_nonnull
        lines = [
            f"fio2: {self.n_nonnull:,} non-null of {self.n_total:,}",
            f"  in fraction band {FRACTION_BAND}: {self.n_in_fraction_band:,} ({pct(self.n_in_fraction_band):.1f}%)",
            f"  in percent band  {PERCENT_BAND}: {self.n_in_percent_band:,} ({pct(self.n_in_percent_band):.1f}%)",
            f"  in neither:       {self.n_neither:,} ({pct(self.n_neither):.1f}%)",
            f"  scale detected: {self.scale}" + ("  -> RESCALED /100" if self.rescaled else ""),
            f"  nulled as out of range: {self.n_nulled:,}",
        ]
        if self.examples_nulled:
            shown = ", ".join(f"{v:g}" for v in self.examples_nulled)
            lines.append(f"  examples nulled: {shown}")
        return "\n".join(lines)


def detect_fio2_scale(values: pd.Series) -> tuple[str, Fio2Report]:
    """Return (scale, report); raises Fio2ScaleError if the column is ambiguous."""
    v = pd.to_numeric(values, errors="coerce")
    nonnull = v.dropna()
    n = len(nonnull)

    in_frac = int(((nonnull >= FRACTION_BAND[0]) & (nonnull <= FRACTION_BAND[1])).sum())
    in_pct = int(((nonnull >= PERCENT_BAND[0]) & (nonnull <= PERCENT_BAND[1])).sum())
    neither = n - in_frac - in_pct

    rep = Fio2Report(
        scale="unknown", n_total=len(v), n_nonnull=n,
        n_in_fraction_band=in_frac, n_in_percent_band=in_pct,
        n_neither=neither, rescaled=False, n_nulled=0,
    )
    if n == 0:
        raise Fio2ScaleError("fio2 column has no non-null values; cannot determine scale")

    if in_frac / n >= MIN_SHARE:
        rep.scale = "fraction"
    elif in_pct / n >= MIN_SHARE:
        rep.scale = "percent"
    else:
        raise Fio2ScaleError(
            f"fio2 column is neither clearly fraction nor clearly percent: "
            f"{100*in_frac/n:.1f}% in {FRACTION_BAND}, {100*in_pct/n:.1f}% in {PERCENT_BAND}, "
            f"{100*neither/n:.1f}% in neither, against a {100*MIN_SHARE:.0f}% threshold. "
            f"Mixed units in one column must be resolved at the site, not by a heuristic."
        )
    return rep.scale, rep


def normalize_fio2(values: pd.Series) -> tuple[pd.Series, Fio2Report]:
    """Return (fio2 as a fraction, report). Rescales whole columns, nulls stray values."""
    scale, rep = detect_fio2_scale(values)
    out = pd.to_numeric(values, errors="coerce")

    if scale == "percent":
        out = out / 100.0
        rep.rescaled = True

    bad = out.notna() & ((out < FRACTION_BAND[0]) | (out > FRACTION_BAND[1]))
    rep.n_nulled = int(bad.sum())
    if rep.n_nulled:
        orig = pd.to_numeric(values, errors="coerce")[bad]
        rep.examples_nulled = sorted(orig.unique().tolist(), reverse=True)[:5]
    out = out.where(~bad)
    return out, rep
