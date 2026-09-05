"""
Phase 0 -- build the analytic tables from CLIF.

Outputs:
    data/intermediate_phi/trajectory_long.parquet   one row per patient-window
    data/intermediate_phi/time_to_event.parquet     one row per patient

Column spec: docs/design_notes.md section 10 (Phase 0).
Per-covariate within-window aggregation rules: section 11.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from utils.paths import site_dirs, provenance  # noqa: E402

CONFIG = json.loads((REPO / "config" / "config.json").read_text())


def main() -> None:
    dirs = site_dirs(REPO)
    prov = provenance(CONFIG)

    print(f"site           {CONFIG['site_name']}")
    print(f"clif_version   {CONFIG['clif_version']}")
    print(f"data_directory {CONFIG['data_directory']}")
    print(f"output         {dirs['out_final']}")

    # TODO Phase 0, in order:
    #  1. Load config["tables_in_use"] with clifpy. config.json is already in
    #     clifpy's own shape (data_directory / filetype / timezone), so it can be
    #     handed over directly -- no translation layer.
    #  2. Cohort: age >= cohort.min_age, IMV present; anchor = first
    #     respiratory_support row with device_category == "IMV".
    #  3. Window grid: cohort.window_hours out to cohort.granular_extent_hours
    #     (and the extended grid for the Phase 1 descriptive view).
    #  4. Dose, per window, in mcg/kg/hr:
    #       infusion -> LOCF between charted changes, then TIME-WEIGHTED mean
    #       bolus    -> SUM within window, NEVER carried forward
    #       total    -> infusion + bolus
    #     Plus other_sedative_categories into their own columns.
    #  5. imv_status per window; dose in extubated windows = 0.
    #  6. Covariates per the section 11 aggregation rules.
    #  7. Table 2 (time-to-event), origin = landmark T:
    #       a) 0=censored 1=successful extubation 2=death 3=tracheostomy
    #       b) 0=censored 1=death 2=discharge alive
    #     Mortality from outcomes.mortality_source; unresolved categories
    #     censored, not counted alive.
    raise NotImplementedError("Phase 0 not implemented -- see the TODO block")


if __name__ == "__main__":
    main()
