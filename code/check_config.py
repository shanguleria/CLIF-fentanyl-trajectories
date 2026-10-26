"""Is config/config.json usable? Answer in a second, not after a table read.

Reports paths and settings only. It opens no CLIF table and prints no directory
listing, so it stays clear of patient data entirely -- table presence is tested
by filename.

    python3 code/check_config.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO / "config" / "config.json"

# Keys clifpy itself requires; the config is handed to it unchanged.
CLIFPY_REQUIRED = ("data_directory", "filetype", "timezone")


# protocol.json names categories by string and every filter matches LITERALLY,
# so a value that is not in the CLIF mCIDE matches zero rows forever. That is not
# hypothetical: "CPOT" was declared here from 2026-09-30 to 2026-10-26 and is not
# a permissible value -- the mCIDE spells it `cpot_total` -- so at a CPOT site the
# pain half of the headline indication metric would have reported a clean, silent
# zero behind a reassuring "not charted here" line. Checked against clifpy's own
# packaged schema rather than a list copied into this file, because a copied list
# is the next thing to drift.
_SCHEMA_FOR_TABLE = {
    "assessments": ("patient_assessments_schema.yaml", "assessment_category"),
    "labs": ("labs_schema.yaml", "lab_category"),
    "vitals": ("vitals_schema.yaml", "vital_category"),
}

# med_category IS mCIDE-controlled -- 75 permissible values for continuous and
# 165 for intermittent -- and the two lists are NOT the same. dexmedetomidine is
# valid continuous and absent from intermittent, so which list a declared drug is
# checked against depends on how this pipeline reads it. Getting that wrong is
# the same silent failure as the CPOT typo: a filter on a non-permissible value
# matches nothing, forever, without complaint.
_MEDICATION_ROLES = {
    "opioid_infusion_categories": ("medication_admin_continuous_schema.yaml",
                                   "continuous infusions"),
    "other_sedative_categories": ("medication_admin_continuous_schema.yaml",
                                  "continuous infusions"),
    "opioid_bolus_categories": ("medication_admin_intermittent_schema.yaml",
                                "intermittent boluses"),
}


def _permissible(schema_file: str, column: str) -> set[str] | None:
    """Permissible values for one mCIDE column, or None if unavailable."""
    try:
        import yaml
        import clifpy
    except ImportError:
        return None
    f = Path(clifpy.__file__).parent / "schemas" / schema_file
    if not f.is_file():
        return None
    for col in yaml.safe_load(f.read_text()).get("columns", []):
        if col.get("name") == column and "permissible_values" in col:
            return set(col["permissible_values"])
    return None


def declared_categories(protocol: dict) -> list[tuple[str, str, str]]:
    """(table, category, where it was declared) for every category in protocol.json.

    Iterates the config's own keys. A hardcoded list here would be the same bug
    one level up: a newly declared instrument would simply not be checked.
    """
    found = []
    for block in ("time_invariant", "time_varying"):
        for name, spec in (protocol.get(block) or {}).items():
            if name.startswith("_") or not isinstance(spec, dict):
                continue
            src = spec.get("source") or {}
            if src.get("table") in _SCHEMA_FOR_TABLE and src.get("category"):
                found.append((src["table"], src["category"], f"{block}.{name}"))
    ind = protocol.get("indication") or {}
    for key in ("pain_instruments", "sedation_instruments"):
        for cat in (ind.get(key) or {}):
            found.append(("assessments", cat, f"indication.{key}"))
    return found


def check_medications(protocol: dict, problems: list[str]) -> int:
    """Every declared med_category, against the schema for the table it is read from."""
    meds = protocol.get("medications") or {}
    checked = 0
    for role, (schema_file, how) in _MEDICATION_ROLES.items():
        cats = meds.get(role) or []
        if not cats:
            continue
        allowed = _permissible(schema_file, "med_category")
        if allowed is None:
            continue
        for cat in cats:
            checked += 1
            if cat not in allowed:
                problems.append(
                    f"medications.{role} declares med_category={cat!r}, which is "
                    f"not a permissible value for {how}. This pipeline reads that "
                    f"role from {schema_file.split('_schema')[0]}, and the two "
                    f"medication tables do not share one vocabulary.")
    return checked


def check_categories(problems: list[str]) -> None:
    protocol_path = REPO / "config" / "protocol.json"
    if not protocol_path.is_file():
        problems.append(f"protocol.json not found at {protocol_path}")
        return
    protocol = json.loads(protocol_path.read_text())

    cache: dict[str, set[str] | None] = {}
    checked = unverifiable = 0
    for table, cat, where in declared_categories(protocol):
        schema_file, column = _SCHEMA_FOR_TABLE[table]
        if schema_file not in cache:
            cache[schema_file] = _permissible(schema_file, column)
        allowed = cache[schema_file]
        if allowed is None:
            unverifiable += 1
            continue
        checked += 1
        if cat not in allowed:
            near = sorted(v for v in allowed if v.lower() == cat.lower()) or \
                   sorted(v for v in allowed if cat.lower() in v.lower())[:4]
            hint = f" Did you mean {near}?" if near else ""
            problems.append(
                f"{where} declares {column}={cat!r}, which is not a permissible "
                f"value in the CLIF {table} mCIDE. The filter matches literally, "
                f"so this matches zero rows at every site.{hint}")
    checked += check_medications(protocol, problems)
    print(f"mCIDE check    {checked} declared categor(ies) validated"
          + (f", {unverifiable} unverifiable (clifpy schema unavailable)"
             if unverifiable else ""))


def main() -> int:
    if not CONFIG_PATH.is_file():
        print(f"ERROR: config not found at {CONFIG_PATH}\n"
              f"       Fix: cp config/config_template.json config/config.json\n"
              f"       then edit it for your site.")
        return 1

    cfg = json.loads(CONFIG_PATH.read_text())
    problems = []

    for key in ("site_name", "clif_version", *CLIFPY_REQUIRED):
        if not cfg.get(key):
            problems.append(f"{key} is missing or empty.")

    data_dir = Path(str(cfg.get("data_directory", "")))
    filetype = cfg.get("filetype", "parquet")
    tables = cfg.get("tables_in_use", [])

    if "/path/to/" in str(data_dir) or str(data_dir) in ("", "."):
        problems.append(f"data_directory is still a placeholder ({data_dir}).")
    elif not data_dir.is_dir():
        problems.append(f"data_directory does not exist: {data_dir}")
    elif tables:
        # is_dir() alone goes green on a parent holding only a versioned
        # subdirectory -- which is the multi-hour failed run this exists to
        # prevent. Names only; nothing is opened.
        missing = [t for t in tables
                   if not (data_dir / f"clif_{t}.{filetype}").is_file()]
        if missing:
            problems.append(
                f"{len(missing)} of {len(tables)} tables in tables_in_use absent "
                f"from data_directory (as clif_<name>.{filetype}): "
                + ", ".join(missing)
                + "\n    If your data ships a versioned subdirectory, "
                  "data_directory must point INTO it.")

    check_categories(problems)

    print(f"config          {CONFIG_PATH}")
    print(f"site_name       {cfg.get('site_name')}")
    print(f"clif_version    {cfg.get('clif_version')}  (CLIF spec)")
    print(f"dataset_version {cfg.get('dataset_version') or '(none)'}")
    print(f"filetype        {filetype}")
    print(f"timezone        {cfg.get('timezone')}")
    print(f"data_directory  {data_dir}")
    print(f"tables_in_use   {len(tables)}")
    print(f"output          {REPO / 'output' / 'final_no_phi'}")

    if problems:
        print("\nNOT READY:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("\nREADY.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
