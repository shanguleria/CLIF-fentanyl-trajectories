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
