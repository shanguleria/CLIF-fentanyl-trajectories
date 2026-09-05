"""Where this run reads and writes, and the provenance stamped on what it shares.

One site, one data source, one output tree -- the standard CLIF layout. The site
is whatever `config/config.json` declares; pointing the pipeline at a different
source (a MIMIC-to-CLIF conversion, say) means editing `site_name`,
`data_directory` and `dataset_version` in that file, not adding a code path.
"""
from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

# Written into intermediate_phi/ when it is created. Created at runtime rather
# than tracked, so it survives `git clean -fdx`.
PHI_LABEL = """# intermediate_phi

Patient-level intermediates: one row per patient, or per patient-window.

**This directory never leaves the site.** It is not part of any export bundle
and must not be committed, copied to a shared drive, or read into an analysis
transcript. Only `output/final_no_phi/` is shareable.
"""


def site_dirs(repo: Path) -> dict[str, Path]:
    """Create and return the four directories this pipeline may write to."""
    d = {
        "data_phi": repo / "data" / "intermediate_phi",
        "out_phi": repo / "output" / "intermediate_phi",
        "out_final": repo / "output" / "final_no_phi",
        "logs": repo / "logs",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)

    for phi in (d["data_phi"], d["out_phi"]):
        label = phi / "README.md"
        if not label.exists():
            label.write_text(PHI_LABEL)

    return d


def provenance(config: dict) -> dict:
    """The block stamped onto every shareable output."""
    try:
        sha = subprocess.check_output(
            ["git", "describe", "--always", "--dirty"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        sha = "unknown"
    return {
        "site_name": config["site_name"],
        "clif_version": config["clif_version"],          # CLIF SPEC version
        "dataset_version": config.get("dataset_version", ""),  # conversion/ETL release
        "code_version": sha,
        "generated": datetime.now(ZoneInfo(config["timezone"])).isoformat(),
    }
