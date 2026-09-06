"""Output directories and the provenance block stamped on shareable outputs."""
from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

# Written into intermediate_phi/ at runtime so it survives `git clean -fdx`.
PHI_LABEL = """# intermediate_phi

Every patient-level artifact this pipeline produces: the Phase 0 analytic tables
and the model fits that later phases read. One row per encounter block, or per
encounter-block-window.

**This directory never leaves the site.** It is not part of any export bundle
and must not be committed, copied to a shared drive, or read into an analysis
transcript. Only `output/final_no_phi/` is shareable.
"""


def site_dirs(repo: Path) -> dict[str, Path]:
    """Create and return the three directories this pipeline may write to."""
    d = {
        "out_phi": repo / "output" / "intermediate_phi",
        "out_final": repo / "output" / "final_no_phi",
        "logs": repo / "logs",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)

    label = d["out_phi"] / "README.md"
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
