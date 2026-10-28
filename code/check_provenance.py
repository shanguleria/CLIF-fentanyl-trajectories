"""Confirm that several sites ran the same analysis.

Run at the coordinating centre over the site bundles: it needs nothing beyond
the `output/final_no_phi/` upload, every directory of which carries a
provenance.json.

    python code/check_provenance.py multisite/ucmc multisite/rush multisite/nu

Protocol divergence (protocol.json, outlier_config.json, definition_version)
exits non-zero. Code divergence is reported but does NOT fail: a site may have
had to extend _to_mcg_hr with a locally charted dose unit to run at all. The
per-file digests name what changed.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

# Must be identical across sites or the pooled estimate is not of one quantity.
PROTOCOL_KEYS = ("protocol_digests", "definition_versions", "clif_version")


def blocks(paths: list[str]) -> list[dict]:
    """Every provenance/manifest block under the given files or directories."""
    out = []
    for p in (Path(x) for x in paths):
        found = (sorted(p.rglob("provenance.json")) + sorted(p.rglob("manifest.json"))
                 if p.is_dir() else [p])
        for f in found:
            try:
                d = json.loads(f.read_text())
            except Exception as e:
                print(f"  SKIP  {f}: {e}")
                continue
            if "site_name" in d:
                d["_file"] = str(f)
                out.append(d)
    return out


def one_per_site(bs: list[dict]) -> dict[str, dict]:
    """Collapse to one block per site, and report any site disagreeing with itself.

    A site whose own blocks differ uploaded two code versions in one bundle,
    which is what a mid-pipeline edit produces.
    """
    by_site: dict[str, list[dict]] = {}
    for b in bs:
        by_site.setdefault(b["site_name"], []).append(b)
    out = {}
    for site, group in sorted(by_site.items()):
        digests = {b.get("code_digest", "?") for b in group}
        if len(digests) > 1:
            print(f"  WARNING {site} uploaded blocks from {len(digests)} different "
                  f"code versions: {sorted(digests)}")
            for b in group:
                print(f"            {b.get('code_digest', '?')}  {b['_file']}")
        out[site] = group[0]
    return out


def fmt(v) -> str:
    if isinstance(v, dict):
        return " / ".join(str(x) for x in v.values())
    return str(v)


def main(paths: list[str]) -> int:
    bs = blocks(paths)
    if not bs:
        print("no provenance blocks found")
        return 1
    sites = one_per_site(bs)
    print(f"\n{len(bs)} provenance block(s) from {len(sites)} site(s)\n")

    # Widths from the content: a 16-char digest in a 14-char column runs into
    # its neighbour.
    head = ("site", "code_digest", "code_version", "protocol", "definitions")
    rows = [(site, b.get("code_digest", "-"), b.get("code_version", "-"),
             fmt(b.get("protocol_digests", "-")),
             fmt(b.get("definition_versions", "-")))
            for site, b in sites.items()]
    w = [max(len(r[i]) for r in (head, *rows)) + 2 for i in range(len(head))]
    for r in (head, *rows):
        print("".join(c.ljust(n) for c, n in zip(r, w)).rstrip())

    bad = False
    print()
    for key in PROTOCOL_KEYS:
        seen = {site: json.dumps(b.get(key), sort_keys=True) for site, b in sites.items()}
        if len(set(seen.values())) > 1:
            bad = True
            print(f"PROTOCOL MISMATCH on {key} -- a pooled result would not be of "
                  f"one quantity:")
            for site, v in seen.items():
                print(f"    {site:<12} {v}")
        else:
            print(f"protocol OK  {key}: all {len(sites)} site(s) agree")

    # Code divergence: reported against the MODAL digest, never fatal.
    counts = Counter(b.get("code_digest", "-") for b in sites.values())
    modal, n_modal = counts.most_common(1)[0]
    print()
    if len(counts) == 1:
        print(f"code OK      all {len(sites)} site(s) ran {modal}")
    else:
        print(f"CODE DIVERGENCE  {n_modal} of {len(sites)} site(s) ran {modal}")
        ref = next(b for b in sites.values() if b.get("code_digest") == modal)
        ref_files = ref.get("file_digests", {})
        for site, b in sites.items():
            if b.get("code_digest") == modal:
                continue
            mine = b.get("file_digests", {})
            if not mine or not ref_files:
                print(f"    {site}: {b.get('code_digest', '-')} "
                      f"(no file_digests -- cannot say what differs)")
                continue
            diff = sorted(k for k in set(ref_files) | set(mine)
                          if ref_files.get(k) != mine.get(k))
            print(f"    {site}: {b.get('code_digest', '-')} differs in "
                  f"{len(diff)} file(s)")
            for k in diff:
                tag = ("only at this site" if k not in ref_files
                       else "missing at this site" if k not in mine else "edited")
                print(f"        {k}  ({tag})")
        print("    Not fatal -- ask each site what changed and whether it moves a "
              "reported number.")

    # site_config_digest is deliberately NOT compared: it covers config.json,
    # which carries site_name and data_directory and so MUST differ.
    print()
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
