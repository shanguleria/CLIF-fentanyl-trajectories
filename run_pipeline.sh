#!/usr/bin/env bash
# Run the full pipeline. Stops on first error.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f config/config.json ]; then
  echo "ERROR: config/config.json not found. Run:" >&2
  echo "  cp config/config_template.json config/config.json" >&2
  echo "then edit it for your site." >&2
  exit 1
fi

PY="${PYTHON:-python3}"
[ -x .venv/bin/python ] && PY=.venv/bin/python

# Preflight. Fails here, in a second, rather than an hour into a table read.
echo "== Preflight =="
if ! "$PY" code/check_config.py; then
  echo "ERROR: config is not usable. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_covariates.py; then
  echo "ERROR: config/covariates.json failed its integrity checks. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_fio2.py; then
  echo "ERROR: fio2 unit handling failed its tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_outliers.py; then
  echo "ERROR: outlier bounds failed their tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_waterfall_cache.py; then
  echo "ERROR: the waterfall cache failed its tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_doses.py; then
  echo "ERROR: dose unit conversion failed its tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_pooling.py; then
  echo "ERROR: the federated-pooling exports are not poolable. Nothing was run." >&2
  exit 1
fi

if ! "$PY" tests/test_paths.py; then
  echo "ERROR: output locations failed their tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_build_cohort.py; then
  echo "ERROR: Phase 0 logic failed its tests. Nothing was run." >&2
  exit 1
fi
echo

echo "== Phase 0: build cohort (Python) =="
"$PY" code/01_build_cohort.py

for s in 02_descriptive_trajectory 03_landmark_cohort 04_gbmt_classes \
         05_lcmm_classes 06_two_indicator 07_outcomes; do
  echo "== ${s} (R) =="
  Rscript "code/${s}.R"
done

echo "Pipeline complete. Review output/final_no_phi/ before sharing."
