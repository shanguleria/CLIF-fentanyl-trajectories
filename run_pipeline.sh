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
if ! "$PY" tests/test_config_parses.py; then
  echo "ERROR: a tracked config is not valid JSON. Nothing was run." >&2
  exit 1
fi
if ! "$PY" code/check_config.py; then
  echo "ERROR: config is not usable. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_protocol.py; then
  echo "ERROR: config/protocol.json failed its integrity checks. Nothing was run." >&2
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
if ! "$PY" tests/test_paths.py; then
  echo "ERROR: output locations failed their tests. Nothing was run." >&2
  exit 1
fi
# The R side of the same machinery. tests/test_paths.py cannot reach it, and the
# config-drift guard in paths.R was INERT from the day it was written until
# 2026-10-28 -- a changed protocol.json never stopped anything.
if ! Rscript tests/test_require_manifest.R; then
  echo "ERROR: the config-drift guard failed its tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_build_cohort.py; then
  echo "ERROR: cohort-build logic failed its tests. Nothing was run." >&2
  exit 1
fi
if ! "$PY" tests/test_vfd.py; then
  echo "ERROR: ventilator-free days failed its tests. Nothing was run." >&2
  exit 1
fi
echo

echo "== Build cohort (Python) =="
"$PY" code/01_build_cohort.py

# The trajectory-modelling scripts (tabled/gbmt_classes, lcmm_classes,
# transition_model, outcomes) are TABLED as of 2026-09-23, and tabled/landmark_cohort
# joined them on 2026-09-30, so they are deliberately not run here. They carry no
# step number: the number line is this list.
# outcomes.R is an
# unconditional stop(), so with `set -e` its presence in this list made the whole
# pipeline exit non-zero every time -- the documented one-command build could not
# succeed. Run a tabled script by hand if you need it: Rscript code/tabled/gbmt_classes.R
for s in 02_descriptive_cohort 03_exemplar 04_delivery_states 05_titration 06_unit_variation 07_outcome_model; do
  echo "== ${s} (R) =="
  Rscript "code/${s}.R"
done

# tests/test_pooling.py validates PRODUCED OUTPUTS, so it runs here and not in
# the preflight. In the preflight it deadlocked the pipeline two ways: on a
# fresh clone the exports do not exist yet, so the suite failed and nothing ever
# ran; and once a check caught a real defect in a shipped file, the only thing
# that could regenerate that file was the run the check was blocking.
# Failing HERE is the useful place to fail -- the outputs exist and are wrong,
# so the actionable instruction is "do not upload", which is what it says.
echo "== Verifying the shareable outputs =="
# Centre-side maths, run here for the same reason tests/test_pooling.py is:
# a site that never pools still catches drift in the engine it ships.
if ! Rscript tests/test_coefficient_meta.R; then
  echo "ERROR: the inverse-variance pool failed its known-answer tests." >&2
  exit 1
fi
if ! "$PY" tests/test_pooling.py; then
  echo "ERROR: the federated-pooling exports did not pass their checks." >&2
  echo "       Do NOT upload output/final_no_phi/ until this is resolved." >&2
  exit 1
fi
# Same reasoning, same place: tests/test_outcome_model.py reads step 7's output.
# Its load-bearing check is that the proportional-odds test produced a VERDICT --
# on 2026-10-27 it shipped every term NA while the console reported zero
# violations, which is what a pass looks like.
if ! "$PY" tests/test_outcome_model.py; then
  echo "ERROR: the outcome model's outputs did not pass their checks." >&2
  echo "       Do NOT upload output/final_no_phi/ until this is resolved." >&2
  exit 1
fi
echo

echo "Pipeline complete. Review output/final_no_phi/ before sharing."
