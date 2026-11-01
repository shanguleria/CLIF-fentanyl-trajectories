# Protocol change history

The changelog that used to live at the top of `config/protocol.json`, where it
was the first thing a reader met and the last thing they needed. Tracked and
shipped beside the protocol, because a site asking "why did this version move"
has the same right to the answer as we do.

`definition_version` is the authority; these entries say what each bump did.

## Version 0.23.0

2026-11-01. FORMATTING ONLY -- no value changed, no key added or removed, key
order preserved, verified by comparing the parsed file against its predecessor
leaf by leaf. Long documentation strings are now stored WRAPPED, as JSON arrays
of lines that rejoin with a single space: JSON has no multi-line string, so a
single `_why` ran to 3,508 characters and the file could only be read by
side-scrolling. Read any documentation value through a join -- `prose()` in
tests/test_protocol.py is the reference. The version changelog moved out of the
protocol and into this file; it was the first thing a reader met and the last
thing they needed. The version still bumps, and every site must still rebuild,
because the file is digested into the manifest and a formatting change
invalidates it exactly as a value change would.

## Version 0.22.0

2026-10-30. PROSE ONLY -- NO DECLARED VALUE CHANGES, no key added or removed,
key order preserved. Measured site results were removed from the _why /
_comment / _DECISION strings: cohort counts, percentages, medians, observed
maxima. Results belong in the output files the pipeline writes, with their
provenance block, where they stay current; a number written into prose is stale
the next time the pipeline runs, and chasing one stale figure through three
documents in a single session is what prompted this. The DECISIONS, their
mechanisms and their warnings are unchanged, and where an output file carries
the measurement the prose now points at it. Published figures from cited
papers, synthetic fixture arithmetic, configured bounds and SG's verbatim
quotes were left intact. Every site must still rebuild: this file is digested
into the manifest, so even a comment edit invalidates it.

## Version 0.21.0

2026-10-30. Adds model.outcome_model.sensitivity_landmark_hours (= [48]) and
titration.adherence_scope_hours (= [24, 48]). 07 refits the PRIMARY model on
the 48h landmark cohort, with exposure windowed to [0,48h) and VFD recomputed
at that landmark, shipped as model = landmark_48h. 24h was swept and DROPPED on
a measured basis -- see _sensitivity_landmark_why. No primary quantity changes;
this is additive.

## Version 0.20.0

2026-10-30. LANDMARK ELIGIBILITY IS NOW CONTINUOUS VENTILATION, not the state
at an instant. cohort.landmark_definition moves state_at_T_observed ->
continuous_imv_to_T_observed and cohort.landmark_max_interruption_hours is
added. A block is landmark-eligible only if it is on IMV at T with an observed
record AND has no charted non-IMV interval longer than the tolerance anywhere
in [0, T). THE LANDMARK COHORT SHRINKS AT EVERY T and every landmark number
moves -- the liberation split, the VFD distribution, Table 1's Landmark column,
the `landmark` rows of the 05 outputs, and the whole of 07. The size is in
07_outcome_model/strobe.csv rather than written here, because it is a site
measurement and this file must read the same at every site.

## Version 0.19.0

2026-10-30. (a) REMOVES cohort.balanced_panel_hours. The balanced-panel
overlays it parameterised were retired from 02_descriptive_cohort.R -- a
trajectory-era figure nothing downstream consumed -- and the key was also being
read into an unused variable in 04_delivery_states.R, so it was dead in both
places. (b) ADDS exposure.sedatives.course_summary_categories: which sedatives
the new per-admission course summary describes (propofol, midazolam), a
deliberate SUBSET of other_sedative_categories -- dexmedetomidine stays loaded
and keeps its window column but is not in that table. (c) Sweeps the dead Phase
1-6 vocabulary out of the prose; the stages are numbered scripts now and some
are retired under code/tabled/. NO MEASURED QUANTITY CHANGES, but every site
must rebuild: this file is digested into the manifest, so even a comment edit
invalidates it.

## What moved out of `config.json` (2026-10-26)

The seven keys below -- cohort, outcomes, medications, model, reporting and
their comments -- lived in config/config.json until 2026-10-26. That file is
NOT tracked in git (.gitignore) and is deliberately excluded from
PROTOCOL_FILES in code/utils/paths.py, because it carries site_name and
data_directory and its digest is EXPECTED to differ between sites. So the
parameters that DEFINE THE COHORT -- imv_episode_gap_hours, min_imv_hours,
min_age, window_hours, landmark_hours, the medication vocabulary, the outcome
categories -- were invisible to the only mechanism that proves two sites ran
the same protocol. A site nudging imv_episode_gap_hours to 4 to make a run work
would build a different cohort and check_provenance.py would still report
agreement. config.json said so itself in its own _comment_covariates --
'nothing estimand-defining may live here alone' -- and then held the estimand.
Moved by SG 2026-10-26 as part of the NU/Rush federation; the file was renamed
covariates.json -> protocol.json in the same change, because it has held cohort
rules, indication, unit_variation, exemplar and disclosure for weeks and the
old name named only a fraction of it.

## Before 0.19.0

`git log --follow config/protocol.json` — the file was `covariates.json`
until 2026-10-26.
