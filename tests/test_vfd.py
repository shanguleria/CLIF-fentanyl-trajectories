"""Ventilator-free days, one test per clause of the definition.

The definition is SG's, written out in config/protocol.json under
outcomes.vfd. Every branch below names the clause it covers, because the
clauses are the specification and a VFD implementation that quietly drops one
still returns plausible numbers -- a patient who was never liberated and a
patient whose liberation was discarded both score 0.

Synthetic frames only; no CLIF data is read.

Run standalone:  .venv/bin/python tests/test_vfd.py
"""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import pandas as pd
import pytest

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
_spec = importlib.util.spec_from_file_location("build", REPO / "code" / "01_build_cohort.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)

T0 = pd.Timestamp("2024-01-01 00:00:00", tz="US/Central")
HORIZON_H = B.COV["outcomes"]["vfd"]["horizon_days"] * 24      # 672
SUCC_H = B.COV["outcomes"]["successful_extubation_hours"]      # 48


def days(hours: float) -> float:
    """VFD in days AT THE SHIPPED PRECISION. vfd_28 is rounded to 2 dp on write,
    so comparing to full precision fails on 572/24 = 23.8333 vs 23.83 -- a unit
    mismatch dressed up as a defect. Compare at the precision that ships."""
    return round(hours / 24.0, 2)


def _frames(blocks):
    """blocks: {name: (discharge_hr, [(hr, device[, scaffold]), ...])}.

    A record is REAL unless a third element says otherwise. Eligibility requires
    a real record within cohort.landmark_record_lookback_hours of T, so a fixture
    of scaffold-only rows is excluded -- which is the point of the rule and the
    reason the flag is settable here.
    """
    co, rs = [], []
    for name, (disc_h, recs) in blocks.items():
        co.append({"encounter_block": name, "patient_id": f"p_{name}",
                   "anchor_dttm": T0,
                   "followup_end_dttm": T0 + pd.Timedelta(hours=disc_h),
                   "block_discharge_dttm": T0 + pd.Timedelta(hours=disc_h)})
        for rec in recs:
            hr, dev = rec[0], rec[1]
            rs.append({"encounter_block": name,
                       "recorded_dttm": T0 + pd.Timedelta(hours=hr),
                       "device_category": dev,
                       "is_scaffold": bool(rec[2]) if len(rec) > 2 else False})
    return pd.DataFrame(co), pd.DataFrame(rs)


def _vfd(blocks, deaths=None):
    cohort, resp = _frames(blocks)
    death_hr = pd.Series(deaths or {}, dtype=float).reindex(
        cohort["encounter_block"])
    death_hr.index = cohort["encounter_block"]
    iv = B._imv_intervals(resp, cohort, HORIZON_H)
    got, _ = B.ventilator_free_days(cohort, iv, death_hr)
    return got.set_index("encounter_block")["vfd_28"].to_dict()


# ------------------------------------------------------------------ the clock
def test_the_clock_starts_at_imv_onset_and_stops_at_the_horizon():
    """'t = 0 is the START OF IMV, horizon 28 days.' A liberation that runs past
    day 28 is credited only to day 28, and one that starts before the anchor
    cannot exist because the anchor IS the first IMV record."""
    got = _vfd({
        "late":  (900, [(0, "IMV"), (600, "Nasal Cannula")]),   # 600 -> 672 = 72h
        "early": (900, [(0, "IMV"), (100, "Room Air")]),        # 100 -> 672 = 572h
    })
    assert got["late"] == days(72)
    assert got["early"] == days(572)
    assert max(got.values()) <= B.COV["outcomes"]["vfd"]["horizon_days"]


def test_the_ceiling_is_25_days_because_accrual_starts_at_the_landmark():
    """THE CEILING IS 25, AND THIS TEST HAS BEEN WRONG ONCE ALREADY.

    An earlier version asserted the opposite -- that a liberation inside the
    first 72h should be credited, making the ceiling 28. SG overruled it on
    2026-10-27: "max VFDs allowed by our definition is 25." The ceiling was only
    the visible symptom. The real problem is that free time inside the first 72h
    sits INSIDE THE EXPOSURE WINDOW: a patient off the ventilator for 51 of the
    first 72 hours had far less opportunity to receive a fentanyl titration, so
    exposure and outcome would be mechanically coupled -- which is precisely what
    the landmark design exists to prevent.

    MEASURED at UCMC before the clamp: 124 blocks exceeded 25 days and the
    maximum was 27.32. After: zero and 25.00.
    """
    # Ventilated without interruption to 72h, then free: the ceiling binds.
    got = _vfd({"b": (900, [(0, "IMV"), (72, "Room Air")])})
    assert got["b"] == days(HORIZON_H - 72) == pytest.approx(25.0)

    # Free from hour 4 onward. The clock starts at intubation, so this is 668h of
    # calendar time off the ventilator -- but only the 600h after the landmark
    # accrue. Before the clamp this returned 27.83.
    got = _vfd({"b": (900, [(0, "IMV"), (4, "Room Air")])})
    assert got["b"] == pytest.approx(25.0), "nothing may exceed the 25-day ceiling"

    # Intubated 4h, liberated 68h, reintubated at 72h and ventilated thereafter.
    # The liberation was real and was successful, but all of it predates the
    # landmark, so it accrues nothing.
    got = _vfd({"b": (900, [(0, "IMV"), (4, "Room Air"), (72, "IMV")])})
    assert got["b"] == 0.0, (
        "a liberation wholly inside the exposure window must accrue nothing")


def test_a_spell_qualifies_on_its_full_length_but_accrues_only_after_T():
    """TWO DIFFERENT WINDOWS. A spell QUALIFIES on its full length anywhere in
    [0, horizon]; it ACCRUES only over [T, horizon]. Testing the 48h threshold on
    the clipped length instead would make a successful extubation stop being one
    because of where the landmark happens to fall -- and at UCMC the grid puts
    160 spells across that boundary."""
    # Liberated at hour 40, reintubated at hour 100. The spell is 60h, so it is a
    # successful extubation; 28 of those hours fall after T and are credited.
    got = _vfd({"b": (900, [(0, "IMV"), (40, "Room Air"), (100, "IMV")])})
    assert got["b"] == days(28), (
        "the spell is 60h and qualifies; only its post-landmark hours accrue")

    # Same shape, but the spell is 40h -- a failed extubation. It accrues nothing
    # even though 12h of it fall after the landmark.
    got = _vfd({"b": (900, [(0, "IMV"), (60, "Room Air"), (100, "IMV")])})
    assert got["b"] == 0.0, "a failed extubation accrues nothing on either side of T"


# ------------------------------------------------------- what counts as earned
def test_a_liberation_shorter_than_48h_contributes_nothing():
    """'Extubated and reintubated within 48h -- that is zero VFDs as the
    extubation was unsuccessful.' Not a prorated 30h; zero."""
    got = _vfd({"fail": (900, [(0, "IMV"), (100, "Room Air"), (130, "IMV")]),
                "pass": (900, [(0, "IMV"), (100, "Room Air"), (200, "IMV")])})
    assert got["fail"] == 0.0
    assert got["pass"] == days(100)


def test_the_threshold_is_closed_at_48h():
    """Exactly 48h qualifies; 47.9h does not. Stated as >= in the protocol, and a
    boundary written with the wrong comparator is invisible on real data."""
    got = _vfd({
        "at":    (900, [(0, "IMV"), (100, "Room Air"), (100 + SUCC_H, "IMV")]),
        "under": (900, [(0, "IMV"), (100, "Room Air"), (100 + SUCC_H - 0.1, "IMV")]),
    })
    assert got["at"] == days(SUCC_H)
    assert got["under"] == 0.0


def test_several_qualifying_liberations_each_count():
    """'Multiple qualifying liberations each count.' A reintubation after a
    qualifying period does not forfeit what was already earned."""
    got = _vfd({"b": (900, [(0, "IMV"), (100, "Room Air"), (200, "IMV"),
                            (300, "Nasal Cannula"), (400, "IMV")])})
    assert got["b"] == days(200)


def test_repeated_charting_of_the_same_state_is_one_interval():
    """A liberation is a RUN, not a record count. Charting 'Room Air' hourly for
    100h must score the same as charting it once -- otherwise VFD would track
    documentation frequency, which varies by unit and by site."""
    dense = [(0, "IMV")] + [(h, "Room Air") for h in range(100, 200)] + [(200, "IMV")]
    got = _vfd({"dense": (900, dense),
                "sparse": (900, [(0, "IMV"), (100, "Room Air"), (200, "IMV")])})
    assert got["dense"] == got["sparse"] == days(100)


# --------------------------------------------------------------------- death
def test_death_at_or_before_day_28_scores_zero_whatever_was_accrued():
    """'Death at or before day 28 gives 0 VFDs.' This is the convention that
    sidesteps the competing risk, so a patient with 20 earned days still gets 0."""
    blocks = {"b": (900, [(0, "IMV"), (100, "Room Air")])}
    assert _vfd(blocks)["b"] == days(572)
    assert _vfd(blocks, {"b": 300.0})["b"] == 0.0
    assert _vfd(blocks, {"b": float(HORIZON_H)})["b"] == 0.0, "at the horizon is in"


def test_death_after_day_28_does_not_zero_anything():
    """Dying on day 40 means being ALIVE on day 28, which is what VFD-28 asks.
    Guards a `death_hr.notna()` test standing in for 'died within the horizon'."""
    got = _vfd({"b": (900, [(0, "IMV"), (100, "Room Air")])}, {"b": 960.0})
    assert got["b"] == days(572)


# ------------------------------------------------------- discharged alive
def test_discharge_alive_off_imv_is_assumed_free_through_day_28():
    """'Discharged alive before day 28 and lost to follow-up -> assume alive and
    IMV-free through day 28.'"""
    got = _vfd({"b": (200, [(0, "IMV"), (100, "Room Air")])})
    assert got["b"] == days(HORIZON_H - 100)


def test_discharge_alive_ON_imv_accrues_nothing_further():
    """SG: a trached patient discharged to an LTACH on a ventilator is still on
    IMV. The naive reading of 'assume alive and free' credits them 24 days they
    did not have, which is the single largest way to get VFD wrong."""
    got = _vfd({"ltach": (200, [(0, "IMV")]),
                "ltach_trach": (200, [(0, "IMV"), (150, "IMV")])})
    assert got["ltach"] == 0.0
    assert got["ltach_trach"] == 0.0


def test_a_discharge_on_imv_scores_zero_even_if_days_were_earned_earlier():
    """SG, 2026-10-27: "most of the censored patients left to a SNF or LTACH on a
    vent. Our definition of VFD should count these patients as 0 VFD28 as they
    almost certainly would not liberate that quickly and are lost to follow up."

    REVERSES 0.9.0, which stopped accrual but let them keep what they had earned.
    The two readings agree whenever the patient was never liberated -- which is
    most of them, and why the difference went unnoticed. They diverge exactly
    here: liberated for 100h, reintubated, then sent out on the ventilator.
    """
    got = _vfd({"b": (300, [(0, "IMV"), (100, "Room Air"), (200, "IMV")])})
    assert got["b"] == 0.0, (
        "lost to follow-up on a ventilator scores zero, like death -- it does not "
        "merely stop accruing")

    # The zeroing is conditional on leaving EARLY. A patient still on IMV whose
    # follow-up reaches day 28 is observed, not lost, and keeps what they earned.
    got = _vfd({"b": (900, [(0, "IMV"), (100, "Room Air"), (200, "IMV")])})
    assert got["b"] == days(100), (
        "a patient observed to day 28 is not lost to follow-up; only an EARLY "
        "discharge on the ventilator zeroes the total")


def test_the_extension_is_applied_before_the_48h_filter_not_after():
    """Liberated at hour 640, discharged home at 660, horizon 672. The whole
    liberation is 32h and qualifies for nothing. Extending AFTER the filter would
    credit the 12 remaining hours instead -- a number with no definition behind
    it. Fixed 2026-10-27 before the first run."""
    got = _vfd({"b": (660, [(0, "IMV"), (640, "Room Air")])})
    assert got["b"] == 0.0

    # And the converse: extended, a 20h liberation that reaches the horizon DOES
    # qualify, because 20h + 52h = 72h >= 48h as one uninterrupted period.
    got = _vfd({"b": (620, [(0, "IMV"), (600, "Room Air")])})
    assert got["b"] == days(HORIZON_H - 600)


def test_a_death_discharge_is_never_extended():
    """Extension is conditional on being alive. A patient who dies at hour 200
    must not be handed the remaining 472h because their last device was non-IMV
    -- which is exactly what a terminal extubation looks like in the series."""
    got = _vfd({"b": (200, [(0, "IMV"), (190, "Room Air")])}, {"b": 200.0})
    assert got["b"] == 0.0


def test_only_a_ventilated_discharge_is_zeroed_not_every_early_discharge():
    """The rule keys on the STATE at discharge, not on leaving early. A patient
    extubated and sent home keeps -- and gains -- days; one sent out on the vent
    gets nothing. Conflating the two would zero most of the cohort."""
    got = _vfd({"home":  (200, [(0, "IMV"), (100, "Room Air")]),
                "ltach": (200, [(0, "IMV"), (100, "Room Air"), (150, "IMV")])})
    assert got["home"] == days(HORIZON_H - 100), "off IMV at discharge: extended"
    assert got["ltach"] == 0.0, "on IMV at discharge: zeroed"


# --------------------------------------------------------- tracheostomy
def test_a_tracheostomy_is_on_imv_until_device_category_says_otherwise():
    """'A trached patient counts as ON IMV while device_category is IMV, and
    begins earning VFDs when it switches to non-IMV for more than 48h -- exactly
    the extubation rule.' So trach and extubation must score identically."""
    got = _vfd({"trach": (900, [(0, "IMV"), (240, "IMV"), (300, "Trach Collar")]),
                "extub": (900, [(0, "IMV"), (240, "IMV"), (300, "Room Air")])})
    assert got["trach"] == got["extub"] == days(HORIZON_H - 300)


# ------------------------------------------------- the intervals themselves
def test_intervals_are_half_open_contiguous_and_within_follow_up():
    """_imv_intervals is the single series both eligibility and the outcome now
    read, so its shape is load-bearing: no gaps, no overlaps, nothing past
    follow-up, and alternating states."""
    cohort, resp = _frames({"b": (500, [(0, "IMV"), (100, "Room Air"),
                                        (200, "IMV"), (300, "Room Air")])})
    iv = B._imv_intervals(resp, cohort, HORIZON_H).sort_values("t_start_hr")
    assert list(iv["t_start_hr"]) == [0, 100, 200, 300]
    assert list(iv["t_end_hr"]) == [100, 200, 300, 500]
    assert list(iv["is_imv"]) == [True, False, True, False]
    assert (iv["t_end_hr"] > iv["t_start_hr"]).all(), "zero-length interval"


def test_records_before_the_anchor_and_past_the_horizon_are_dropped():
    """The anchor is the first IMV record, so a negative hour is a data fault;
    a record past the horizon cannot change a 28-day outcome."""
    cohort, resp = _frames({"b": (900, [(0, "IMV"), (100, "Room Air")])})
    extra = pd.DataFrame([
        {"encounter_block": "b", "recorded_dttm": T0 - pd.Timedelta(hours=5),
         "device_category": "Room Air"},
        {"encounter_block": "b", "recorded_dttm": T0 + pd.Timedelta(hours=800),
         "device_category": "IMV"}])
    iv = B._imv_intervals(pd.concat([resp, extra]), cohort, HORIZON_H)
    assert iv["t_start_hr"].min() == 0
    assert iv["t_end_hr"].max() <= HORIZON_H


# =====================================================================
# build_time_to_event -- the liberation and the death are separate facts
# =====================================================================
def _tte(blocks, categories, death_dttm_hr=None):
    """blocks as above; categories: {block: discharge_category}."""
    cohort, resp = _frames(blocks)
    cohort["discharge_category"] = cohort["encounter_block"].map(categories)
    long = pd.DataFrame([
        {"encounter_block": b, "window_start_hr": B.COV["cohort"]["landmark_hours"]
         - B.WINDOW_H, "alive_admitted": True, "imv_status": 1}
        for b in cohort["encounter_block"]])
    patient = pd.DataFrame({
        "patient_id": cohort["patient_id"],
        "death_dttm": [T0 + pd.Timedelta(hours=h) if h is not None else pd.NaT
                       for h in [(death_dttm_hr or {}).get(b)
                                 for b in cohort["encounter_block"]]]})
    cohort = cohort.merge(B.resolve_death(cohort, patient),
                          on="encounter_block", how="left")
    return B.build_time_to_event(cohort, long, resp)[0].set_index("encounter_block")


def test_eligibility_is_the_state_at_T_not_any_record_in_the_window():
    """SG, 2026-10-27: "landmark eligibility must be the state at 72h."

    Until 0.10.0 eligibility read `imv_status == 1` off the window grid, which is
    max(device_category == "IMV") OVER [T-4, T). A patient charted on IMV at hour
    69 and extubated at hour 70 was therefore "ventilated at T" while already off
    the ventilator. MEASURED at UCMC: 238 blocks (3.6% of the old 6,728), and
    they were the source of a VFD-28 ceiling of 27.3 days where the definition
    allows 25.
    """
    T = B.COV["cohort"]["landmark_hours"]
    got = _tte({
        # On IMV at T -- in the cohort.
        "on":  (500, [(0, "IMV"), (T - 3, "IMV")]),
        # IMV charted inside the window, extubated BEFORE T. The old rule let this
        # block in; the state at T is Room Air, so it must now be excluded.
        "off": (500, [(0, "IMV"), (T - 3, "IMV"), (T - 2, "Room Air")]),
    }, {"on": "Home", "off": "Home"})
    assert "on" in got.index
    assert "off" not in got.index, (
        "a block extubated inside the window ending at T is not ventilated AT T")


def test_a_block_interrupted_inside_the_window_is_not_landmark_eligible():
    """SG, 2026-10-30: "in order to be in the landmark cohort you are
    CONTINUOUSLY alive and on IMV for that time window."

    Being on IMV at the INSTANT T was not enough. A patient extubated at hour 12
    and reintubated at hour 60 was eligible despite 48h off the ventilator, and
    did not have the exposure the landmark window exists to guarantee. MEASURED
    at UCMC before the rule: 838 of 5,808 landmark blocks (14.4%) carried some
    charted non-IMV interval inside [0, 72h).
    """
    T = B.COV["cohort"]["landmark_hours"]
    got = _tte({
        "straight": (500, [(0, "IMV"), (T - 3, "IMV")]),
        # Off the ventilator for most of the window, back on before T.
        "gap":      (500, [(0, "IMV"), (12, "Room Air"), (60, "IMV")]),
    }, {"straight": "Home", "gap": "Home"})
    assert "straight" in got.index
    assert "gap" not in got.index, (
        "a block off IMV for 48h inside the window is not continuously ventilated")


def test_an_interruption_within_the_tolerance_keeps_the_block():
    """cohort.landmark_max_interruption_hours. The tolerance is NOT cosmetic:
    clifpy's waterfall fills a NULL device_category but does not denoise a
    CHARTED one, so a single contradictory entry survives as a real interval. At
    UCMC the shortest such interruption is ONE MINUTE, and 140 blocks (2.4% of
    the cohort) would have been excluded on an interruption under an hour."""
    T = B.COV["cohort"]["landmark_hours"]
    tol = float(B.COV["cohort"]["landmark_max_interruption_hours"])
    got = _tte({
        # A blip strictly shorter than the tolerance, then straight back on IMV.
        "blip": (500, [(0, "IMV"), (30, "Nasal Cannula"),
                       (30 + tol / 2, "IMV"), (T - 3, "IMV")]),
        # The same shape, but longer than the tolerance: a real interruption.
        "real": (500, [(0, "IMV"), (30, "Nasal Cannula"),
                       (30 + tol + 1, "IMV"), (T - 3, "IMV")]),
    }, {"blip": "Home", "real": "Home"})
    assert "blip" in got.index, (
        "a sub-tolerance charting blip must not disqualify a ventilated patient")
    assert "real" not in got.index, (
        "an interruption longer than the tolerance is a real break in exposure")


def test_the_state_at_T_must_be_observed_not_forward_filled():
    """cohort.landmark_record_lookback_hours. SG chose a lookback on 2026-10-27
    after seeing what the bare state-at-T rule admitted: 444 blocks whose "on IMV
    at T" was the waterfall forward-filling across a MEDIAN 34.5h of no real
    charting, of which 0.2% were ever liberated and 93.7% scored zero VFDs as
    lost-on-a-ventilator. Two imputations compounding."""
    T = B.COV["cohort"]["landmark_hours"]
    L = B.COV["cohort"]["landmark_record_lookback_hours"]
    got = _tte({
        # Real IMV record comfortably inside the lookback.
        "fresh": (500, [(0, "IMV"), (T - L + 1, "IMV")]),
        # Last REAL record well outside it; the waterfall carries IMV to T, but
        # nothing observed supports that.
        "stale": (500, [(0, "IMV"), (T - L - 1, "IMV")]),
        # Scaffold INSIDE the window does not count as observation -- this is the
        # case that made the old grid rule admit 638 blocks on imputed rows alone.
        "scaffold_only": (500, [(0, "IMV"), (T - L - 1, "IMV"), (T - 1, "IMV", True)]),
    }, {k: "Home" for k in ("fresh", "stale", "scaffold_only")})
    assert "fresh" in got.index
    assert "stale" not in got.index, (
        f"the last real record is more than {L}h before T; the state at T is "
        f"forward-filled, not observed")
    assert "scaffold_only" not in got.index, (
        "a scaffold row is the waterfall's own imputation and cannot be the "
        "evidence that the imputation is correct")


def test_the_lookback_reuses_the_episode_gap_threshold_on_purpose():
    """L and imv_episode_gap_hours answer the same question -- how long may
    charting be silent before the last value stops being believable -- so they
    are the same number by design. If someone changes one, this fails and makes
    them look at the other."""
    assert (B.COV["cohort"]["landmark_record_lookback_hours"]
            == B.COV["cohort"]["imv_episode_gap_hours"]), (
        "the lookback and the episode gap have diverged; that may be right, but "
        "it must be a decision -- see cohort._landmark_record_lookback_why")


def test_a_successful_liberation_survives_a_later_death():
    """THE REGRESSION THIS GUARDS. Before 2026-10-27 extubation_event = 2 erased
    the liberation entirely: MEASURED at UCMC, 399 blocks (5.9% of the landmark
    cohort) were liberated and died a median 9.8 days later, and every one was
    coded 'death' with the extubation discarded. Both facts are now kept."""
    got = _tte({"b": (500, [(0, "IMV"), (71, "IMV"), (100, "Room Air")])}, {"b": "Expired"})
    assert got.loc["b", "extubation_event"] == 1, "the liberation was erased again"
    assert got.loc["b", "mortality_event"] == 1, "the death was lost"
    assert got.loc["b", "vfd_28"] == 0.0, "death inside the horizon still scores 0"


def test_a_death_before_any_qualifying_liberation_is_the_competing_event():
    """Extubated at 100h and dead at 120h: the 20h off the vent never becomes a
    qualifying liberation, so there is nothing to protect and the block is a
    death. The competing-event rule is now implicit in the 48h filter rather
    than a separate `died within succ_h` branch -- one definition, not two."""
    got = _tte({"b": (120, [(0, "IMV"), (71, "IMV"), (100, "Room Air")])}, {"b": "Expired"})
    assert pd.isna(got.loc["b", "ext_hr"])
    assert got.loc["b", "extubation_event"] == 2
    assert got.loc["b", "vfd_28"] == 0.0


def test_the_three_event_codes_are_exhaustive_and_distinct():
    got = _tte({"lib":  (900, [(0, "IMV"), (71, "IMV"), (100, "Room Air")]),
                "dead": (300, [(0, "IMV"), (71, "IMV")]),
                "cens": (300, [(0, "IMV"), (71, "IMV")])},
               {"lib": "Home", "dead": "Expired", "cens": "Still Admitted"})
    assert got.loc["lib", "extubation_event"] == 1
    assert got.loc["dead", "extubation_event"] == 2
    assert got.loc["cens", "extubation_event"] == 0
    assert got.loc["cens", "mortality_event"] == 2, "censored is not dead"


def test_death_dttm_establishes_death_without_a_terminal_discharge():
    """'death_dttm or, if unavailable, discharge_dttm when the category is
    hospice or expired.' MEASURED at UCMC: 3,066 patients carry a death_dttm with
    NO death discharge, so this branch is not hypothetical."""
    got = _tte({"b": (900, [(0, "IMV"), (71, "IMV"), (100, "Room Air")])}, {"b": "Home"},
               death_dttm_hr={"b": 300})
    assert bool(got.loc["b", "died"]) is True
    assert got.loc["b", "mortality_event"] == 1
    assert got.loc["b", "mortality_time"] == pytest.approx(
        300 - B.COV["cohort"]["landmark_hours"])
    assert got.loc["b", "extubation_event"] == 1, "liberated first, still a success"
    assert got.loc["b", "vfd_28"] == 0.0


def test_death_dttm_is_preferred_over_the_discharge_timestamp():
    """Both sources present and disagreeing: death_dttm wins, because it is the
    timestamp and the discharge is only a date of separation."""
    got = _tte({"b": (500, [(0, "IMV"), (71, "IMV")])}, {"b": "Expired"}, death_dttm_hr={"b": 400})
    assert got.loc["b", "death_hr"] == pytest.approx(400.0)


def test_a_site_without_death_dttm_falls_back_to_the_discharge():
    """CLIF 2.1 makes death_dttm optional. A missing COLUMN must not raise -- the
    pipeline ships to sites whose patient table does not carry it."""
    cohort, resp = _frames({"b": (500, [(0, "IMV"), (71, "IMV")])})
    cohort["discharge_category"] = "Expired"
    long = pd.DataFrame([{"encounter_block": "b",
                          "window_start_hr": B.COV["cohort"]["landmark_hours"] - B.WINDOW_H,
                          "alive_admitted": True, "imv_status": 1}])
    cohort = cohort.merge(B.resolve_death(cohort, cohort[["patient_id"]]),
                          on="encounter_block", how="left")
    got = B.build_time_to_event(cohort, long, resp)[0].set_index("encounter_block")
    assert bool(got.loc["b", "died"]) is True
    assert got.loc["b", "death_hr"] == pytest.approx(500.0)


# ------------------------------------------- the exported intervals decompose
def test_the_exported_intervals_sum_to_the_vfd_total_they_describe():
    """ventilator_free_days() now also returns the intervals it summed, so a
    figure can DECOMPOSE VFD rather than re-derive it from ext_hr. If the two
    can drift the figure becomes a second, silently different, estimate."""
    cohort, resp = _frames({
        # liberated at 100h for 100h, reintubated, liberated again at 300h
        "two": (672, [(0, "IMV"), (100, "NIPPV"), (200, "IMV"), (300, "NIPPV")]),
        # a 20h interruption: a FAILED extubation, credited nothing
        "failed": (672, [(0, "IMV"), (100, "NIPPV"), (120, "IMV")]),
        "never": (672, [(0, "IMV")]),
    })
    death_hr = pd.Series(dtype=float).reindex(cohort["encounter_block"])
    death_hr.index = cohort["encounter_block"]
    iv = B._imv_intervals(resp, cohort, HORIZON_H)
    out, annot = B.ventilator_free_days(cohort, iv, death_hr)

    per_block = annot.groupby("encounter_block")["accrued_h"].sum()
    for blk, hours in out.set_index("encounter_block")["vfd_hours"].items():
        assert per_block.get(blk, 0.0) == pytest.approx(hours), blk
    # EVERY interval is exported, ventilated ones included -- a figure has to
    # draw the IMV time too, and only the off-IMV ones carry accrued_h.
    assert annot["is_imv"].any() and (~annot["is_imv"]).any()
    assert (annot.loc[annot["is_imv"], "accrued_h"] == 0).all()
    assert not annot.loc[annot["encounter_block"] == "failed", "qualifies"].any()


# ------------------------------------ the raw-gap rule, applied as a diagnostic
def test_the_raw_gap_rule_closes_a_carried_forward_run():
    """imv_episodes() closes an episode on a gap between RAW IMV records longer
    than imv_episode_gap_hours; _imv_intervals() does not. _close_imv_on_raw_gap
    applies that rule to the intervals so the difference can be MEASURED. It is
    a diagnostic: nothing in the production path calls it."""
    gap = float(B.EPISODE_GAP_H)
    cohort, resp = _frames({"b": (672, [(0, "IMV")])})
    iv = B._imv_intervals(resp, cohort, HORIZON_H)
    assert len(iv) == 1 and iv.iloc[0]["t_end_hr"] == pytest.approx(672.0), \
        "the waterfalled interval should run to follow-up; that is the premise"

    raw = pd.DataFrame({"encounter_block": ["b", "b"], "hr": [0.0, 10.0]})
    closed = B._close_imv_on_raw_gap(iv, raw, gap).sort_values("t_start_hr")
    assert list(closed["is_imv"]) == [True, False]
    assert closed.iloc[0]["t_end_hr"] == pytest.approx(10.0 + gap)
    assert closed.iloc[1]["t_end_hr"] == pytest.approx(672.0)


def test_the_raw_gap_rule_leaves_a_charted_run_alone():
    """A run whose raw records reach its end is not carry-forward, and an IMV
    run with no raw record inside it has no anchor to measure a gap from -- a
    rule that fired there would be inventing an extubation."""
    gap = float(B.EPISODE_GAP_H)
    cohort, resp = _frames({"b": (200, [(0, "IMV")])})
    iv = B._imv_intervals(resp, cohort, HORIZON_H)

    charted = pd.DataFrame({"encounter_block": ["b"], "hr": [200.0 - gap / 2]})
    assert B._close_imv_on_raw_gap(iv, charted, gap)["t_end_hr"].max() == \
        pytest.approx(200.0)

    none_inside = pd.DataFrame({"encounter_block": ["other"], "hr": [1.0]})
    left = B._close_imv_on_raw_gap(iv, none_inside, gap)
    assert len(left) == 1 and left.iloc[0]["t_end_hr"] == pytest.approx(200.0)


def test_the_raw_gap_rule_merges_the_span_it_frees_with_the_next_off_imv_run():
    """A released span butting up against an existing off-IMV interval must
    become ONE spell: ventilator_free_days() qualifies a liberation on a single
    interval's length, so two 30h halves would each fail the 48h test that the
    merged 60h spell passes."""
    gap = float(B.EPISODE_GAP_H)
    cohort, resp = _frames({"b": (672, [(0, "IMV"), (200.0, "NIPPV")])})
    iv = B._imv_intervals(resp, cohort, HORIZON_H)
    raw = pd.DataFrame({"encounter_block": ["b"], "hr": [100.0]})
    closed = B._close_imv_on_raw_gap(iv, raw, gap).sort_values("t_start_hr")

    free = closed[~closed["is_imv"]]
    assert len(free) == 1, "the freed span and the charted NIPPV run must merge"
    assert free.iloc[0]["t_start_hr"] == pytest.approx(100.0 + gap)
    assert free.iloc[0]["t_end_hr"] == pytest.approx(672.0)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
