"""Between-group dispersion, in two families that must not be conflated.

`overdispersion()` in code/utils/pooling.R returns both a FUNNEL family (null
binomial variances, the scale the funnel envelope is drawn on) and a
META-ANALYTIC family (each group's own variance, prefixed `ma_`). They share most
of their algebra, which is why they are tested apart: this repo has already
shipped one defect of exactly that shape -- a within-episode design effect
applied under a citation for between-unit over-dispersion.

There is no metafor in renv.lock and adding it would breach the pinned-stack
discipline, so the expected values here are computed **longhand in Python** from
the published formulas. That makes this a cross-implementation check rather than
the function agreeing with itself.

Run standalone:  .venv/bin/python tests/test_heterogeneity.py
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
OD_CSV = REPO / "output" / "final_no_phi" / "06_unit_variation" / "overdispersion.csv"


def call_r(pct, n, p0, trim=0.1, min_k=5):
    """Run overdispersion() in R and return its single row as a dict."""
    vec = lambda v: "c(" + ", ".join(repr(float(x)) for x in v) + ")" if v else "numeric(0)"
    script = f"""
      suppressMessages(library(here)); suppressMessages(library(jsonlite))
      source(here("code", "utils", "pooling.R"))
      r <- overdispersion({vec(pct)}, {vec(n)}, {p0!r}, {trim!r}, {min_k}L)
      cat(toJSON(as.list(r), auto_unbox = TRUE, na = "null", digits = 12))
    """
    r = subprocess.run(["Rscript", "-e", script], cwd=REPO,
                       capture_output=True, text=True)
    assert r.returncode == 0, f"Rscript failed:\n{r.stderr[-800:]}"
    return json.loads(r.stdout[r.stdout.index("{"):])


# ------------------------------------------------------------ expected, longhand
def expect_meta(pct, n):
    """Fixed-effect Q, I^2 and DL tau^2 from the published formulas, by hand."""
    v = [p * (100 - p) / ni for p, ni in zip(pct, n)]
    w = [1 / vi for vi in v]
    theta = sum(wi * p for wi, p in zip(w, pct)) / sum(w)
    q = sum(wi * (p - theta) ** 2 for wi, p in zip(w, pct))
    df = len(pct) - 1
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    tau2 = max(0.0, (q - df) / (sum(w) - sum(x * x for x in w) / sum(w)))
    return {"theta_fe": theta, "q": q, "df": df, "i2": i2, "tau2": tau2}


# --------------------------------------------------------------- the known case
#
# Two groups, 10% and 30% of 100. Chosen because every quantity is exact:
#   v  = (9, 21);  w = (1/9, 1/21);  theta_fe = 16.0 exactly
#   Q  = (1/9)(36) + (1/21)(196) = 4 + 9.3333... = 13.3333...
#   I2 = (Q - 1)/Q = 0.925 exactly
KNOWN_PCT, KNOWN_N = [10.0, 30.0], [100.0, 100.0]


def test_the_known_two_group_case_matches_hand_arithmetic():
    got = call_r(KNOWN_PCT, KNOWN_N, p0=20.0)
    exp = expect_meta(KNOWN_PCT, KNOWN_N)
    assert exp["theta_fe"] == pytest.approx(16.0, abs=1e-9), "the hand case drifted"
    assert exp["q"] == pytest.approx(40 / 3, abs=1e-9)
    assert exp["i2"] == pytest.approx(0.925, abs=1e-9)
    assert got["ma_q"] == pytest.approx(exp["q"], rel=1e-9)
    assert got["ma_i2"] == pytest.approx(exp["i2"], rel=1e-9)
    assert got["ma_tau2"] == pytest.approx(exp["tau2"], rel=1e-9)


def test_homogeneous_groups_have_no_heterogeneity():
    """Identical rates must give Q = 0, I2 = 0, tau2 = 0 -- the floor case. A
    function that cannot return zero heterogeneity will never report it."""
    got = call_r([10.0] * 5, [100.0] * 5, p0=10.0)
    assert got["ma_q"] == pytest.approx(0.0, abs=1e-9)
    assert got["ma_i2"] == pytest.approx(0.0, abs=1e-12)
    assert got["ma_tau2"] == pytest.approx(0.0, abs=1e-12)
    assert got["ma_re_pct"] == pytest.approx(10.0, abs=1e-9), (
        "with no heterogeneity the random-effects estimate is the common rate")


def test_the_i2_h2_and_q_identities_hold():
    """Three routes to I^2 must agree, or one of them is wrong."""
    for pct, n in (([5.0, 12.0, 30.0, 44.0], [200.0, 900.0, 150.0, 400.0]),
                   ([2.0, 3.0, 2.5], [1000.0, 1200.0, 800.0])):
        g = call_r(pct, n, p0=sum(pct) / len(pct))
        df_ma = g["ma_k"] - 1
        assert g["ma_i2"] == pytest.approx(max(0.0, (g["ma_q"] - df_ma) / g["ma_q"]), rel=1e-10)
        assert g["ma_h2"] == pytest.approx(g["ma_q"] / df_ma, rel=1e-10)
        h2 = g["ma_h2"]
        if h2 > 1:
            assert g["ma_i2"] == pytest.approx((h2 - 1) / h2, rel=1e-10), (
                "I2 = (H2-1)/H2 and (Q-df)/Q must be the same number")
        assert 0.0 <= g["ma_i2"] <= 1.0
        assert g["ma_i2_lo"] <= g["ma_i2"] <= g["ma_i2_hi"], (
            "the point estimate must sit inside its own interval")
        assert min(pct) <= g["ma_re_pct"] <= max(pct), (
            "a weighted mean cannot fall outside the values it averages")


def test_tau2_is_zero_whenever_q_does_not_exceed_df():
    """The DL estimator is truncated at zero; a negative tau^2 is not a variance."""
    g = call_r([9.9, 10.0, 10.1], [100.0, 100.0, 100.0], p0=10.0)
    assert g["ma_q"] <= g["ma_k"] - 1 + 1e-9
    assert g["ma_tau2"] == pytest.approx(0.0, abs=1e-12)
    assert g["ma_i2"] == pytest.approx(0.0, abs=1e-12)


def test_hartung_knapp_is_wider_than_dersimonian_laird_at_small_k():
    """The whole reason HKSJ is here. At k = 3 the t multiplier is 4.30 against
    the normal's 1.96, and that difference is the honesty."""
    g = call_r([5.0, 20.0, 40.0], [300.0, 300.0, 300.0], p0=21.7)
    hk = g["ma_re_hi"] - g["ma_re_lo"]
    dl = g["ma_dl_hi"] - g["ma_dl_lo"]
    assert hk > dl, f"HKSJ {hk:.3f} should exceed DL {dl:.3f} at k = 3"
    assert g["ma_re_lo"] < g["ma_re_pct"] < g["ma_re_hi"]


def test_reportable_flips_exactly_at_the_declared_k():
    """Computed always, reportable above a threshold. The numbers are never
    hidden -- only their use is governed -- so this flag is the whole gate."""
    for k in (2, 3, 4, 5, 6):
        g = call_r([10.0 + i for i in range(k)], [200.0] * k, p0=12.0, min_k=5)
        assert g["k"] == k
        assert g["reportable"] is (k >= 5), f"k={k} reportable={g['reportable']}"
        assert g["ma_i2"] is not None, (
            f"k={k} must still COMPUTE I2 -- gating reporting is not gating "
            f"computation, or we cannot watch it stabilise as sites arrive")


def test_a_single_group_estimates_nothing():
    """k = 1 means df = 0. Returning a confident-looking number there is the
    failure this guards."""
    g = call_r([10.0], [100.0], p0=10.0)
    assert g["k"] == 1
    assert g["reportable"] is False
    for key in ("q_stat", "ma_q", "ma_i2", "ma_tau2", "ma_re_pct"):
        assert g[key] is None, f"{key} must be NA at k = 1, got {g[key]}"


def test_i2_is_not_share_not_sampling():
    """They are different statistics and the names must not be swapped.

    share_not_sampling is tau2/(tau2 + mean(v)) on NULL variances -- the funnel's
    own question. I^2 is (Q-df)/Q on each group's OWN variance. Reporting the
    former as I^2 in a manuscript would be wrong, so a refactor that quietly
    unified them has to fail here.
    """
    # p0 MUST be the n-weighted mean of the group rates, because that is what it
    # is in real use: 06 passes the overall rate, total events / total n. Took
    # two attempts to get this fixture right -- with an arbitrary p0 the two Qs
    # differ even WHEN conflated (they are centred on different points), so the
    # test passed the broken code twice. Verified by breaking it.
    pct = [5.0, 12.0, 30.0, 44.0]
    n = [200.0, 900.0, 150.0, 400.0]
    p0 = sum(pi * ni for pi, ni in zip(pct, n)) / sum(n)   # 20.5454...
    g = call_r(pct, n, p0=p0)
    assert g["share_not_sampling"] != pytest.approx(g["ma_i2"], abs=1e-6), (
        "share_not_sampling and ma_i2 collapsed to the same value")

    # THE DIRECT SIGNATURE OF THE CONFLATION, and the assertion that actually
    # catches it. The funnel Q uses null variances and is centred on p0; the
    # meta-analytic Q uses each group's own variance and is centred on the
    # weighted mean. If someone computes the ma_ family on null variances the
    # two Qs become the same number -- and the random-effects "estimate" becomes
    # p0 by construction, i.e. the fixed-effect estimate wearing another name.
    #
    # Checked here because the I2-vs-share_not_sampling assertion above is NOT
    # sufficient: those two differ by 0.2-1.7 points even when the families HAVE
    # been conflated, so it passes the broken version. Under conflation the
    # weights become proportional to n, the weighted mean becomes p0 exactly, and
    # ma_q becomes q_stat exactly -- which is why this is the assertion that
    # bites, and why p0 above has to be the n-weighted mean for it to.
    assert g["ma_q"] != pytest.approx(g["q_stat"], rel=1e-9), (
        "ma_q equals the funnel q_stat, so the meta-analytic family is being "
        "computed on NULL variances rather than each group's own")


def test_the_funnel_family_still_describes_the_shipped_figure():
    """The funnel columns are consumed by 06's envelope and captions, so they
    must keep their meaning after the 2026-10-26 lift into utils/."""
    g = call_r([5.0, 12.0, 30.0], [200.0, 900.0, 150.0], p0=14.0)
    assert g["q_df"] == g["k"] - 1
    assert g["phi"] == pytest.approx(g["q_stat"] / g["k"], rel=1e-12), (
        "phi carries divisor k while Q is chi-square on k-1 df; the critical "
        "value must carry the same divisor or the comparison is off by k/(k-1)")
    assert 0.0 <= g["share_not_sampling"] <= 1.0
    assert g["tau2"] >= 0


@pytest.mark.skipif(not OD_CSV.is_file(), reason="06 has not run")
def test_the_shipped_diagnostic_still_carries_the_funnel_family():
    """Guards the rename: the function returns `k`, the file ships `n_units`,
    because in that file the groups are units and the captions say so."""
    import csv
    rows = list(csv.DictReader(OD_CSV.open()))
    assert rows, "overdispersion.csv is empty"
    cols = set(rows[0])
    assert "n_units" in cols and "k" not in cols, (
        "overdispersion.csv must ship n_units; `k` is the shared function's "
        "grouping-agnostic name and is renamed on write")
    assert not [c for c in cols if c.startswith("ma_")], (
        "the meta-analytic family must NOT ship in the within-site file -- a "
        "`reportable` flag about a count of UNITS would be read as a statement "
        "about a count of SITES")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
