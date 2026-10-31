# Disclosure rules on Table 1's tidy categorical export, in code/utils/table1.R.
#
# Every scenario here is built by hand because the branches that matter do not
# fire on UCMC's data: no race level at the `overall` scope is below the
# threshold, so the cross-view rule is inert at the only scope that publishes
# both views. It is live at a smaller site, which is where it has to be right.
#
# `recoverable()` is the attack, not the implementation: it re-derives the two
# routes a withheld cell comes back by -- the row's own denominator, and a
# published collapsed level whose sources are published but one -- and every
# scenario is held to returning nothing.
#
# Run standalone:  Rscript tests/test_table1_disclosure.R
suppressPackageStartupMessages(library(here))
source(here("code", "utils", "pooling.R"))
source(here("code", "utils", "table1.R"))

pass <- 0L; fail <- 0L
check <- function(name, ok, detail = "") {
  if (isTRUE(ok)) { pass <<- pass + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else { fail <<- fail + 1L; cat(sprintf("  FAIL  %s\n        %s\n", name, detail)) }
}

MAP <- list(A = "G1", B = "G1", C = "G2", D = "G2", E = "G3", F = "G3",
            White = "White", Unknown = "Unknown", `_default` = "G1")
MIN <- 10L

# A scope is a named vector of level -> count, expanded to the raw vector the
# builders take.
vec <- function(...) {
  k <- c(...)
  rep(names(k), unname(k))
}
emit <- function(scopes, full_scope = "all")
  t1_categorical("V", scopes, collapse = MAP, min_cell = MIN,
                 scope_order = names(scopes), full_scope = full_scope)$pool_cat

held <- function(d, var, st) sort(d$level[d$variable == var & d$stratum == st &
                                            d$n_suppressed_small_cell == 1L])
shown_n <- function(d, var, st, lv) d$n[d$variable == var & d$stratum == st &
                                          d$level == lv]

recoverable <- function(d) {
  bad <- character(0)
  for (k in split(d, list(d$variable, d$stratum), drop = TRUE)) {
    w <- k$n_suppressed_small_cell == 1L
    if (sum(w) == 1L)
      bad <- c(bad, sprintf("%s/%s/%s from its denominator",
                            k$variable[1], k$stratum[1], k$level[w]))
  }
  for (st in unique(d$stratum)) {
    src <- d[d$variable == "V" & d$stratum == st, , drop = FALSE]
    col <- d[d$variable == "V (collapsed)" & d$stratum == st, , drop = FALSE]
    if (!nrow(src) || !nrow(col)) next
    tgt <- collapse_levels(src$level, map = MAP)
    s_all <- sum(src$n_suppressed_small_cell == 1L)
    for (i in seq_len(nrow(col))) {
      if (col$n_suppressed_small_cell[i] == 1L) next
      m <- sum(src$n_suppressed_small_cell[tgt == col$level[i]] == 1L)
      if (m == 1L)
        bad <- c(bad, sprintf("V/%s/%s from published %s", st,
                              src$level[tgt == col$level[i] &
                                          src$n_suppressed_small_cell == 1L],
                              col$level[i]))
      if (m >= 1L && s_all - m == 1L)
        bad <- c(bad, sprintf("V/%s the one source outside %s", st, col$level[i]))
    }
  }
  bad
}
check_safe <- function(name, d) {
  bad <- recoverable(d)
  check(name, !length(bad), paste(bad, collapse = "; "))
}

# ---- 1. all the suppressed sources inside ONE collapsed level ---------------
# The 2026-10-31 regression. Withholding G1 told a lie -- the collapsed row sums
# to its denominator, so G1 came straight back -- while publishing it reveals
# only A + B, which the full row's own denominator already gives.
s1 <- list(overall = vec(A = 5, B = 6, White = 100, Unknown = 50))
d1 <- emit(s1)
check("a collapsed level holding every suppressed source is published",
      identical(held(d1, "V (collapsed)", "overall"), character(0)) &&
        shown_n(d1, "V (collapsed)", "overall", "G1") == 11,
      paste(held(d1, "V (collapsed)", "overall"), collapse = ","))
check("and its sources stay withheld",
      identical(held(d1, "V", "overall"), c("A", "B")))
check_safe("scenario 1 leaks nothing", d1)

# ---- 2. one suppressed source inside, one outside ---------------------------
# Both collapsed levels isolate a cell, so both go. Their sum is then all the
# denominator gives, and that sum is public already.
s2 <- list(overall = vec(A = 5, B = 40, C = 6, D = 40, White = 100, Unknown = 50))
d2 <- emit(s2)
check("a collapsed level with exactly one suppressed source is withheld",
      identical(held(d2, "V (collapsed)", "overall"), c("G1", "G2")),
      paste(held(d2, "V (collapsed)", "overall"), collapse = ","))
check_safe("scenario 2 leaks nothing", d2)

# ---- 3. the min_suppressed_cells_per_row partner ----------------------------
# G1 isolates A and goes; G2 and G3 each hold two, so the rule alone would leave
# exactly one withheld cell in the row -- recoverable by subtraction. The
# partner is the smallest level that itself holds two or more.
s3 <- list(overall = vec(A = 5, B = 40, C = 5, D = 6, G = 40, E = 5, F = 6,
                         H = 50, White = 100))
d3 <- emit(s3)
check("a lone withheld cell takes the smallest two-source level with it",
      identical(held(d3, "V (collapsed)", "overall"), c("G1", "G2")),
      paste(held(d3, "V (collapsed)", "overall"), collapse = ","))
check_safe("scenario 3 leaks nothing", d3)

# ---- 4. nothing small: nothing withheld, both views ship --------------------
s4 <- list(overall = vec(A = 40, C = 40, White = 100, Unknown = 50))
d4 <- emit(s4)
check("no suppression means no withholding",
      !any(d4$n_suppressed_small_cell == 1L) && !anyNA(d4$n))
check("both views ship when full_scope is 'all'",
      setequal(unique(d4$variable), c("V", "V (collapsed)")))

# ---- 5. full_levels_strata_scope ---------------------------------------------
s5 <- list(overall = vec(A = 5, B = 40, C = 6, D = 40, White = 100),
           sub     = vec(A = 5, B = 40, C = 6, D = 40, White = 100))
d5 <- emit(s5, full_scope = "overall")
check("the full view ships in the named scope only",
      identical(sort(unique(d5$stratum[d5$variable == "V"])), "overall"),
      paste(unique(d5$stratum[d5$variable == "V"]), collapse = ","))
check("the collapsed view still ships in every scope",
      setequal(unique(d5$stratum[d5$variable == "V (collapsed)"]),
               c("overall", "sub")))
check("a scope with no full view withholds nothing across views",
      identical(held(d5, "V (collapsed)", "sub"), character(0)) &&
        shown_n(d5, "V (collapsed)", "sub", "G1") == 45)
check("while the scope that publishes both still does",
      identical(held(d5, "V (collapsed)", "overall"), c("G1", "G2")))
check_safe("scenario 5 leaks nothing", d5)
d5all <- emit(s5, full_scope = "all")
check("flipping the lever to 'all' restores the per-scope full view",
      setequal(unique(d5all$stratum[d5all$variable == "V"]), c("overall", "sub")))
check_safe("scenario 5 leaks nothing with the lever at 'all'", d5all)

# ---- 6. the lever must reach an operation ------------------------------------
check("a scope name that no scope carries is refused",
      inherits(try(emit(s5, full_scope = "nope"), silent = TRUE), "try-error"))
check("a full_scope on a variable with no collapsed twin is refused",
      inherits(try(t1_categorical("W", list(overall = c("a", "b")),
                                  min_cell = MIN, full_scope = "overall"),
                   silent = TRUE), "try-error"))

cat(sprintf("\n%d/%d passed\n", pass, pass + fail))
if (fail) quit(status = 1)
