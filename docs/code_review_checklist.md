# Code review checklist

**Any generated code gets read by a human before it is run or committed.** Code
you did not read is code you cannot defend to a co-author or a reviewer.

Work through this for every `code/*.R` and `code/*.py` implementation.

## Readability

- [ ] Could I explain every line to a co-author without re-deriving it?
- [ ] Do comments explain **why**, not **what**? Delete any comment that
      restates the code (`# loop over patients` above a `for` loop).
- [ ] Are variable names the ones a person would choose? Short and concrete
      beats long and descriptive-of-its-own-type.
- [ ] Does it match the house style in `~/Desktop/Research/CLIF/R_setup.R` —
      banner, package loop, `here()` paths, `sessionInfo()` at the end?

## LLM fluff to cut

- [ ] **Over-commenting.** One comment per logical block, not per line.
- [ ] **Defensive branches for impossible conditions.** If the input cannot be
      `NULL` there, do not check for `NULL`.
- [ ] **Helpers used once.** Inline them unless the name genuinely earns its keep.
- [ ] **Boilerplate `tryCatch` / `try` that swallows errors.** A pipeline should
      fail loudly. Only catch what you can actually handle.
- [ ] **Restating the config in comments.** A comment claiming "same as the
      config" is a claim to verify, not documentation. Read the file instead.
- [ ] **Unused arguments, dead parameters, speculative generality.**

## Correctness traps specific to this project

- [ ] Infusion uses **LOCF then a time-weighted mean**; boluses are **summed and
      never carried forward**. Check the rule is applied per column.
- [ ] `scaling = 0` reaches `gbmt()`. Print the effective value, do not assume.
- [ ] Every config key is actually consumed. A key that no code path reaches is
      worse than none: it reads as policy, so the next reader trusts it.
- [ ] Dose in extubated windows is `0`, driven by `imv_status`.
- [ ] Nothing patient-level is written outside `data/intermediate_phi/` or
      `output/intermediate_phi/`. Anything in `output/final_no_phi/` is
      aggregate and safe to share.
- [ ] Effective `d` was not silently capped (`gbmt` warns rather than errors).

## Before committing

- [ ] Script runs top-to-bottom in a **fresh** R session.
- [ ] `sessionInfo()` log written.
- [ ] No absolute paths, no `setwd()`, no `rm(list = ls())`.
- [ ] `git status` shows nothing under `data/`, `output/*/intermediate_phi/`, or
      `logs/`.
