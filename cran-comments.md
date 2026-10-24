# CRAN submission comments — G2Gcov 0.1.0

## R CMD check results

0 errors | 0 warnings | 1 note

The single note is from the CRAN incoming feasibility check:

    Possibly misspelled words in DESCRIPTION:
      Fader (13:37)
      Grassia (13:47)
      Hardie's (13:47)

These are not misspellings. "Grassia" is the surname of the statistician who
introduced the Grassia(II) distribution (Grassia, 1977). "Fader" and
"Hardie's" refer to Peter Fader and Bruce Hardie, the authors of the
methodology this package implements.

## Platform testing

All platforms returned `Status: OK` with 0 errors, 0 warnings, 0 notes.

| Platform                        | R version | Result |
|---------------------------------|-----------|--------|
| Windows 11 x64 (local)          | R 4.5.2   | OK     |
| Linux x86_64 Ubuntu 24.04 LTS   | R-devel   | OK     |
| macOS arm64 (macos-latest)      | R-devel   | OK     |
| Windows (windows-latest)        | R-devel   | OK     |

Linux, macOS, and Windows R-devel checks were run via R-hub v2 GitHub Actions
(`rhub::rhub_check()`). Local Windows check run via `devtools::check()`.

## Notes for CRAN reviewers

* This is a new submission.
* The package implements the Grassia(II)-Geometric (G2G) discrete-time survival
  model with support for both time-invariant and time-varying covariates,
  following the methodology of Fader & Hardie (1997).
* All \donttest{} examples run successfully; they are wrapped in \donttest{}
  because numerical optimization on the bundled kb_data takes up to about a
  minute per call.
* During BFGS optimization, R may internally produce "NaNs produced" warnings
  as the optimizer explores regions of the parameter space where the
  log-likelihood is undefined. These are expected, benign, and are suppressed
  within the package using withCallingHandlers(). They do not appear to end
  users and did not produce any warnings in R CMD check across any platform.
* survival is listed under Suggests and used only in vignettes and tests.
  Model formulas are written as Surv(time, status) ~ ..., but the package
  reads the variable names from the formula without calling Surv(), so the
  core functions work without survival installed.
* The time-varying vignette fits several models and takes a few minutes to
  run, so it is precomputed: its output is generated locally from a source
  file kept in the GitHub repository, and the package builds it without
  refitting any models. The time-invariant vignette is fast and runs live.
* URL and BugReports fields point to https://github.com/kaloklee/G2G_cov.

## Downstream dependencies

None — this is a new package with no reverse dependencies.
