# Regenerate the precomputed time-varying vignette.
#
# time-varying-introduction.Rmd.orig is the live source. Knitting it runs all
# of its code (several model fits, a few minutes) and writes
# time-varying-introduction.Rmd with the printed output already filled in,
# plus its figures in vignettes/figures-tv/. CRAN then builds the vignette
# from that file without refitting any models.
#
# Re-run this after any change to the package code or to the .Rmd.orig file,
# from the package root:
#
#   source("vignettes/precompile.R")

# The vignette calls library(G2Gcov), so install the current source first.
devtools::install(upgrade = "never", build_vignettes = FALSE, quiet = TRUE)

local({
  old_wd <- setwd("vignettes")
  on.exit(setwd(old_wd))
  unlink("figures-tv", recursive = TRUE)
  knitr::knit("time-varying-introduction.Rmd.orig",
              "time-varying-introduction.Rmd",
              envir = new.env())
})
