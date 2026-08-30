#!/usr/bin/env Rscript
# Run the established circadian rhythm-detection methods on the simulated series:
# JTK_CYCLE, ARSER and Lomb-Scargle via MetaCycle::meta2d, and RAIN.
.libPaths("~/Rlibs")
suppressPackageStartupMessages({
  library(MetaCycle)
  library(rain)
})

sim_dir <- file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))), "sim")
out_dir <- file.path(dirname(sim_dir), "results")
dir.create(out_dir, showWarnings = FALSE)

sampling_h <- 2
ns <- c(12, 18, 24, 48)

for (n in ns) {
  infile <- file.path(sim_dir, sprintf("sim_n%d.csv", n))
  dat <- read.csv(infile, check.names = FALSE)
  times <- seq(0, by = sampling_h, length.out = n)

  # MetaCycle wants a file; ARS needs evenly spaced data, which we have.
  meta_out <- file.path(out_dir, sprintf("meta_n%d", n))
  dir.create(meta_out, showWarnings = FALSE)
  methods <- if (n >= 24) c("JTK", "LS", "ARS") else c("JTK", "LS")
  suppressWarnings(
    meta2d(infile = infile, filestyle = "csv", outdir = meta_out, timepoints = times,
           minper = 20, maxper = 28, cycMethod = methods, outIntegration = "noIntegration",
           outputFile = TRUE, ARSdefaultPer = 24)
  )

  # RAIN: independent implementation, non-parametric, handles asymmetric waveforms
  mat <- as.matrix(dat[, -1, drop = FALSE])
  rownames(mat) <- dat[[1]]
  r <- rain(t(mat), deltat = sampling_h, period = 24, verbose = FALSE)
  write.csv(data.frame(series_id = rownames(r), rain_p = r$pVal, rain_period = r$period),
            file.path(out_dir, sprintf("rain_n%d.csv", n)), row.names = FALSE)
  cat(sprintf("n=%d done (methods: %s + RAIN)\n", n, paste(methods, collapse = ",")))
}
