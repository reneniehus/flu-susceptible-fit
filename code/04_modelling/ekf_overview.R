# ekf_overview.R
#
# Fit the EKF method across the configured countries (params$susc_countries) and draw a SIDE-BY-SIDE
# overview to see how the fitted curves translate into susceptibility:
#   LEFT  : the fitted EKF curves, one row per country x one column per season (observed points +
#           filtered red curve), exactly the per-method figure but stacked over several countries.
#   RIGHT : the fitted S0 (susceptibility) across countries and seasons -- one line per country plus
#           the cross-country mean -- so the country x season structure is readable at a glance.
# Base R only -- no data pipeline, no tidyverse.
#
# Run from the repo root:  Rscript code/04_modelling/ekf_overview.R

source("code/02_settings/settings_version0.R"); params = settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_sir_deterministic.R")
source("code/01_main_supporting/methods/method_sir_ekf.R")
source("code/01_main_supporting/methods_registry.R")

countries = params$susc_countries
pal = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(countries)]

# ---- |-fit the EKF for every country (independent fits -> run them in parallel) ----
fits = parallel::mclapply(countries, function(cc)
         run_method("ekf", load_flu_iliplus_slim(cc), params, n_starts = 4),
         mc.cores = min(length(countries), parallel::detectCores()))
names(fits) = countries
summ    = do.call(rbind, lapply(fits, summarise_method_fit))
seasons = sort(unique(summ$season))

# ---- |-side-by-side figure: EKF fit grid (left) | S0 summary (right) ----
dir.create("output", showWarnings = FALSE)
png("output/ekf_overview.png", width = 1600, height = 950)
nc = length(countries)
n_panels = 5                                       # season columns drawn per country (LEFT grid); later seasons are not shown
layout(cbind(matrix(1:(nc*n_panels), nrow = nc, byrow = TRUE), nc*n_panels + 1), widths = c(rep(1, n_panels), 2.0))
par(mar = c(3, 3, 2.4, 1), mgp = c(1.8, 0.6, 0))

# LEFT: row per country, column per season -- observed points + fitted EKF (red) curve
for (i in seq_along(countries)){
  fit = fits[[countries[i]]]
  for (j in 1:n_panels){
    if (j <= length(fit$seasons)){
      y = fit$ylist[[j]]; wk = fit$season_week[[j]]; mu = fit$mu[[j]]
      plot(wk, y, pch = 19, cex = 0.5, col = "grey45",
           xlab = "season week", ylab = if (j == 1) paste0(countries[i], "  ILI+") else "",
           main = sprintf("%s %s  S0=%.2f", countries[i], fit$seasons[j], fit$params$S0[j]), cex.main = 0.95)
      lines(wk, mu, col = pal[i], lwd = 2)
    } else plot.new()
  }
}

# RIGHT: fitted S0 across seasons, one line per country + the cross-country mean
par(mar = c(8, 4, 2.4, 1))
plot(NA, xlim = c(1, length(seasons)), ylim = c(0.5, 1), xaxt = "n",
     xlab = "", ylab = "EKF S0 (susceptibility)", main = "Fitted S0 across countries & seasons")
axis(1, at = seq_along(seasons), labels = seasons, las = 2, cex.axis = 0.9)
for (i in seq_along(countries)){
  d = summ[summ$country == countries[i], ]
  lines(match(d$season, seasons), d$S0, col = pal[i], lwd = 2, type = "b", pch = 19)
}
mean_s0 = tapply(summ$S0, factor(summ$season, seasons), mean)
lines(seq_along(seasons), as.numeric(mean_s0), col = "black", lwd = 3, lty = 2)   # cross-country mean
abline(h = 1 / params$susc_R0, col = "grey70", lty = 3)                           # epidemic threshold S0=1/R0
legend("bottomleft", c(countries, "mean", "1/R0"), col = c(pal, "black", "grey70"),
       lwd = c(rep(2, nc), 3, 1), lty = c(rep(1, nc), 2, 3),
       pch = c(rep(19, nc), NA, NA), bty = "n", cex = 0.95)
dev.off()

cat("figure written to output/ekf_overview.png\n")
print(summ[, c("country", "season", "S0", "R_eff", "peak_week", "onset_week", "cor")],
      row.names = FALSE, digits = 3)
