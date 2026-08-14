# descriptive_overview.R
#
# Smart side-by-side overview of the descriptive method across the configured countries
# (params$susc_countries):
#   LEFT  : the DK/FR/IE/HU x season grid of smoothed ILI+ curves, with each descriptive feature
#           drawn ON the curve so you can see what it measures -- AUC as the shaded area above the
#           baseline, peak height as a dot, onset week as a vertical line.
#   RIGHT : the four features across seasons, one line per country (no average). AUC and peak height
#           use a LOG y-axis so the very different per-country reporting scales are legible together
#           (a reminder that those two are within-country, across-season comparisons, not cross-country);
#           onset week and steepness are in natural units.
# Smoothing is cheap, so this runs in a second or two. Base R only.
#
# Run from the repo root:  Rscript code/04_modelling/descriptive_overview.R

source("code/02_settings/settings_version0.R"); params = settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_descriptive.R")   # only the method actually run here
source("code/01_main_supporting/methods_registry.R")

countries = params$susc_countries
country_cols = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(countries)]
features = c("auc", "peak_height", "onset_week", "steepness")
feature_labels = c(auc = "AUC (burden, log)", peak_height = "peak height (log)",
         onset_week = "onset week", steepness = "steepness (growth /wk)")
feature_log_axis = c(auc = "y", peak_height = "y", onset_week = "", steepness = "")   # log y for the scale features

# ---- |-smooth every country's seasons and collect the feature table ----
fits = lapply(countries, function(cc) run_method("descriptive", load_flu_iliplus_slim(cc), params))
names(fits) = countries
fit_summary = do.call(rbind, lapply(fits, summarise_method_fit))
seasons = sort(unique(fit_summary$season))

# ---- |-side-by-side figure: feature-annotated curve grid (left) | feature trajectories (right) ----
dir.create("output", showWarnings = FALSE)
n_countries = length(countries)
n_panels = 5                                       # season columns drawn per country (LEFT grid); later seasons are not shown
png("output/descriptive_overview.png", width = 1800, height = 1000)
grid_ids  <- matrix(1:(n_countries*n_panels), n_countries, n_panels, byrow = TRUE)
right_ids <- n_countries*n_panels + seq_len(n_countries)          # right column: one cell per country ROW
layout(cbind(grid_ids, right_ids), widths = c(rep(1, n_panels), 1.7))
# the right column holds the FOUR FEATURE panels stacked into those country-row cells, so it fits
# only while there are at least as many country rows as features -- fail loudly rather than
# silently dropping a feature panel if susc_countries is ever shrunk below 4
stopifnot(length(features) <= n_countries)
par(mar = c(3, 3, 2.2, 1), mgp = c(1.8, 0.6, 0))

# LEFT: row per country, column per season -- smoothed curve with AUC shaded, peak dot, onset line
for (i in seq_len(n_countries)){
  fit = fits[[countries[i]]]
  for (j in 1:n_panels){
    if (j <= length(fit$seasons)){
      y  = fit$ylist[[j]]; wk = fit$season_week[[j]]; mu = fit$mu[[j]]
      baseline = .curve_baseline(mu); finite = is.finite(mu)
      plot(wk, y, type = "n", xlab = "season week",
           ylab = if (j == 1) paste0(countries[i], "  ILI+") else "",
           main = sprintf("%s %s", countries[i], fit$seasons[j]), cex.main = 0.95,
           ylim = range(c(y, mu, baseline), na.rm = TRUE))
      polygon(c(wk[finite], rev(wk[finite])), c(pmax(mu[finite], baseline), rep(baseline, sum(finite))),   # AUC = area above baseline
              col = adjustcolor(country_cols[i], 0.18), border = NA)
      points(wk, y, pch = 19, cex = 0.4, col = "grey55")                     # raw observations
      lines(wk, mu, col = country_cols[i], lwd = 2)                          # smoothed curve
      abline(h = baseline, col = "grey70", lty = 3)                          # baseline
      abline(v = .onset_week(mu, wk, baseline), col = "grey35", lty = 2)     # onset week
      peak_i = which.max(mu); points(wk[peak_i], mu[peak_i], pch = 19, col = country_cols[i], cex = 1.2)  # peak
    } else plot.new()
  }
}

# RIGHT: one panel per feature, one line per country (log y for the scale-dependent features)
par(mar = c(6, 4, 2.2, 1))
for (k in seq_len(n_countries)){
  if (k <= length(features)){
    feature = features[k]
    plot(NA, xlim = c(1, length(seasons)), ylim = range(fit_summary[[feature]], na.rm = TRUE), log = feature_log_axis[[feature]],
         xaxt = "n", xlab = "", ylab = "", main = feature_labels[[feature]], cex.main = 1.05)
    axis(1, at = seq_along(seasons), labels = seasons, las = 2, cex.axis = 0.8)
    for (i in seq_len(n_countries)){
      country_rows = fit_summary[fit_summary$country == countries[i], ]
      lines(match(country_rows$season, seasons), country_rows[[feature]], col = country_cols[i], lwd = 2, type = "b", pch = 19)
    }
    if (k == 1) legend("topright", countries, col = country_cols, lwd = 2, pch = 19, bty = "n", cex = 0.95)
  } else plot.new()
}
dev.off()

cat("figure written to output/descriptive_overview.png\n")
cat("left-panel marks: shaded = AUC (area above baseline), dot = peak, dashed vertical = onset week\n")
print(fit_summary[, c("country", "season", "auc", "peak_height", "onset_week", "steepness")],
      row.names = FALSE, digits = 3)
