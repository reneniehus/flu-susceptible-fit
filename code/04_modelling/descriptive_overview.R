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
source("code/01_main_supporting/methods/method_sir_deterministic.R")
source("code/01_main_supporting/methods/method_sir_ekf.R")
source("code/01_main_supporting/methods/method_descriptive.R")
source("code/01_main_supporting/methods_registry.R")

countries = params$susc_countries
pal = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(countries)]
features = c("auc", "peak_height", "onset_week", "steepness")
flab = c(auc = "AUC (burden, log)", peak_height = "peak height (log)",
         onset_week = "onset week", steepness = "steepness (growth /wk)")
flog = c(auc = "y", peak_height = "y", onset_week = "", steepness = "")   # log y for the scale features

# ---- |-smooth every country's seasons and collect the feature table ----
fits = lapply(countries, function(cc) run_method("descriptive", load_flu_iliplus_slim(cc), params))
names(fits) = countries
summ    = do.call(rbind, lapply(fits, summarise_method_fit))
seasons = sort(unique(summ$season))

# ---- |-side-by-side figure: feature-annotated curve grid (left) | feature trajectories (right) ----
dir.create("output", showWarnings = FALSE)
nc = length(countries)
n_panels = 5                                       # season columns drawn per country (LEFT grid); later seasons are not shown
png("output/descriptive_overview.png", width = 1800, height = 1000)
layout(cbind(matrix(1:(nc*n_panels), nc, n_panels, byrow = TRUE), nc*n_panels + seq_len(nc)), widths = c(rep(1, n_panels), 1.7))
par(mar = c(3, 3, 2.2, 1), mgp = c(1.8, 0.6, 0))

# LEFT: row per country, column per season -- smoothed curve with AUC shaded, peak dot, onset line
for (i in seq_len(nc)){
  fit = fits[[countries[i]]]
  for (j in 1:n_panels){
    if (j <= length(fit$seasons)){
      y  = fit$ylist[[j]]; wk = fit$season_week[[j]]; mu = fit$mu[[j]]
      b  = .curve_baseline(mu); ok = is.finite(mu)
      plot(wk, y, type = "n", xlab = "season week",
           ylab = if (j == 1) paste0(countries[i], "  ILI+") else "",
           main = sprintf("%s %s", countries[i], fit$seasons[j]), cex.main = 0.95,
           ylim = range(c(y, mu, b), na.rm = TRUE))
      polygon(c(wk[ok], rev(wk[ok])), c(pmax(mu[ok], b), rep(b, sum(ok))),   # AUC = area above baseline
              col = adjustcolor(pal[i], 0.18), border = NA)
      points(wk, y, pch = 19, cex = 0.4, col = "grey55")                     # raw observations
      lines(wk, mu, col = pal[i], lwd = 2)                                   # smoothed curve
      abline(h = b, col = "grey70", lty = 3)                                 # baseline
      abline(v = .onset_week(mu, wk, b), col = "grey35", lty = 2)            # onset week
      pk = which.max(mu); points(wk[pk], mu[pk], pch = 19, col = pal[i], cex = 1.2)  # peak
    } else plot.new()
  }
}

# RIGHT: one panel per feature, one line per country (log y for the scale-dependent features)
par(mar = c(6, 4, 2.2, 1))
for (k in seq_len(nc)){
  if (k <= length(features)){
    ft = features[k]
    plot(NA, xlim = c(1, length(seasons)), ylim = range(summ[[ft]], na.rm = TRUE), log = flog[[ft]],
         xaxt = "n", xlab = "", ylab = "", main = flab[[ft]], cex.main = 1.05)
    axis(1, at = seq_along(seasons), labels = seasons, las = 2, cex.axis = 0.8)
    for (i in seq_len(nc)){
      d = summ[summ$country == countries[i], ]
      lines(match(d$season, seasons), d[[ft]], col = pal[i], lwd = 2, type = "b", pch = 19)
    }
    if (k == 1) legend("topright", countries, col = pal, lwd = 2, pch = 19, bty = "n", cex = 0.95)
  } else plot.new()
}
dev.off()

cat("figure written to output/descriptive_overview.png\n")
cat("left-panel marks: shaded = AUC (area above baseline), dot = peak, dashed vertical = onset week\n")
print(summ[, c("country", "season", "auc", "peak_height", "onset_week", "steepness")],
      row.names = FALSE, digits = 3)
