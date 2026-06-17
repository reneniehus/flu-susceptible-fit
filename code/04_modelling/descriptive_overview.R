# descriptive_overview.R
#
# Run the descriptive method across the configured countries (params$susc_countries) and draw a
# side-by-side overview:
#   LEFT  : raw weekly ILI+ (points) + the smoothed curve, one row per country x one column per season.
#   RIGHT : the descriptive summary statistics across seasons -- AUC, peak height, onset week and
#           steepness -- one panel each, one line per country (no average line).
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
flab = c(auc = "AUC (burden)", peak_height = "peak height",
         onset_week = "onset week", steepness = "steepness (growth /wk)")

# ---- |-smooth every country's seasons and collect the feature table ----
fits = lapply(countries, function(cc) run_method("descriptive", load_flu_iliplus_slim(cc), params))
names(fits) = countries
summ    = do.call(rbind, lapply(fits, summarise_method_fit))
seasons = sort(unique(summ$season))

# ---- |-side-by-side figure: smoothed fits (left grid) | feature summaries (right column) ----
dir.create("output", showWarnings = FALSE)
nc = length(countries)
png("output/descriptive_overview.png", width = 1700, height = 950)
layout(cbind(matrix(1:(nc*5), nc, 5, byrow = TRUE), nc*5 + seq_len(nc)), widths = c(rep(1, 5), 1.7))
par(mar = c(3, 3, 2.4, 1), mgp = c(1.8, 0.6, 0))

# LEFT: row per country, column per season -- raw points + smoothed curve
for (i in seq_len(nc)){
  fit = fits[[countries[i]]]
  for (j in 1:5){
    if (j <= length(fit$seasons)){
      y = fit$ylist[[j]]; wk = fit$season_week[[j]]; mu = fit$mu[[j]]
      plot(wk, y, pch = 19, cex = 0.5, col = "grey45",
           xlab = "season week", ylab = if (j == 1) paste0(countries[i], "  ILI+") else "",
           main = sprintf("%s %s", countries[i], fit$seasons[j]), cex.main = 0.95)
      lines(wk, mu, col = pal[i], lwd = 2)              # smoothed curve
    } else plot.new()
  }
}

# RIGHT: one panel per descriptive feature, one line per country (no average line)
par(mar = c(6, 4, 2.4, 1))
for (k in seq_len(nc)){
  if (k <= length(features)){
    ft = features[k]
    plot(NA, xlim = c(1, length(seasons)), ylim = range(summ[[ft]], na.rm = TRUE), xaxt = "n",
         xlab = "", ylab = "", main = flab[[ft]], cex.main = 1.0)
    axis(1, at = seq_along(seasons), labels = seasons, las = 2, cex.axis = 0.75)
    for (i in seq_len(nc)){
      d = summ[summ$country == countries[i], ]
      lines(match(d$season, seasons), d[[ft]], col = pal[i], lwd = 2, type = "b", pch = 19)
    }
    if (k == 1) legend("topleft", countries, col = pal, lwd = 2, pch = 19, bty = "n", cex = 0.9)
  } else plot.new()
}
dev.off()

cat("figure written to output/descriptive_overview.png\n")
print(summ[, c("country", "season", "auc", "peak_height", "onset_week", "steepness", "S0")],
      row.names = FALSE, digits = 3)
