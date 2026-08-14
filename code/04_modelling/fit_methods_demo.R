# fit_methods_demo.R
#
# Run EVERY registered susceptibility method (methods_registry.R) on the committed slim panel,
# plot each method's fitted curves season by season, and print the standard cross-method summary
# table (S0, R_eff, reporting c, peak/onset week, fit cor). Adding a method to the registry makes
# it appear here automatically. Base R only -- no data pipeline, no tidyverse.
#
# Run from the repo root:  Rscript code/04_modelling/fit_methods_demo.R

source("code/02_settings/settings_version0.R"); params = settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_sir_deterministic.R")
source("code/01_main_supporting/methods/method_sir_ekf.R")
source("code/01_main_supporting/methods/method_descriptive.R")
source("code/01_main_supporting/methods_registry.R")

# ---- |-settings ----
countries = params$susc_countries                      # the slim panel's countries (settings)
panels    = lapply(countries, load_flu_iliplus_slim); names(panels) = countries
# one column per season so EVERY season is drawn (the demo's whole point). A fixed column count
# smaller than the season count would overflow the mfrow grid and base graphics would silently
# keep only the LAST page of the png.
n_panels  = max(vapply(panels, function(p) length(p$seasons), integer(1)))

# ---- |-one figure + summary table per registered method ----
dir.create("output", showWarnings = FALSE)
for (method in names(sir_methods())){
  png(sprintf("output/fit_%s.png", method), width = 260 * n_panels, height = 140 * length(countries))
  par(mfrow = c(length(countries), n_panels), mar = c(4, 4, 3, 1))

  for (cc in countries){
    country_panel = panels[[cc]]
    elapsed = system.time(fit <- run_method(method, country_panel, params, n_starts = 4))[["elapsed"]]
    pars = fit$params
    qI = if (is.na(pars$qI[1])) "-" else sprintf("%.1e", pars$qI[1])
    cat(sprintf("\n[%s] %s  (R0=%.2f, I0=%.0e)  shared: b=%.1f  phi=%.1f  qI=%s  | conv=%d  %.0fs\n",
                method, cc, fit$R0, fit$seed_i0, pars$b, pars$phi, qI, fit$convergence, elapsed))
    print(summarise_method_fit(fit)[, c("season", "S0", "R_eff", "c", "peak_week", "onset_week", "cor")],
          row.names = FALSE, digits = 3)

    for (s in seq_len(n_panels)){                      # fill the whole row: season panel or blank
      if (s <= length(country_panel$seasons)){
        y = country_panel$ylist[[s]]; wk = country_panel$season_week[[s]]; mu = fit$mu[[s]]
        plot(wk, y, pch = 19, col = "grey30", xlab = "season week", ylab = "flu ILI+",
             main = sprintf("%s %s  S0=%.2f", cc, country_panel$seasons[s], pars$S0[s]))
        lines(wk, mu, col = "red", lwd = 2)             # fitted curve (red), as plotted per method
        abline(h = pars$b, col = "grey60", lty = 3)     # shared off-season baseline
      } else plot.new()
    }
  }
  dev.off()
  cat(sprintf("figure written to output/fit_%s.png\n", method))
}
