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
source("code/01_main_supporting/methods_registry.R")

# ---- |-settings ----
countries = params$susc_countries                      # the slim panel's countries (settings)
n_panels  = 5                                          # max season panels per country row

# ---- |-one figure + summary table per registered method ----
dir.create("output", showWarnings = FALSE)
for (m in names(sir_methods())){
  png(sprintf("output/fit_%s.png", m), width = 1300, height = 560)
  par(mfrow = c(length(countries), n_panels), mar = c(4, 4, 3, 1))

  for (cc in countries){
    sl = load_flu_iliplus_slim(cc)
    t  = system.time(fit <- run_method(m, sl, params, n_starts = 4))[["elapsed"]]
    p  = fit$params
    qI = if (is.na(p$qI[1])) "-" else sprintf("%.1e", p$qI[1])
    cat(sprintf("\n[%s] %s  (R0=%.2f, I0=%.0e)  shared: b=%.1f  phi=%.1f  qI=%s  | conv=%d  %.0fs\n",
                m, cc, fit$R0, fit$seed_i0, p$b, p$phi, qI, fit$convergence, t))
    print(summarise_method_fit(fit)[, c("season", "S0", "R_eff", "c", "peak_week", "onset_week", "cor")],
          row.names = FALSE, digits = 3)

    for (s in seq_along(sl$seasons)){
      y = sl$ylist[[s]]; wk = sl$season_week[[s]]; mu = fit$mu[[s]]
      plot(wk, y, pch = 19, col = "grey30", xlab = "season week", ylab = "flu ILI+",
           main = sprintf("%s %s  S0=%.2f", cc, sl$seasons[s], p$S0[s]))
      lines(wk, mu, col = "red", lwd = 2)                # fitted curve (red), as plotted per method
      abline(h = p$b, col = "grey60", lty = 3)           # shared off-season baseline
    }
    if (length(sl$seasons) < n_panels)                   # keep the grid aligned across rows
      for (k in seq_len(n_panels - length(sl$seasons))) plot.new()
  }
  dev.off()
  cat(sprintf("figure written to output/fit_%s.png\n", m))
}
