# fit_kalman_sir_multiseason_demo.R
#
# Demo: the deterministic per-season SUSCEPTIBILITY fit (code/01_main_supporting/model_kalman_sir.R,
# fit_sir_susceptibility). R0, the infectious period and the seed I0 are FIXED (from settings); for
# each season the model fits the susceptibility S0 and the reporting fraction c, with a baseline b
# and overdispersion phi shared across that country's seasons. No Kalman filter -- each season's S0
# is identified by the wave's rise rate, so the headline output is the per-season S0 (read as
# RELATIVE susceptibility across seasons).
#
# Runs on the committed slim panel (data/slim_flu_iliplus.csv: DK & FR flu sentinel ILI+, each
# seeded at week 1 = season start), so it needs only base R -- no data pipeline, no tidyverse.
#
# Run from the repo root:  Rscript code/04_modelling/fit_kalman_sir_multiseason_demo.R

source("code/02_settings/settings_version0.R"); params <- settings()
source("code/01_main_supporting/model_kalman_sir.R")

# ---- |-settings ----
countries = c("DK", "FR")
n_panels  = 5                                          # max season panels per country row

# ---- |-fit each country across its seasons (shared b, phi), then plot every season ----
dir.create("output", showWarnings = FALSE)
png("output/kalman_sir_multiseason_fit.png", width = 1300, height = 560)
par(mfrow = c(length(countries), n_panels), mar = c(4, 4, 3, 1))

for (cc in countries){
  sl = load_flu_iliplus_slim(cc)
  t  = system.time(
    fit <- fit_sir_susceptibility(sl$ylist, R0 = params$susc_R0,
                                  infectious_period_days = params$susc_infectious_period_days,
                                  seed_i0 = params$susc_seed_i0, n_starts = 4))[["elapsed"]]
  p = fit$params
  cat(sprintf("\n%s  (R0=%.2f, I0=%.0e, %d seasons -> %d params)  shared: b=%.1f  phi=%.1f  | conv=%d  %.0fs\n",
              cc, fit$R0, fit$seed_i0, length(sl$seasons), 2*length(sl$seasons) + 2, p$b, p$phi,
              fit$convergence, t))
  ord = order(p$S0, decreasing = TRUE)                 # susceptibility ranking = the headline output
  cat("   S0 ranking:  ", paste(sprintf("%s=%.2f", sl$seasons[ord], p$S0[ord]), collapse = "  "), "\n")

  for (s in seq_along(sl$seasons)){
    y  = sl$ylist[[s]]; wk = sl$season_week[[s]]; mu = fit$mu[[s]]
    ok = is.finite(y) & is.finite(mu)
    corr = if (sum(ok) > 2) cor(y[ok], mu[ok]) else NA_real_
    cat(sprintf("   %s  S0=%.2f  c=%.0f  | cor=%.3f\n", sl$seasons[s], p$S0[s], p$c[s], corr))

    plot(wk, y, pch = 19, col = "grey30", xlab = "season week", ylab = "flu ILI+",
         main = sprintf("%s %s\nS0=%.2f  cor=%.2f", cc, sl$seasons[s], p$S0[s], corr))
    lines(wk, mu, col = "red", lwd = 2)                 # deterministic SIR fit
    abline(h = p$b, col = "grey60", lty = 3)            # shared off-season baseline
  }
  if (length(sl$seasons) < n_panels)                   # keep the grid aligned across rows
    for (k in seq_len(n_panels - length(sl$seasons))) plot.new()
}
legend("topright", c("observed", "SIR fit", "baseline"), bty = "n",
       col = c("grey30", "red", "grey60"), pch = c(19, NA, NA), lty = c(NA, 1, 3))
dev.off()
cat("\nfigure written to output/kalman_sir_multiseason_fit.png\n")
