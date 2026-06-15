# fit_kalman_sir_multiseason_demo.R
#
# Demo: the Stan-matched re-parameterisation of the EKF-SIR (code/01_main_supporting/
# model_kalman_sir.R). Unlike fit_kalman_sir_demo.R -- which fits ONE season and lets R0 float --
# here R0 is FIXED from the literature and the model fits the initial condition (S0, I0)
# SEPARATELY FOR EACH SEASON while SHARING the reporting fraction c and overdispersion phi across
# seasons. That is exactly the fitted-parameter set of stan/SIR_multiseason_age_vax_2.stan.
#
# Runs on the committed slim panel (data/slim_flu_iliplus.csv: DK & FR flu sentinel ILI+), so it
# needs only base R -- no data pipeline, no tidyverse.
#
# Run from the repo root:  Rscript code/04_modelling/fit_kalman_sir_multiseason_demo.R

source("code/01_main_supporting/model_kalman_sir.R")

# ---- |-settings ----
R0_flu    = 1.3                                        # fixed seasonal-influenza R0 (literature)
countries = c("DK", "FR")
n_panels  = 5                                          # max season panels per country row

# ---- |-fit each country jointly across its seasons, then plot every season ----
dir.create("output", showWarnings = FALSE)
png("output/kalman_sir_multiseason_fit.png", width = 1300, height = 560)
par(mfrow = c(length(countries), n_panels), mar = c(4, 4, 3, 1))

for (cc in countries){
  sl  = load_flu_iliplus_slim(cc)
  t   = system.time(
    fit <- fit_kalman_sir_multiseason(sl$ylist, R0 = R0_flu, infectious_period_days = 3,
                                      n_sub = 7, n_starts = 4))[["elapsed"]]
  p = fit$params
  cat(sprintf("\n%s  (R0 fixed = %.2f, %d seasons -> %d params)  shared: c=%.0f  phi=%.1f  | conv=%d  %.0fs\n",
              cc, fit$R0, length(sl$seasons), 2*length(sl$seasons) + 2, p$c, p$phi, fit$convergence, t))

  for (s in seq_along(sl$seasons)){
    y    = sl$ylist[[s]]; wk = sl$season_week[[s]]
    mu   = fit$fits[[s]]$mu_pred                        # one-step-ahead filtered mean (tracking)
    traj = kalman_sir_ms_trajectory(fit, s)            # deterministic SIR (susceptible recon.)
    ok   = is.finite(y) & is.finite(mu)
    cor_filt = if (sum(ok) > 2) cor(y[ok], mu[ok])   else NA_real_
    cor_det  = if (sum(ok) > 2) cor(y[ok], traj[ok]) else NA_real_
    cat(sprintf("  %s  S0=%.2f  I0=%.1e  | cor_filt=%.3f  cor_det=%.3f\n",
                sl$seasons[s], p$S0[s], p$I0[s], cor_filt, cor_det))

    plot(wk, y, pch = 19, col = "grey30", xlab = "season week", ylab = "flu ILI+",
         main = sprintf("%s %s\nS0=%.2f  cor=%.2f", cc, sl$seasons[s], p$S0[s], cor_det))
    lines(wk, mu,   col = "red",  lwd = 2)                                  # EKF filtered fit
    lines(wk, traj, col = "blue", lwd = 1, lty = 2)                         # deterministic SIR
  }
  if (length(sl$seasons) < n_panels)                   # keep the grid aligned across rows
    for (k in seq_len(n_panels - length(sl$seasons))) plot.new()
}
legend("topright", c("observed", "EKF fit", "SIR mean"), bty = "n",
       col = c("grey30", "red", "blue"), pch = c(19, NA, NA), lty = c(NA, 1, 2))
dev.off()
cat("\nfigure written to output/kalman_sir_multiseason_fit.png\n")
