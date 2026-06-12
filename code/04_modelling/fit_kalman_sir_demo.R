# fit_kalman_sir_demo.R
#
# Demo: fit the EKF-SIR model (code/01_main_supporting/model_kalman_sir.R) to flu ILI+ for two
# countries and save a fit figure. This is the R analogue of stan/SIR_multiseason_age_vax_2.stan,
# using an Extended Kalman Filter for the likelihood and optim() (MAP) for fitting.
#
# Run from the repo root:  Rscript code/04_modelling/fit_kalman_sir_demo.R

source("code/01_main_supporting/setup.R")
source("code/02_settings/settings_version0.R"); params <- settings()
for (f in c("flu_functions", "validate", "load_data", "gen_model_input", "model_kalman_sir"))
  source(paste0("code/01_main_supporting/", f, ".R"))

data <- load_data_epi(load_data(params, regenerate = FALSE, new_from_online = FALSE),
                      params, regenerate = TRUE, new_from_online = FALSE)
models_in <- gen_model_input(params, data)

# two contrasting countries (clean, full single-peak seasons; very different magnitudes)
targets <- list(c("DK", "2023/2024"), c("GR", "2023/2024"))

dir.create("output", showWarnings = FALSE)
png("output/kalman_sir_fit.png", width = 1100, height = 520); par(mfrow = c(1, 2))
for (tg in targets) {
  cc <- tg[1]; ss <- tg[2]
  s   <- kalman_sir_series(models_in, cc, ss)
  fit <- fit_kalman_sir(s$value, infectious_period_days = 3, n_sub = 7, n_starts = 8)
  ok  <- is.finite(s$value) & is.finite(fit$mu_pred)
  corr <- cor(s$value[ok], fit$mu_pred[ok])
  cat(sprintf("%s %s: R0=%.2f S0=%.2f I0=%.1e c=%.0f b=%.1f phi=%.1f | cor=%.3f\n",
              cc, ss, fit$params$R0, fit$params$S0, fit$params$I0,
              fit$params$c, fit$params$b, fit$params$phi, corr))
  plot(s$season_week, s$value, pch = 19, col = "grey30", xlab = "season week", ylab = "flu ILI+",
       main = sprintf("%s %s  (R0=%.2f, cor=%.2f)", cc, ss, fit$params$R0, corr))
  lines(s$season_week, fit$mu_pred, col = "red", lwd = 2)                 # one-step-ahead filtered mean
  lines(s$season_week, kalman_sir_trajectory(fit), col = "blue", lwd = 1, lty = 2) # deterministic SIR
  abline(h = fit$params$b, col = "grey60", lty = 3)                      # baseline
  legend("topright", c("observed", "EKF fit", "SIR mean"), bty = "n",
         col = c("grey30", "red", "blue"), pch = c(19, NA, NA), lty = c(NA, 1, 2))
}
dev.off()
cat("figure written to output/kalman_sir_fit.png\n")
