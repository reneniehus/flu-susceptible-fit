# methods_registry.R
#
# The registry that turns the individual fitting methods (code/01_main_supporting/methods/) into
# one swappable framework. Each method is a file that defines a fit_*() returning the common
# method-fit contract:
#   list(method, R0, gamma, seed_i0, params = list(S0, c, b, phi, qI), mu = <per-season curves>,
#        ylist, convergence, negloglik)
# and shares the signature fit_*(ylist, R0, infectious_period_days, seed_i0, n_sub, n_starts, seed).
#
# To ADD A METHOD: write methods/method_<name>.R, source it, and add ONE line to sir_methods().
# Everything downstream (run, summarise, plot, correlate) is method-agnostic.
#
# Requires: sir_core.R and the method files to be sourced first.

# ---- |-method registry: name -> (human label, fit function) ----
# Built lazily so it only references fit functions that have been sourced.
sir_methods = function(){
  list(
    deterministic = list(label = "Deterministic SIR (no process noise)", fit = fit_sir_deterministic),
    ekf           = list(label = "EKF SIR (fitted process noise)",        fit = fit_sir_ekf)
  )
}

# ---- |-fit one method on one country's panel; attach the country/season metadata ----
# panel: output of load_flu_iliplus_slim(). params: settings() list (susc_* fields).
run_method = function(method_name, panel, params, n_sub = 7, n_starts = 4, seed = 1){
  spec = sir_methods()[[method_name]]
  if (is.null(spec)) stop(sprintf("run_method: unknown method '%s'", method_name))
  fit = spec$fit(panel$ylist, R0 = params$susc_R0,
                 infectious_period_days = params$susc_infectious_period_days,
                 seed_i0 = params$susc_seed_i0, n_sub = n_sub, n_starts = n_starts, seed = seed)
  fit$method_name  = method_name
  fit$method_label = spec$label
  fit$country      = panel$country
  fit$seasons      = panel$seasons
  fit$season_week  = panel$season_week
  fit
}

# ---- |-onset week: first week the fitted curve rises to a fraction of its above-baseline peak ----
# Simple, timing-of-curve threshold read off the red curve (mu), robust to the baseline level.
.onset_week = function(mu, season_week, b, onset_frac = 0.1){
  above = mu - b
  pk = max(above, na.rm = TRUE)
  if (!is.finite(pk) || pk <= 0) return(NA_integer_)
  idx = which(above >= onset_frac * pk)[1]
  if (is.na(idx)) NA_integer_ else season_week[idx]
}

# ---- |-standard per-season summary statistics from a method fit (the cross-method schema) ----
# One row per (country, season): the susceptibility (S0, R_eff), reporting fraction c, the fitted
# curve's timing (peak + onset week), the process noise, and fit diagnostics. This is the single
# table the downstream cross-season / cross-country analysis consumes.
summarise_method_fit = function(fit, onset_frac = 0.1){
  pn = if (is.null(fit$params$qI)) NA_real_ else fit$params$qI[1]
  rows = lapply(seq_along(fit$seasons), function(s){
    mu = fit$mu[[s]]; wk = fit$season_week[[s]]; y = fit$ylist[[s]]
    ok = is.finite(y) & is.finite(mu)
    data.frame(
      method        = fit$method_name,
      country       = fit$country,
      season        = fit$seasons[s],
      S0            = fit$params$S0[s],
      R_eff         = fit$R0 * fit$params$S0[s],
      c             = fit$params$c[s],
      peak_week     = wk[which.max(mu)],
      onset_week    = .onset_week(mu, wk, fit$params$b, onset_frac),
      cor           = if (sum(ok) > 2) cor(y[ok], mu[ok]) else NA_real_,
      process_noise = pn,
      convergence   = fit$convergence,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- |-run every (country x method) and bind the per-season summaries into one tidy table ----
run_all_methods = function(countries, params, methods = names(sir_methods()),
                           n_starts = 4, verbose = TRUE){
  out = list()
  for (cc in countries){
    panel = load_flu_iliplus_slim(cc)
    for (m in methods){
      t = system.time(fit <- run_method(m, panel, params, n_starts = n_starts))[["elapsed"]]
      if (verbose) cat(sprintf("  %-4s %-14s %2d seasons  conv=%d  %4.0fs\n",
                               cc, m, length(panel$seasons), fit$convergence, t))
      out[[paste(cc, m, sep = "_")]] = summarise_method_fit(fit)
    }
  }
  do.call(rbind, out)
}
