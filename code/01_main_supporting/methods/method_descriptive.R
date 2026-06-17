# method_descriptive.R
#
# METHOD: descriptive curve-feature method (registry name "descriptive").
# Requires sir_core.R (the .curve_* feature helpers) to be sourced first.
#
# A deliberately NON-mechanistic, parsimonious method: it does not fit an SIR at all. It just
# SMOOTHS each observed ILI+ curve (centered moving average, NA-aware) and reads descriptive features
# off the smooth:
#   area under the curve, peak height, onset week, and a STEEPNESS (the exponential growth rate over
#   the rise). Steepness maps to susceptibility via the SIR rise-rate relation r = gamma*(R0*S0 - 1),
#   so the method also reports an implied S0 = (r/gamma + 1)/R0 comparable to the mechanistic methods.
# The features themselves (auc, peak_height, onset_week, steepness) are computed uniformly for every
# method by summarise_method_fit(); this method's job is to provide the smoothed curve as fit$mu.
#
# Key smoothing choice: a centered moving average (settings: params$susc_smooth_window, default 4)
# was preferred over loess -- on these regular, single-wave weekly curves it gives essentially the
# same features, is simpler/transparent, and the even "2x4" window preserves the peak while staying
# smooth. NA-aware, so it bridges short gaps in the series.
#
# Returns the common method-fit contract (mu = the smoothed curves; params$S0 = implied S0; the
# data-scale SIR parameters c/b/phi/qI are NA -- this method does not have them).

# ---- |-centered, NA-aware moving-average smooth of one season's weekly ILI+ curve ----
# Even windows use half-weights at the two ends (a symmetric "2xk" average, so there is no timing
# shift). The mean is taken over whatever is OBSERVED in the window, so short NA gaps are bridged;
# a week whose whole window is missing stays NA and is simply skipped by the feature helpers. (When
# more countries bring longer internal gaps, this is where a fill / interpolation step would go.)
.smooth_curve = function(y, window){
  n = length(y)
  if (window %% 2 == 1){ offs = -((window-1)/2):((window-1)/2); wts = rep(1, window) }
  else                 { offs = -(window/2):(window/2);          wts = c(0.5, rep(1, window-1), 0.5) }
  vapply(seq_len(n), function(i){
    idx = i + offs; ok = idx >= 1 & idx <= n
    iv = idx[ok]; wv = wts[ok]; keep = !is.na(y[iv])
    if (any(keep)) sum(wv[keep] * y[iv][keep]) / sum(wv[keep]) else NA_real_
  }, numeric(1))
}

# ---- |-fit the descriptive method for one country (smooth every season; implied S0 from steepness) ----
# Same call shape as the other methods (extra arguments are ignored via ...), so the registry can
# dispatch it uniformly. Only R0, the infectious period (for the S0 mapping) and smooth_window are used.
fit_descriptive = function(ylist, R0 = 1.5, infectious_period_days = 3, smooth_window = 4, ...){
  gamma = 7 / infectious_period_days
  K = length(ylist)
  mu = vector("list", K); S0 = numeric(K)
  for (s in seq_len(K)){
    y  = ylist[[s]]; wk = seq_along(y)                   # contiguous weekly grid (week index from season start)
    m  = pmax(.smooth_curve(y, smooth_window), 0); mu[[s]] = m   # ILI+ cannot be negative
    b  = .curve_baseline(m)
    r  = .curve_steepness(m, wk, b)                      # empirical exponential growth rate (per week)
    S0[s] = if (is.na(r)) NA_real_ else (r / gamma + 1) / R0   # S0 implied by the rise rate
  }
  list(method = "descriptive", R0 = R0, gamma = gamma, seed_i0 = NA_real_,
       params = list(S0 = S0, c = rep(NA_real_, K), b = NA_real_, phi = NA_real_, qI = NA_real_,
                     smooth_window = smooth_window),
       mu = mu, ylist = ylist, convergence = 0L, negloglik = NA_real_)
}
