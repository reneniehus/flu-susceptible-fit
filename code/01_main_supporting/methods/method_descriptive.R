# method_descriptive.R
#
# METHOD: descriptive curve-feature method (registry name "descriptive").
# Requires sir_core.R (the .curve_* feature helpers) to be sourced first.
#
# A deliberately NON-mechanistic, parsimonious method: it does not fit an SIR at all. It just
# SMOOTHS each observed ILI+ curve (loess, degree 2) and reads descriptive features off the smooth:
#   area under the curve, peak height, onset week, and a STEEPNESS (the exponential growth rate over
#   the rise). Steepness maps to susceptibility via the SIR rise-rate relation r = gamma*(R0*S0 - 1),
#   so the method also reports an implied S0 = (r/gamma + 1)/R0 comparable to the mechanistic methods.
# The features themselves (auc, peak_height, onset_week, steepness) are computed uniformly for every
# method by summarise_method_fit(); this method's job is to provide the smoothed curve as fit$mu.
#
# Key smoothing choice: loess span (settings: params$susc_smooth_span, default 0.3) -- chosen to
# preserve the peak while staying smooth, and robust (unlike a GCV spline) on these short series.
#
# Returns the common method-fit contract (mu = the smoothed curves; params$S0 = implied S0; the
# data-scale SIR parameters c/b/phi/qI are NA -- this method does not have them).

# ---- |-smooth one season's weekly ILI+ curve (loess, degree 2); NA outside the observed range ----
.smooth_curve = function(week, y, span){
  d = data.frame(week = week, y = y); dok = d[is.finite(d$y), ]
  fit = loess(y ~ week, data = dok, span = span, degree = 2)
  pmax(as.numeric(predict(fit, newdata = data.frame(week = week))), 0)   # ILI+ cannot be negative
}

# ---- |-fit the descriptive method for one country (smooth every season; implied S0 from steepness) ----
# Same call shape as the other methods (extra arguments are ignored via ...), so the registry can
# dispatch it uniformly. Only R0, the infectious period (for the S0 mapping) and smooth_span are used.
fit_descriptive = function(ylist, R0 = 1.5, infectious_period_days = 3, smooth_span = 0.3, ...){
  gamma = 7 / infectious_period_days
  K = length(ylist)
  mu = vector("list", K); S0 = numeric(K)
  for (s in seq_len(K)){
    y  = ylist[[s]]; wk = seq_along(y)                   # contiguous weekly grid (week index from season start)
    m  = .smooth_curve(wk, y, smooth_span); mu[[s]] = m
    b  = .curve_baseline(m)
    r  = .curve_steepness(m, wk, b)                      # empirical exponential growth rate (per week)
    S0[s] = if (is.na(r)) NA_real_ else (r / gamma + 1) / R0   # S0 implied by the rise rate
  }
  list(method = "descriptive", R0 = R0, gamma = gamma, seed_i0 = NA_real_,
       params = list(S0 = S0, c = rep(NA_real_, K), b = NA_real_, phi = NA_real_, qI = NA_real_,
                     smooth_span = smooth_span),
       mu = mu, ylist = ylist, convergence = 0L, negloglik = NA_real_)
}
