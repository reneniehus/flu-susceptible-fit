# method_descriptive.R
#
# METHOD: descriptive curve-feature method (registry name "descriptive").
# Requires sir_core.R (the .curve_* feature helpers) to be sourced first.
#
# A deliberately NON-mechanistic, parsimonious method: it does not fit an SIR at all. It just
# SMOOTHS each observed ILI+ curve (centered moving average, NA-aware); the descriptive features --
# area under the curve, peak height, onset week, and STEEPNESS (the rate of rise) -- are then read off
# the smoothed curve uniformly for every method by summarise_method_fit(). This method's only job is
# to provide the smoothed curve as fit$mu; it reports no mechanistic parameters (S0/c/b/phi/qI are NA).
#
# Note on interpretation: steepness is an OBSERVED-SHAPE feature, not a clean susceptibility. The
# observed national rise reflects not only how fast the virus spreads but also how a country's many
# local epidemics overlay (staggered starts, travel) and how reporting is delayed in space and time.
# So steepness -- like AUC and peak height -- is compared across seasons WITHIN a country, not across
# countries. Mapping steepness onto an SIR susceptibility is a separate mechanistic-vs-phenomenological
# question and is deliberately out of scope here.
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

# ---- |-fit the descriptive method for one country (just smooth every season) ----
# Same call shape as the other methods (extra arguments are ignored via ...), so the registry can
# dispatch it uniformly. Only smooth_window is used; the features are read off fit$mu by the summary.
fit_descriptive = function(ylist, R0 = 1.5, infectious_period_days = 3, smooth_window = 4, ...){
  K = length(ylist)
  mu = lapply(ylist, function(y) pmax(.smooth_curve(y, smooth_window), 0))   # ILI+ cannot be negative
  list(method = "descriptive", R0 = R0, gamma = 7 / infectious_period_days, seed_i0 = NA_real_,
       params = list(S0 = rep(NA_real_, K), c = rep(NA_real_, K), b = NA_real_, phi = NA_real_,
                     qI = NA_real_, smooth_window = smooth_window),
       mu = mu, ylist = ylist, convergence = 0L, negloglik = NA_real_)
}
