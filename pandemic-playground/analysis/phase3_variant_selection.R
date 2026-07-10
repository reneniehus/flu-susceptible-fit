# phase3_variant_selection.R
#
# PHASE 3 -- sustained transmission, later waves, a new variant.
# QUESTION: "Is a new variant taking over, and how fast?"
#
# THE LEAN METHOD. Under a constant transmissibility advantage, the LOG-ODDS of the variant frequency
# grow linearly in time -- the hallmark of selection. So fit a binomial GLM of the sequenced counts,
#     logit( P(variant) )  =  a + s * day
# and read off the SELECTION COEFFICIENT s (the per-day increase in log-odds) with its confidence
# interval. The fitted line gives the frequency trajectory and the day the variant passes 50%. For
# more than two co-circulating variants the same idea generalises to a multinomial GLM
# (nnet::multinom) -- provided here as an optional path.
#
# WHY LOG-ODDS. Two strains growing exponentially at rates r_variant and r_wild have a frequency ratio
# that grows at exactly r_variant - r_wild; that difference IS the selection coefficient s, and it is
# linear on the logit scale even while both strains rise and fall. This is why s is a cleaner, more
# portable summary of a variant's advantage than any single-strain growth rate.
#
# References
#   Chevallier et al. / Volz et al. -- selection coefficients from variant frequency dynamics.
#   nnet::multinom -- multinomial logistic regression for >2 variants.
#
# Requires: analysis_common.R.

# ---- |-selection coefficient from sequenced variant counts (binomial GLM on the log-odds) ----
# input : as_analysis_input(sim). location. window : day indices over which to fit (the takeover
# window). Returns s (log-odds/day) with CI, the fitted frequency line, and the 50% crossover day.
variant_selection <- function(input, location, window = NULL, level = 0.95) {
  vc <- loc_series(input$variant_cases, location)
  if (!is.null(window)) vc <- vc[vc$day %in% window, ]
  vc <- vc[vc$sequenced > 0, ]                           # only days with typing data
  if (nrow(vc) < 5) stop("variant_selection: too few sequenced days to fit")

  fit <- stats::glm(cbind(variant, sequenced - variant) ~ day, family = stats::binomial(), data = vc)
  co  <- stats::coef(fit); se <- sqrt(diag(stats::vcov(fit))); z <- stats::qnorm(1 - (1 - level) / 2)
  s   <- unname(co["day"]); s_se <- unname(se["day"]); intercept <- unname(co["(Intercept)"])
  crossover <- if (s != 0) -intercept / s else NA_real_  # day the fitted frequency reaches 50%

  freq_line <- data.frame(day = vc$day,
                          observed = vc$variant / vc$sequenced,
                          fitted = stats::plogis(intercept + s * vc$day))
  list(location = location, s = s, s_ci = c(s - z * s_se, s + z * s_se),
       crossover_day = crossover, n_days = nrow(vc), freq_line = freq_line)
}

# ---- |-score the estimated selection coefficient against the truth ----
# Two truths: (a) the REALIZED selection coefficient -- the logit-slope of the true variant frequency
# over the same window (what actually happened); (b) the design-implied advantage from the config R
# step and fitness. We score against (a); (b) is reported for context.
score_variant_selection <- function(sim, vs, window = NULL) {
  loc <- vs$location; vf <- truth_variant_freq(sim, loc)
  if (!is.null(window)) vf <- vf[vf$day %in% window, ]
  vf <- vf[is.finite(vf$variant_freq) & vf$variant_freq > 1e-4 & vf$variant_freq < 1 - 1e-4, ]
  realized_s <- if (nrow(vf) > 3) unname(stats::coef(stats::lm(stats::qlogis(variant_freq) ~ day, data = vf))[2]) else NA_real_

  list(location = loc, s_estimate = vs$s, s_realized = realized_s,
       error = vs$s - realized_s,
       in_ci = !is.na(realized_s) && realized_s >= vs$s_ci[1] && realized_s <= vs$s_ci[2],
       crossover_day = vs$crossover_day)
}
