# phase0_importation_risk.R
#
# PHASE 0 -- before local introduction.
# QUESTION: "Which of our countries import next, and who is already missing cases?"
#
# THE LEAN METHOD (De Salazar et al. 2020). Expected imported infections into a country are
# proportional to its air-travel volume from the source times the source's prevalence -- so on any
# given window, detected imports should scale with flight volume. Fit that relationship on the
# WELL-SURVEILLED countries only (a Poisson regression of detected imports on log flight volume), then
# use it to PREDICT how many imports each country should have detected. Countries whose detected count
# falls below the lower prediction interval are detecting fewer than their connectivity implies -- they
# are the likely under-detectors, missing imported cases. Countries above/along the line are the ones
# to watch for the next introductions.
#
# WHY REGRESS ON THE WELL-SURVEILLED ONLY. They define what "detecting your fair share" looks like;
# extrapolating that line to everyone else exposes the gap. Using ALL countries would let the
# under-detectors drag the line down and hide themselves.
#
# References
#   De Salazar PM, Niehus R, Taylor A, Buckee CO, Lipsitch M. Identifying locations with possible
#     undetected imported SARS-CoV-2 cases by using importation predictions. Emerg Infect Dis. 2020.
#
# Requires: analysis_common.R.

# ---- |-regress detected imports on flight volume (well-surveilled anchor) and flag under-detectors ----
# input : as_analysis_input(sim). window : aggregation window (days). min_surveillance : the anchor
# threshold. level : prediction-interval level. Returns a per-country table with the expected imports,
# the prediction interval, and an under-detection flag.
importation_risk <- function(input, window, min_surveillance = 0.8, level = 0.95) {
  imp <- input$detected_imports; vol <- input$flight_volumes
  in_win <- function(d) d$day %in% window

  # aggregate detected imports and flight volume per country over the window
  agg_imp <- stats::aggregate(detected_imports ~ country, data = imp[in_win(imp), ], sum)
  agg_vol <- stats::aggregate(volume ~ country, data = vol[in_win(vol), ], sum)
  sq <- tapply(imp$surveillance_quality, imp$country, function(x) x[1])
  d  <- merge(agg_imp, agg_vol, by = "country")
  d$surveillance_quality <- as.numeric(sq[d$country])
  d$log_volume <- log(d$volume)

  # fit the imports-vs-volume line on the well-surveilled anchor countries only
  ref <- d[d$surveillance_quality >= min_surveillance, ]
  if (nrow(ref) < 3) stop("importation_risk: fewer than 3 anchor countries; lower min_surveillance")
  fit <- stats::glm(detected_imports ~ log_volume, family = stats::poisson(), data = ref)

  # predict expected imports for EVERY country, with a prediction interval (estimation + Poisson noise)
  pr <- stats::predict(fit, newdata = d, type = "link", se.fit = TRUE)
  z  <- stats::qnorm(1 - (1 - level) / 2)
  mu       <- exp(pr$fit)
  mu_lower <- exp(pr$fit - z * pr$se.fit)                 # parameter uncertainty on the mean
  mu_upper <- exp(pr$fit + z * pr$se.fit)
  # prediction interval: fold in Poisson sampling variability around the (uncertain) mean
  d$expected  <- mu
  d$pi_lower  <- stats::qpois((1 - level) / 2, mu_lower)
  d$pi_upper  <- stats::qpois(1 - (1 - level) / 2, mu_upper)
  d$under_detecting <- d$detected_imports < d$pi_lower
  d$residual <- (d$detected_imports - d$expected) / sqrt(pmax(d$expected, 1))   # Pearson-ish residual

  d[order(d$residual), c("country", "volume", "surveillance_quality", "detected_imports",
                         "expected", "pi_lower", "pi_upper", "under_detecting", "residual")]
}

# ---- |-score the flagging against the truth (which countries really are poorly surveilled) ----
# The truth is surveillance_quality. A good flagger puts low-surveillance countries below the line:
# residuals should correlate with surveillance quality, and flagged countries should have lower
# quality than unflagged ones.
score_importation_risk <- function(ir) {
  cor_rq <- stats::cor(ir$residual, ir$surveillance_quality)
  flagged <- ir$under_detecting
  list(
    cor_residual_surveillance = cor_rq,                  # expect strongly positive
    n_flagged = sum(flagged),
    mean_surv_flagged   = if (any(flagged))  mean(ir$surveillance_quality[flagged])  else NA_real_,
    mean_surv_unflagged = if (any(!flagged)) mean(ir$surveillance_quality[!flagged]) else NA_real_
  )
}
