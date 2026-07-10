# observe.R
#
# Stage 5: the OBSERVATION MODEL -- the deliberate wall between the latent truth and what an analyst
# actually sees. Nothing here changes the epidemic; it only degrades our view of it, in the specific
# ways real surveillance does:
#   - infections are never seen; we see ONSETS (infection + incubation), and only a time-varying
#     FRACTION of them (case ascertainment rho_case,t) as reported cases;
#   - each reported case appears in the record only after a REPORTING DELAY, so recent onset days are
#     right-truncated -- this is the reporting triangle the nowcast has to invert;
#   - deaths are a small fraction (IFR) of infections, land after an onset->death delay, and are
#     themselves right-censored -- the reason naive CFR is biased early;
#   - importations are detected only in proportion to a country's SURVEILLANCE QUALITY;
#   - only a SEQUENCED subsample of cases is typed, giving the variant-frequency series;
#   - flight volumes are seen through measurement noise (already produced in simulate_flights()).
#
# The output is two parallel worlds: `truth` (everything, exactly) and `observed` (the degraded,
# analyst-facing data). The whole point of the playground is that every analysis tool consumes only
# `observed` and is then scored against `truth`.
#
# Requires: R/utils.R sourced first.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Observe one location ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-degrade one location's latent infections into observed cases, deaths and typing ----
# infections: n x K matrix (columns = strains). Returns the location's observed & latent series.
.observe_location <- function(infections, incub, o2r, inf2death,
                              rho_case_vec, death_rho, ifr, seq_prob,
                              inf2adm = NULL, admission_rate = NULL) {
  n <- nrow(infections); K <- ncol(infections)
  infections_total <- rowSums(infections)

  # ---- onsets: infections shifted by the incubation period (expected onsets by onset day) ----
  onsets_total   <- shift_by_delay(infections_total, incub)
  onsets_strain  <- vapply(seq_len(K), function(k) shift_by_delay(infections[, k], incub), numeric(n))
  if (is.null(dim(onsets_strain))) onsets_strain <- matrix(onsets_strain, ncol = K)

  # ---- eventual reported cases by onset day: a Poisson draw thinned by time-varying ascertainment ----
  # This is the number that will EVENTUALLY be reported for each onset day (the nowcast target, kept
  # in `truth`). What an analyst sees now is the right-truncated version below.
  cases_by_onset_eventual <- stats::rpois(n, onsets_total * rho_case_vec)

  # ---- reporting triangle: scatter each onset day's cases across report days by the delay PMF ----
  # triangle[o, d] = cases with onset on day (o-1) that are reported at delay (d-1). The eventual
  # report day is (o-1)+(d-1); cells whose report day exceeds the horizon are simply never seen.
  D <- length(o2r) - 1
  triangle <- matrix(0L, nrow = n, ncol = D + 1)
  for (o in seq_len(n)) {
    if (cases_by_onset_eventual[o] > 0)
      triangle[o, ] <- as.integer(stats::rmultinom(1, cases_by_onset_eventual[o], o2r))
  }
  # cases by report day = anti-diagonal sums of the triangle, within the horizon (ALWAYS complete to
  # as_of -- this is the by-report-date epidemic curve). cases by onset day OBSERVED = row sums within
  # the horizon (RIGHT-TRUNCATED -- recent onsets look artificially low until their reports arrive).
  cases_by_report    <- numeric(n)
  cases_by_onset_obs <- numeric(n)
  for (o in seq_len(n)) for (d in 0:D) {
    if (triangle[o, d + 1] == 0) next
    r <- o + d                                        # report day index (1-based)
    if (r <= n) {
      cases_by_report[r]    <- cases_by_report[r]    + triangle[o, d + 1]
      cases_by_onset_obs[o] <- cases_by_onset_obs[o] + triangle[o, d + 1]
    }
  }

  # ---- variant typing: a sequenced subsample of cases, labelled by strain ----
  # Only a fraction seq_prob of reported cases are sequenced; among those the variant share follows
  # the onset-day strain composition. Gives the realistic, small-sample frequency series.
  seq_total <- stats::rbinom(n, cases_by_onset_eventual, seq_prob)
  if (K == 2L) {
    p_variant <- ifelse(onsets_total > 0, onsets_strain[, 2] / onsets_total, 0)
    seq_variant <- stats::rbinom(n, seq_total, p_variant)
  } else seq_variant <- integer(n)

  # ---- deaths: a fraction (IFR) of infections, dated by the infection->death delay, then detected ----
  # [REFLECTION] deaths are thinned from INFECTIONS independently of case detection, so a death need not
  # be a confirmed case. Consequence: the confirmed CFR (detected deaths / detected cases) equals
  # IFR * death_rho / case_rho -- it overstates the IFR whenever ascertainment < death detection, which
  # is exactly the early-pandemic bias the CFR tool then exposes. A line-list mode (deaths subset of
  # cases) is the alternative and would change that interpretation. See documentation/decisions.md.
  expected_deaths <- ifr * shift_by_delay(infections_total, inf2death)
  deaths_by_date  <- stats::rpois(n, expected_deaths * death_rho)      # deaths registered on each day
  # deaths also attributed to the fatal case's ONSET day (for onset-based CFR): IFR of onsets, detected
  deaths_by_onset <- stats::rpois(n, ifr * onsets_total * death_rho)

  # ---- hospital admissions: a fraction of infections, dated by the infection->admission delay ----
  # Treated as well-ascertained (hospital census), so no detection thinning -- the Phase 2 forecast
  # target. Absent config admissions, this stays zero (older configs remain valid).
  admissions_by_date <- if (!is.null(inf2adm) && !is.null(admission_rate))
    stats::rpois(n, admission_rate * shift_by_delay(infections_total, inf2adm)) else integer(n)

  list(
    infections_total = infections_total,
    onsets_total     = onsets_total,
    onsets_strain    = onsets_strain,
    admissions_by_date = admissions_by_date,
    cases_by_onset_eventual = cases_by_onset_eventual,   # eventual detected (nowcast target -> truth)
    cases_by_onset   = cases_by_onset_obs,               # right-truncated (analyst-facing -> observed)
    cases_by_report  = cases_by_report,
    triangle         = triangle,
    seq_total        = seq_total,
    seq_variant      = seq_variant,
    deaths_by_date   = deaths_by_date,
    deaths_by_onset  = deaths_by_onset
  )
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Observe the whole run ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-run the observation model over the source and every country; assemble truth + observed ----
# Returns list(truth = <latent everything>, observed = <degraded, analyst-facing>). Tidy long data
# frames throughout, so REAL surveillance data in the same schema can be dropped in unchanged
# (see documentation/real_data.md).
observe <- function(cfg, par, source, flights, imports, local) {
  n     <- cfg$n_days
  cc    <- cfg$countries
  day   <- par$grid$day
  date  <- par$grid$date
  K     <- par$n_strains

  incub     <- par$pmf$incubation
  o2r       <- par$pmf$onset_to_report
  o2d       <- par$pmf$onset_to_death
  inf2death <- .convolve_pmf(incub, o2d)                            # infection -> death delay
  inf2adm   <- if (!is.null(par$pmf$onset_to_admission)) .convolve_pmf(incub, par$pmf$onset_to_admission) else NULL
  adm_rate  <- cfg$admission_rate                                   # NULL in older configs -> admissions stay 0
  rho_case_vec <- step_at(cfg$ascertainment$case_rho, day)         # length n
  death_rho <- cfg$ascertainment$death_rho
  ifr       <- par$ifr_eff
  seq_prob  <- cfg$variant$sequencing_prob %||% 0.1

  # ---- observe the source X (its own case curve feeds the Phase-0 growth estimate) ----
  src_obs <- .observe_location(source$infections, incub, o2r, inf2death,
                               rho_case_vec, death_rho, ifr, seq_prob, inf2adm, adm_rate)

  # ---- observe every country ----
  loc_obs <- lapply(seq_along(cc), function(j)
    .observe_location(matrix(local$infections[, j, ], nrow = n, ncol = K),
                      incub, o2r, inf2death, rho_case_vec, death_rho, ifr, seq_prob, inf2adm, adm_rate))
  names(loc_obs) <- cc

  # ---- detected imports: true imports thinned by per-country surveillance quality ----
  # The importation SIGNAL: a well-surveilled country reveals near-true imports; a poor one misses
  # most. This is what the De Salazar regression exploits.
  detected_imports <- matrix(0L, nrow = n, ncol = length(cc), dimnames = list(NULL, cc))
  for (j in seq_along(cc))
    detected_imports[, j] <- stats::rbinom(n, imports$imports_total[, j], cfg$surveillance_quality[cc[j]])

  # ---- assemble tidy long data frames --------------------------------------------------------
  truth <- .assemble_truth(cfg, par, source, flights, imports, local, src_obs, loc_obs, ifr)
  observed <- .assemble_observed(cfg, par, flights, detected_imports, src_obs, loc_obs)

  list(truth = truth, observed = observed)
}

# ---- |-convolve two PMFs (both indexed from delay 0) into their sum's PMF ----
.convolve_pmf <- function(a, b) {
  out <- stats::convolve(a, rev(b), type = "open")                 # full linear convolution
  out[out < 0] <- 0
  out / sum(out)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Assemble the two worlds into tidy tables ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-latent TRUTH: infections, R_t, rho_t, IFR/CFR, imports, flights, variant frequency ----
.assemble_truth <- function(cfg, par, source, flights, imports, local, src_obs, loc_obs, ifr) {
  n <- cfg$n_days; cc <- cfg$countries; day <- par$grid$day; date <- par$grid$date

  # true reproduction number by day: source + each country (the step schedules, evaluated). We keep
  # BOTH the nominal per-strain (wild-type) schedule `Rt` and the REALIZED instantaneous reproduction
  # number `Rt_effective` -- the strain-mix- and depletion-weighted R the epidemic actually runs at,
  # which is what a Cori/EpiEstim estimate off the case curve actually targets.
  gi <- par$pmf$gi
  rt_list_src <- .strain_rt_list(cfg$rt_source, cfg, par)
  reff_src <- .realized_rt(source$infections, source$susceptible, cfg$source$population, rt_list_src, gi)
  rt_source <- data.frame(location = cfg$source$code, day = day, date = date,
                          Rt = step_at(cfg$rt_source, day), Rt_effective = reff_src,
                          stringsAsFactors = FALSE)
  rt_country <- do.call(rbind, lapply(seq_along(cc), function(j) {
    x <- cc[j]
    rl <- .strain_rt_list(cfg$rt_country[[x]], cfg, par)
    reff <- .realized_rt(matrix(local$infections[, j, ], nrow = n, ncol = par$n_strains),
                         local$susceptible[, j], par$country_pop[[x]], rl, gi)
    data.frame(location = x, day = day, date = date, Rt = step_at(cfg$rt_country[[x]], day),
               Rt_effective = reff, stringsAsFactors = FALSE)
  }))

  # true infections (long): source + countries, total across strains
  inf_source <- data.frame(location = cfg$source$code, day = day, date = date,
                           infections = source$infections_total, stringsAsFactors = FALSE)
  inf_country <- do.call(rbind, lapply(seq_along(cc), function(j)
    data.frame(location = cc[j], day = day, date = date,
               infections = local$infections_total[, j], stringsAsFactors = FALSE)))

  # true variant frequency (source + countries), among new infections
  vf_source <- data.frame(location = cfg$source$code, day = day, date = date,
                          variant_freq = source$variant_freq, stringsAsFactors = FALSE)
  vf_country <- do.call(rbind, lapply(seq_along(cc), function(j) {
    li <- local$infections; tot <- li[, j, 1] + (if (par$n_strains == 2L) li[, j, 2] else 0)
    vf <- if (par$n_strains == 2L) ifelse(tot > 0, li[, j, 2] / tot, NA_real_) else rep(0, n)
    data.frame(location = cc[j], day = day, date = date, variant_freq = vf, stringsAsFactors = FALSE)
  }))

  # true imports per country per day
  imports_true <- do.call(rbind, lapply(seq_along(cc), function(j)
    data.frame(country = cc[j], day = day, date = date, imports = imports$imports_total[, j],
               stringsAsFactors = FALSE)))

  # true flight volumes (long)
  flights_true <- do.call(rbind, lapply(seq_along(cc), function(j)
    data.frame(country = cc[j], day = day, date = date, volume = flights$true[, j],
               stringsAsFactors = FALSE)))

  # true onsets and the eventual-detected onset curve (the nowcast target), source + countries
  all_locs <- c(cfg$source$code, cc)
  all_obs  <- c(list(src_obs), loc_obs)
  onsets_true <- do.call(rbind, Map(function(loc, o)
    data.frame(location = loc, day = day, date = date, onsets = o$onsets_total,
               stringsAsFactors = FALSE), all_locs, all_obs))
  cases_onset_eventual <- do.call(rbind, Map(function(loc, o)
    data.frame(location = loc, day = day, date = date, cases = o$cases_by_onset_eventual,
               stringsAsFactors = FALSE), all_locs, all_obs))
  rownames(onsets_true) <- rownames(cases_onset_eventual) <- NULL

  list(
    infections       = rbind(inf_source, inf_country),
    onsets           = onsets_true,
    cases_by_onset   = cases_onset_eventual,             # eventual detected cases by onset (nowcast target)
    Rt               = rbind(rt_source, rt_country),
    case_ascertainment = data.frame(day = day, date = date,
                                    rho_case = step_at(cfg$ascertainment$case_rho, day),
                                    stringsAsFactors = FALSE),
    ifr              = ifr,
    death_ascertainment = cfg$ascertainment$death_rho,
    imports          = imports_true,
    flight_volumes   = flights_true,
    variant_freq     = rbind(vf_source, vf_country),
    capacity         = data.frame(country = cc,
                                  capacity = (cfg$capacity_per_million %||% NA_real_) * par$country_pop[cc] / 1e6,
                                  stringsAsFactors = FALSE),
    prevalence       = data.frame(location = cfg$source$code, day = day, date = date,
                                  prevalence = source$prevalence, prevalence_frac = source$prevalence_frac,
                                  stringsAsFactors = FALSE),
    source_size      = data.frame(day = day, date = date,
                                  infections = source$infections_total,
                                  cumulative = cumsum(source$infections_total),
                                  stringsAsFactors = FALSE)
  )
}

# ---- |-OBSERVED, analyst-facing data: cases (triangle + aggregations), deaths, imports, typing, flights ----
.assemble_observed <- function(cfg, par, flights, detected_imports, src_obs, loc_obs) {
  n <- cfg$n_days; cc <- cfg$countries; day <- par$grid$day; date <- par$grid$date
  src <- cfg$source$code

  onset_df <- function(loc, o) data.frame(location = loc, day = day, date = date,
                                          cases = o$cases_by_onset, stringsAsFactors = FALSE)
  report_df <- function(loc, o) data.frame(location = loc, day = day, date = date,
                                           cases = o$cases_by_report, stringsAsFactors = FALSE)
  deaths_df <- function(loc, o) data.frame(location = loc, day = day, date = date,
                                           deaths_by_date = o$deaths_by_date,
                                           deaths_by_onset = o$deaths_by_onset, stringsAsFactors = FALSE)
  adm_df    <- function(loc, o) data.frame(location = loc, day = day, date = date,
                                           admissions = o$admissions_by_date, stringsAsFactors = FALSE)
  variant_df <- function(loc, o) data.frame(location = loc, day = day, date = date,
                                            sequenced = o$seq_total, variant = o$seq_variant,
                                            stringsAsFactors = FALSE)

  all_locs   <- c(src, cc)
  all_obs    <- c(list(src_obs), loc_obs); names(all_obs) <- all_locs

  cases_by_onset  <- do.call(rbind, Map(onset_df,  all_locs, all_obs))
  cases_by_report <- do.call(rbind, Map(report_df, all_locs, all_obs))
  deaths          <- do.call(rbind, Map(deaths_df, all_locs, all_obs))
  admissions      <- do.call(rbind, Map(adm_df,    all_locs, all_obs))
  variant_cases   <- do.call(rbind, Map(variant_df, all_locs, all_obs))
  rownames(cases_by_onset) <- rownames(cases_by_report) <- rownames(deaths) <-
    rownames(admissions) <- rownames(variant_cases) <- NULL

  # reporting triangles: one long table (location, onset_day, report_day, cases) -- the nowcast input
  triangle <- do.call(rbind, Map(function(loc, o) .triangle_long(loc, o$triangle, n), all_locs, all_obs))

  # detected imports per country + the surveillance-quality covariate
  imports_obs <- do.call(rbind, lapply(seq_along(cc), function(j)
    data.frame(country = cc[j], day = day, date = date,
               detected_imports = detected_imports[, j],
               surveillance_quality = unname(cfg$surveillance_quality[cc[j]]),
               stringsAsFactors = FALSE)))

  # observed (noisy) flight volumes
  flights_obs <- do.call(rbind, lapply(seq_along(cc), function(j)
    data.frame(country = cc[j], day = day, date = date, volume = flights$observed[, j],
               stringsAsFactors = FALSE)))

  list(
    cases_by_onset    = cases_by_onset,
    cases_by_report   = cases_by_report,
    reporting_triangle = triangle,
    deaths            = deaths,
    admissions        = admissions,
    detected_imports  = imports_obs,
    variant_cases     = variant_cases,
    flight_volumes    = flights_obs,
    as_of             = max(day)                                   # the data cutoff (see everything to the horizon)
  )
}

# ---- |-turn a location's onset x delay triangle into a long (onset_day, report_day, cases) table ----
# Only nonzero cells within the horizon are emitted (sparse) -- the exact right-truncated record an
# analyst would hold at the data cutoff.
.triangle_long <- function(loc, tri, n) {
  D <- ncol(tri) - 1
  idx <- which(tri > 0, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(location = character(0), onset_day = integer(0),
                                    report_day = integer(0), cases = integer(0)))
  onset_day  <- idx[, 1] - 1
  report_day <- onset_day + (idx[, 2] - 1)
  keep <- report_day <= (n - 1)                                    # within the observation horizon
  data.frame(location = loc, onset_day = onset_day[keep], report_day = report_day[keep],
             cases = tri[idx][keep], stringsAsFactors = FALSE)
}

# ---- |-per-strain R schedules for a location (wild type + variant x (1+fitness)) ----
# Rebuilds the same strain R list the renewal engine used, so the realized-Rt truth is exact.
.strain_rt_list <- function(base_rt, cfg, par) {
  rl <- list(base_rt)
  if (par$n_strains == 2L) { vr <- base_rt; vr$value <- vr$value * (1 + cfg$variant$fitness); rl[[2]] <- vr }
  rl
}

# ---- |-realized instantaneous reproduction number: expected new infections / total infectiousness ----
# R_eff(t) = [ sum_k R_k(t) (S_{t-1}/N) Lambda_{k,t} ] / [ sum_k Lambda_{k,t} ], with
# Lambda_{k,t} = sum_{s>=1} g(s) I_{k,t-s}. This is the instantaneous R the epidemic actually runs at,
# accounting for the strain mix and susceptible depletion -- the quantity a Cori estimate targets.
.realized_rt <- function(infections, susceptible, N, rt_list, gi_pmf) {
  n <- nrow(infections); K <- ncol(infections); D <- length(gi_pmf) - 1
  reff <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    smax <- min(D, t - 1); if (smax < 1) next
    Sfrac <- if (t == 1) 1 else susceptible[t - 1] / N
    num <- 0; den <- 0
    for (k in seq_len(K)) {
      Lk  <- sum(gi_pmf[2:(smax + 1)] * infections[(t - 1):(t - smax), k])
      num <- num + step_at(rt_list[[k]], t - 1) * Sfrac * Lk
      den <- den + Lk
    }
    if (den > 0) reff[t] <- num / den
  }
  reff
}
