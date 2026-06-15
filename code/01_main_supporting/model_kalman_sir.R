# model_kalman_sir.R
#
# Base-R SIR fits for flu ILI+, sharing one generative story with the Stan model
# (stan/SIR_multiseason_age_vax_2.stan), single-population form:
#   - latent SIR in population proportions:  dS/dt = -beta S I,  dI/dt = beta S I - gamma I
#   - weekly ILI+ is proportional to weekly new infections:  E[y_t] = c * (new infections in week t)
#     (c folds in the infection->ILI fraction "prop_ili", positivity scaling and the per-100k rate)
#   - observation noise is neg-binomial-like:  Var(y) = mu + mu^2 / phi   (matches the Stan obs model)
#
# This file holds two fits on that model:
#   1. Extended Kalman Filter (single season, fit_kalman_sir): a state-space TRACKING fit that
#      fits R0, S0, I0, c, b, phi and a process-noise sd q_I by maximising the EKF log-likelihood.
#      Good for nowcasting; the filter can track a wave even where the dynamics are mis-specified.
#   2. Deterministic per-season SUSCEPTIBILITY fit (fit_sir_susceptibility): R0, gamma and the seed
#      I0 are FIXED (from settings); each season's susceptibility S0 and reporting fraction c are
#      fitted, with a shared baseline b and overdispersion phi. No filter -- S0 is identified by the
#      wave's rise rate, r = gamma*(R0*S0 - 1). This is the path for reading per-season susceptibility.

# ---- |-one sub-step (Euler) of the SIR in proportions, tracking cumulative incidence C ----
.sir_substep = function(x, beta, gamma, dt){
  S = x[1]; I = x[2]
  newinf = beta * S * I            # infection rate (per capita)
  dS = -newinf
  dI =  newinf - gamma * I
  c(S + dt*dS, I + dt*dI, x[3] + dt*newinf)
}

# ---- |-integrate one week (n_sub Euler steps); C enters at 0 so output C = weekly incidence ----
.sir_week = function(SI, beta, gamma, n_sub){
  x = c(SI[1], SI[2], 0)           # reset cumulative-incidence accumulator
  dt = 1 / n_sub
  for (k in seq_len(n_sub)) x = .sir_substep(x, beta, gamma, dt)
  x                                # (S_end, I_end, C_week)
}

# ---- |-Jacobian of the one-week map wrt (S, I) via finite differences (cheap, 2x perturbations) ----
.sir_week_jac = function(SI, beta, gamma, n_sub, eps=1e-6){
  f0 = .sir_week(SI, beta, gamma, n_sub)
  J  = matrix(0, 3, 2)             # d(S_end,I_end,C_week) / d(S0,I0)
  for (j in 1:2){
    dSI = SI; dSI[j] = dSI[j] + eps
    J[, j] = (.sir_week(dSI, beta, gamma, n_sub) - f0) / eps
  }
  list(f0 = f0, J = J)
}

# ---- |-unpack the (transformed) parameter vector into natural-scale params ----
.kf_unpack = function(par){
  list(
    R0  = exp(par[["log_R0"]]),                       # basic reproduction number (>0)
    S0  = plogis(par[["logit_S0"]]),                  # initial susceptible proportion (0,1)
    I0  = exp(par[["log_I0"]]),                       # initial infectious proportion (>0, small)
    c   = exp(par[["log_c"]]),                        # infection -> observed-ILI+ scale
    b   = exp(par[["log_b"]]),                        # off-season ILI+ baseline (sporadic detections)
    phi = exp(par[["log_phi"]]),                      # neg-binomial-like overdispersion
    qI  = exp(par[["log_qI"]])                        # process-noise sd on I
  )
}

# ---- |-weakly-informative log-priors on the scale-free parameters ----
# These regularise the fit (MAP estimation) and prevent the SIR identifiability degeneracies
# (e.g. R0 -> Inf with S0 -> 0, or phi -> Inf removing all observation noise). They mirror the
# priors that the Stan model declares but currently leaves commented out. Only the scale-free
# parameters are penalised; the data-scale parameters c and b are left free.
.kf_logprior = function(par){
  dnorm(par[["log_R0"]],   log(1.4),     0.5, log = TRUE) +   # R0 ~ roughly [0.7, 3]
  dnorm(par[["logit_S0"]], qlogis(0.55), 1.2, log = TRUE) +   # S0 not pinned to 0 or 1
  dnorm(par[["log_I0"]],   log(1e-5),    1.5, log = TRUE) +   # plausible epidemic seed
  dnorm(par[["log_phi"]],  log(15),      0.8, log = TRUE) +   # keep observation noise finite
  dnorm(par[["log_qI"]],   log(1e-4),    1.5, log = TRUE)
}

# ---- |-core EKF pass over one weekly series: returns log-lik + filtered trajectory ----
# The state-space TRACKING likelihood used by the single-season fit (ekf_sir_negloglik). All
# parameters arrive on the natural scale; beta = R0 * gamma.
# y: numeric vector (NA allowed). c, b: obs scale + baseline. qI: process-noise sd on I.
.ekf_filter = function(y, S0, I0, c, b, phi, beta, gamma, n_sub, qI){
  Tn = length(y)

  # state x = (S, I, C); start with C = 0
  x = c(S0, I0, 0)
  P = diag(c((0.3*S0)^2, (1.0*I0)^2, 0))              # initial state covariance
  Q = diag(c(0, qI^2, 0))                             # process noise (on I)

  ll = 0
  mu_pred = rep(NA_real_, Tn); I_filt = rep(NA_real_, Tn); S_filt = rep(NA_real_, Tn)

  for (t in seq_len(Tn)){
    # reset within-week cumulative incidence (deterministically known = 0)
    x[3] = 0; P[3, ] = 0; P[, 3] = 0

    # ---- predict one week ----
    wj   = .sir_week_jac(x[1:2], beta, gamma, n_sub)
    xpr  = wj$f0
    A    = cbind(wj$J, c(0, 0, 0))                     # 3x3: only S,I propagate (C reset each week)
    Ppr  = A %*% P %*% t(A) + Q
    xpr[1:2] = pmin(pmax(xpr[1:2], 1e-12), 1)          # keep S,I in valid range
    xpr[3]   = max(xpr[3], 1e-12)

    mu = c * xpr[3] + b                                # expected observed ILI+ (epidemic + baseline)
    mu_pred[t] = mu

    if (!is.na(y[t])){
      H   = matrix(c(0, 0, c), nrow = 1)
      Rt  = mu + mu^2 / phi                            # neg-binomial-like observation variance
      Sinn = max(as.numeric(H %*% Ppr %*% t(H)) + Rt, 1e-12)  # innovation var (guard tiny/neg)
      K   = (Ppr %*% t(H)) / Sinn
      innov = y[t] - mu
      x = as.numeric(xpr + K * innov)
      P = (diag(3) - K %*% H) %*% Ppr
      ll = ll + dnorm(y[t], mean = mu, sd = sqrt(Sinn), log = TRUE)
    } else {
      x = xpr; P = Ppr
    }
    x[1:2] = pmin(pmax(x[1:2], 1e-12), 1)              # guard the filtered state
    S_filt[t] = x[1]; I_filt[t] = x[2]
  }
  list(loglik = ll, mu_pred = mu_pred, S = S_filt, I = I_filt)
}

# ---- |-EKF negative log-(posterior) for one weekly ILI+ series ----
# y: numeric vector (NA allowed for missing weeks). gamma in per-week units. Returns -loglik
# (plus -logprior when use_prior=TRUE, i.e. the negative log-posterior used for MAP fitting).
ekf_sir_negloglik = function(par, y, gamma, n_sub = 7, use_prior = TRUE, return_fit = FALSE){
  p = .kf_unpack(par)
  beta = p$R0 * gamma
  f = .ekf_filter(y, p$S0, p$I0, p$c, p$b, p$phi, beta, gamma, n_sub, p$qI)

  lp = if (use_prior) .kf_logprior(par) else 0

  if (return_fit) return(list(negloglik = -(f$loglik + lp), loglik = f$loglik, mu_pred = f$mu_pred,
                              S = f$S, I = f$I, params = p))
  if (!is.finite(f$loglik + lp)) return(1e10)
  -(f$loglik + lp)
}

# ---- |-fit the EKF-SIR to one series via multi-start optim ----
# infectious_period_days sets gamma (per week). Returns best fit + the filtered trajectory.
fit_kalman_sir = function(y, infectious_period_days = 3, n_sub = 7, n_starts = 6, seed = 1){
  gamma = 7 / infectious_period_days                   # per-week recovery rate
  ypos  = y[is.finite(y) & y > 0]
  ymax  = max(ypos); base_guess = max(as.numeric(quantile(ypos, 0.1)), 1e-3)

  # a sensible starting point, then jittered restarts
  base = c(log_R0 = log(1.3), logit_S0 = qlogis(0.6), log_I0 = log(1e-5),
           log_c = log(ymax / 0.02), log_b = log(base_guess), log_phi = log(10), log_qI = log(1e-4))
  set.seed(seed)
  best = NULL
  for (s in seq_len(n_starts)){
    start = base + if (s == 1) 0 else rnorm(length(base), 0, c(0.4, 0.8, 1.0, 0.8, 0.6, 0.6, 1.0))
    fit = tryCatch(
      optim(start, ekf_sir_negloglik, y = y, gamma = gamma, n_sub = n_sub,
            method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-9)),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && (is.null(best) || fit$value < best$value)) best = fit
  }
  if (is.null(best)) stop("fit_kalman_sir: all optim starts failed")

  fitted = ekf_sir_negloglik(best$par, y = y, gamma = gamma, n_sub = n_sub, return_fit = TRUE)
  list(par = best$par, negloglik = best$value, convergence = best$convergence,
       gamma = gamma, infectious_period_days = infectious_period_days,
       params = fitted$params, mu_pred = fitted$mu_pred, S = fitted$S, I = fitted$I, y = y)
}

# ---- |-deterministic SIR trajectory from fitted params (smooth curve / projection) ----
# Runs the fitted SIR forward without filter updates, giving the model's own expected ILI+ curve
# (the analogue of the Stan model's generated-quantities mean). Set n_weeks > length(y) to project.
kalman_sir_trajectory = function(fit, n_weeks = length(fit$y), n_sub = 7){
  p = fit$params; beta = p$R0 * fit$gamma
  SI = c(p$S0, p$I0); mu = numeric(n_weeks)
  for (t in seq_len(n_weeks)){
    w = .sir_week(SI, beta, fit$gamma, n_sub)
    mu[t] = p$c * w[3] + p$b
    SI = w[1:2]
  }
  mu
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Per-season susceptibility fit: deterministic SIR, fixed R0 and seed ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Goal: read a per-season "sense of susceptibility" off the observed ILI+ wave. Susceptibility
# (pre-existing exposure, vaccine uptake/match, circulating strains) is summarised by the initial
# susceptible fraction S0, which the wave's RISE RATE identifies:  r = gamma*(R0*S0 - 1).  This is
# a DETERMINISTIC SIR (no Kalman filter): the per-season S0 has to be earned by the dynamics, not
# laundered through state updates. It is the same deterministic-SIR -> overdispersed-ILI generative
# model the Stan file fits by HMC (stan/SIR_multiseason_age_vax_2.stan), in single-population form.
#
# Fixed by assumption (from settings): R0 (literature), gamma (infectious period) and the seed I0
# (a small constant import -- "flu is always seeded from the southern hemisphere"), planted at week
# 1 of every season (the season start, ~August). Because I0 is fixed and shared across seasons, the
# ABSOLUTE S0 is conditional on R0/I0, but the RELATIVE S0 across seasons (which season is more
# susceptible) is robust -- that ranking is the quantity of interest.
#
# Fitted, per country (K seasons):
#   per season : S0[1..K] (susceptibility) and c[1..K] (reporting fraction = detection likelihood)
#   shared     : b (additive off-season baseline) and phi (observation overdispersion)

# ---- |-deterministic weekly new-infection incidence for one season (proportions) ----
# Seeds (S0, I0) at week 1 and integrates the SIR forward; returns the weekly new-infection
# proportion (the C accumulator of .sir_week), i.e. the quantity the reporting fraction c scales.
.sir_season_incidence = function(S0, I0, beta, gamma, n_weeks, n_sub){
  SI = c(S0, I0); inc = numeric(n_weeks)
  for (t in seq_len(n_weeks)){
    w = .sir_week(SI, beta, gamma, n_sub)
    inc[t] = w[3]; SI = w[1:2]
  }
  inc
}

# ---- |-unpack the susceptibility parameter vector into natural-scale params (K seasons) ----
.susc_unpack = function(par, K){
  list(
    S0  = plogis(par[1:K]),                              # per-season susceptibility (0,1)
    c   = exp(par[(K+1):(2*K)]),                          # per-season reporting fraction / obs scale
    b   = exp(par[[2*K + 1]]),                            # shared additive off-season baseline
    phi = exp(par[[2*K + 2]])                             # shared observation overdispersion
  )
}

# ---- |-weakly-informative log-priors (regularise S0 and phi; c, b left free) ----
# With R0 fixed the wave only grows when R0*S0 > 1, so S0 is centred in the plausible epidemic
# range; phi is kept finite. The data-scale reporting c and baseline b are unpenalised.
.susc_logprior = function(par, K){
  sum(dnorm(par[1:K], qlogis(0.75), 1.0, log = TRUE)) +  # logit_S0 per season
  dnorm(par[[2*K + 2]], log(15), 0.8, log = TRUE)        # shared log_phi
}

# ---- |-deterministic negative log-(posterior) over a country's seasons (fixed R0, I0) ----
# ylist: list of weekly ILI+ vectors, one per season (NA allowed), each seeded at week 1. The
# expected ILI+ is  mu = c[s] * (weekly new infections) + b,  with overdispersed observation noise
# Var = mu + mu^2/phi (neg-binomial-like, the Gaussian analogue used for these rate data).
sir_susceptibility_negll = function(par, ylist, R0, gamma, I0, n_sub = 7,
                                    use_prior = TRUE, return_fit = FALSE){
  K = length(ylist); beta = R0 * gamma
  p = .susc_unpack(par, K)

  ll = 0; mu_list = vector("list", K)
  for (s in seq_len(K)){
    y   = ylist[[s]]
    inc = .sir_season_incidence(p$S0[s], I0, beta, gamma, length(y), n_sub)
    mu  = p$c[s] * inc + p$b
    mu_list[[s]] = mu
    ok = !is.na(y)
    if (any(ok)){
      v  = mu[ok] + mu[ok]^2 / p$phi                     # neg-binomial-like observation variance
      ll = ll + sum(dnorm(y[ok], mean = mu[ok], sd = sqrt(v), log = TRUE))
    }
  }

  lp = if (use_prior) .susc_logprior(par, K) else 0

  if (return_fit) return(list(negloglik = -(ll + lp), loglik = ll, params = p, mu = mu_list))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
}

# ---- |-fit the per-season susceptibility model for one country via multi-start optim ----
# ylist: list of weekly ILI+ vectors (one per season, seeded at week 1). R0, infectious period and
# seed I0 come from settings. Returns per-season (S0, c) + shared (b, phi) and the fitted curves.
fit_sir_susceptibility = function(ylist, R0 = 1.5, infectious_period_days = 3, seed_i0 = 1e-5,
                                  n_sub = 7, n_starts = 4, seed = 1){
  gamma = 7 / infectious_period_days                     # per-week recovery rate
  K = length(ylist)
  pos    = lapply(ylist, function(y) y[is.finite(y) & y > 0])
  ymax   = max(vapply(pos, function(y) if (length(y)) max(y)                  else NA_real_, numeric(1)), na.rm = TRUE)
  bguess = max(vapply(pos, function(y) if (length(y)) as.numeric(quantile(y, 0.1)) else NA_real_, numeric(1)), na.rm = TRUE)

  # start in the GROWING regime (S0 high) so the optimiser finds the wave, not a decaying solution
  base = c(rep(qlogis(0.8), K), rep(log(ymax / 0.02), K), log(max(bguess, 1e-3)), log(15))
  names(base) = c(paste0("logit_S0_", seq_len(K)), paste0("log_c_", seq_len(K)), "log_b", "log_phi")
  jit_sd = c(rep(0.8, K), rep(0.5, K), 0.5, 0.5)        # wider spread on S0 to explore growth rates
  set.seed(seed)
  best = NULL
  for (s in seq_len(n_starts)){
    start = base + if (s == 1) 0 else rnorm(length(base), 0, jit_sd)
    fit = tryCatch(
      optim(start, sir_susceptibility_negll, ylist = ylist, R0 = R0, gamma = gamma, I0 = seed_i0,
            n_sub = n_sub, method = "BFGS", control = list(maxit = 200, reltol = 1e-8)),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && (is.null(best) || fit$value < best$value)) best = fit
  }
  if (is.null(best)) stop("fit_sir_susceptibility: all optim starts failed")

  fitted = sir_susceptibility_negll(best$par, ylist = ylist, R0 = R0, gamma = gamma, I0 = seed_i0,
                                    n_sub = n_sub, return_fit = TRUE)
  list(par = best$par, negloglik = best$value, convergence = best$convergence,
       R0 = R0, gamma = gamma, infectious_period_days = infectious_period_days, seed_i0 = seed_i0,
       params = fitted$params, mu = fitted$mu, ylist = ylist)
}

# ---- |-deterministic ILI+ trajectory for one fitted season (smooth curve / projection) ----
# Identical to the fitted mu within the data window; set n_weeks > length(y) to project forward.
sir_susceptibility_trajectory = function(fit, s, n_weeks = length(fit$ylist[[s]]), n_sub = 7){
  inc = .sir_season_incidence(fit$params$S0[s], fit$seed_i0, fit$R0 * fit$gamma, fit$gamma, n_weeks, n_sub)
  fit$params$c[s] * inc + fit$params$b
}

# ---- |-load the committed slim flu ILI+ panel (base R) -> per-season weekly series ----
# Reads data/slim_flu_iliplus.csv (country_short, season, week, season_week, date, value) and
# returns, for one country, an ordered list of weekly ILI+ vectors (one per season, indexed from
# the season start at week 1) plus the season labels and season-week indices -- ready for
# fit_sir_susceptibility(). No tidyverse needed.
load_flu_iliplus_slim = function(country, path = "data/slim_flu_iliplus.csv"){
  d = read.csv(path, stringsAsFactors = FALSE)
  d = d[d$country_short == country, ]
  if (!nrow(d)) stop(sprintf("load_flu_iliplus_slim: no rows for country '%s'", country))
  seasons = sort(unique(d$season))
  by_season = function(col) lapply(seasons, function(ss){ di = d[d$season == ss, ]; di[[col]][order(di$week)] })
  ylist = by_season("value"); names(ylist) = seasons
  weeks = by_season("season_week"); names(weeks) = seasons
  list(country = country, seasons = seasons, ylist = ylist, season_week = weeks)
}

# ---- |-extract one country x season weekly flu ILI+ series from models_in ----
kalman_sir_series = function(models_in, country, season,
                             stream = "ili_plus_sentinel", pathogen = "Influenza"){
  models_in$data_timeseries_long %>%
    filter(indicator == "ili_plus", stream == !!stream, pathogen == !!pathogen,
           agegroup == "age_total", country_short == !!country, season == !!season) %>%
    arrange(date) %>%
    transmute(date, season_week, value)
}
