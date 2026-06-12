# model_kalman_sir.R
#
# An R re-implementation of the Stan SIR model (stan/SIR_multiseason_age_vax_2.stan), fitting
# flu ILI+ with an Extended Kalman Filter (EKF) instead of HMC. Everything is in base R.
#
# Idea (same generative story as the Stan model, single-population form):
#   - latent SIR in population proportions:  dS/dt = -beta S I,  dI/dt = beta S I - gamma I
#   - weekly ILI+ is proportional to weekly new infections:  E[y_t] = c * (new infections in week t)
#     (c folds in the infection->ILI fraction "prop_ili", positivity scaling and the per-100k rate)
#   - observation noise is neg-binomial-like:  Var(y) = mu + mu^2 / phi   (matches the Stan obs model)
# The EKF linearises the one-week SIR map to propagate a Gaussian state (S, I, C), where C is the
# within-week cumulative incidence (reset to 0 each week so the weekly increment IS the observation).
# Parameters are fit by maximising the EKF log-likelihood with optim().
#
# Fixed by choice (kept simple / identifiable): gamma (recovery rate) -> mean infectious period.
# Fitted: R0, initial susceptible S0, initial infectious I0, obs scale c, overdispersion phi,
#         process-noise sd on I (q_I).

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

# ---- |-EKF negative log-(posterior) for one weekly ILI+ series ----
# y: numeric vector (NA allowed for missing weeks). gamma in per-week units. Returns -loglik
# (plus -logprior when use_prior=TRUE, i.e. the negative log-posterior used for MAP fitting).
ekf_sir_negloglik = function(par, y, gamma, n_sub = 7, use_prior = TRUE, return_fit = FALSE){
  p = .kf_unpack(par)
  beta = p$R0 * gamma
  Tn = length(y)

  # state x = (S, I, C); start with C = 0
  x = c(p$S0, p$I0, 0)
  P = diag(c((0.3*p$S0)^2, (1.0*p$I0)^2, 0))          # initial state covariance
  Q = diag(c(0, p$qI^2, 0))                            # process noise (on I)

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

    mu = p$c * xpr[3] + p$b                            # expected observed ILI+ (epidemic + baseline)
    mu_pred[t] = mu

    if (!is.na(y[t])){
      H   = matrix(c(0, 0, p$c), nrow = 1)
      Rt  = mu + mu^2 / p$phi                          # neg-binomial-like observation variance
      Sinn = as.numeric(H %*% Ppr %*% t(H)) + Rt
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

  lp = if (use_prior) .kf_logprior(par) else 0

  if (return_fit) return(list(negloglik = -(ll + lp), loglik = ll, mu_pred = mu_pred,
                              S = S_filt, I = I_filt, params = p))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
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

# ---- |-extract one country x season weekly flu ILI+ series from models_in ----
kalman_sir_series = function(models_in, country, season,
                             stream = "ili_plus_sentinel", pathogen = "Influenza"){
  models_in$data_timeseries_long %>%
    filter(indicator == "ili_plus", stream == !!stream, pathogen == !!pathogen,
           agegroup == "age_total", country_short == !!country, season == !!season) %>%
    arrange(date) %>%
    transmute(date, season_week, value)
}
