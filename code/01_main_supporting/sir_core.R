# sir_core.R
#
# Shared, method-agnostic SIR engine + data loaders for the flu ILI+ fits. Every fitting method
# (code/01_main_supporting/methods/) builds on this; the methods themselves only add a likelihood,
# a parameterisation and an optimiser.
#
# Generative backbone (single-population SIR in proportions, shared with the Stan model
# stan/SIR_multiseason_age_vax_2.stan):
#   dS/dt = -beta S I,  dI/dt = beta S I - gamma I,  beta = R0 * gamma
#   weekly ILI+ is proportional to weekly new infections:  E[y_t] = c * (new infections) + b
#   observation noise is neg-binomial-like:  Var(y) = mu + mu^2 / phi
# Everything here is base R.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### SIR integrator ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Extended Kalman Filter pass (state-space methods) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-core EKF pass over one weekly series: returns log-lik + filtered trajectory ----
# Propagates a Gaussian state (S, I, C) through the linearised one-week SIR map with process noise
# qI on I, updating at each observed week. C (within-week cumulative incidence) is reset to 0 each
# week so its weekly value IS the modelled new-infection signal. All parameters arrive on the
# natural scale; beta = R0 * gamma.
#   p0_S, p0_I : initial state sd as a fraction of (S0, I0). Small -> the filter trusts the fitted
#                initial condition (used by the susceptibility EKF, where S0 is the parameter of
#                interest); large -> the filter lets early data move the state freely (tracking).
.ekf_filter = function(y, S0, I0, c, b, phi, beta, gamma, n_sub, qI, p0_S = 0.3, p0_I = 1.0){
  Tn = length(y)

  # state x = (S, I, C); start with C = 0
  x = c(S0, I0, 0)
  P = diag(c((p0_S*S0)^2, (p0_I*I0)^2, 0))            # initial state covariance
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Curve features (method-agnostic descriptors of a weekly ILI+ curve) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Computed from any method's curve mu (a fitted SIR mean OR a smoothed observation), so every method
# reports the same descriptive features. baseline = off-season floor (a low quantile of the curve).

# ---- |-off-season baseline (low quantile of the curve) ----
.curve_baseline = function(mu) as.numeric(quantile(mu[is.finite(mu)], 0.10, names = FALSE))

# ---- |-area under the curve above baseline (seasonal burden) ----
.curve_auc = function(mu, b) sum(pmax(mu - b, 0), na.rm = TRUE)

# ---- |-absolute peak height of the curve ----
.curve_peak_height = function(mu) max(mu, na.rm = TRUE)

# ---- |-steepness = exponential growth rate over the rise (slope of log(above-baseline), onset->peak) ----
# A regression over the whole rising limb (not a single-week jump), so it is robust to local noise
# in the curve. (Phenomenological: it is NOT converted to an SIR S0 -- see method_descriptive.R.)
.curve_steepness = function(mu, season_week, b, onset_frac = 0.1){
  above = mu - b; pkv = max(above, na.rm = TRUE)
  if (!is.finite(pkv) || pkv <= 0) return(NA_real_)
  pk = which.max(mu); on = which(above >= onset_frac * pkv)[1]
  if (is.na(on) || pk - on < 2) return(NA_real_)
  idx = on:pk; z = log(pmax(above[idx], 1e-9)); x = season_week[idx]
  ok = is.finite(z) & is.finite(x)
  if (sum(ok) < 3) return(NA_real_)
  unname(coef(lm(z[ok] ~ x[ok]))[2])
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Data loaders ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-load the committed slim flu ILI+ panel (base R) -> per-season weekly series ----
# Reads data/slim_flu_iliplus.csv (country_short, season, week, season_week, date, value) and
# returns, for one country, an ordered list of weekly ILI+ vectors (one per season, indexed from
# the season start at week 1) plus the season labels and season-week indices -- the panel every
# fitting method consumes. No tidyverse needed.
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
