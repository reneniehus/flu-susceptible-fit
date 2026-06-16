# method_sir_ekf.R
#
# METHOD: Extended Kalman Filter per-season susceptibility fit (registry name "ekf").
# Requires sir_core.R (.ekf_filter, the SIR integrator) to be sourced first.
#
# Same generative model and fixed assumptions as the deterministic method (fixed R0, gamma and seed
# I0; per-season S0 and reporting fraction c; shared baseline b and overdispersion phi), but the
# latent epidemic is allowed to WANDER off the deterministic SIR: a state-space model with PROCESS
# NOISE q_I on the infectious compartment, fitted as ONE shared value per country. The likelihood is
# the EKF prediction-error decomposition (the state is integrated out by the filter).
#
# Why this still identifies S0 (unlike a loose tracking filter): R0 and the seed I0 are FIXED, the
# initial state covariance is TIGHT (p0_S, p0_I small -> the filter trusts the fitted S0), and q_I is
# regularised small. So the deterministic backbone driven by S0 (the wave's rise rate) carries the
# signal, while q_I only makes modest weekly corrections. Epidemics genuinely have process noise --
# this method states that honestly while keeping S0 interpretable.
#
# Returns the common method-fit contract (see methods_registry.R). The fitted curve fit$mu is the
# one-step-ahead FILTERED mean (the flexible "red curve"); fit$params$qI is the fitted process noise.

# ---- |-tight initial-state covariance (fraction of S0, I0): trust the fitted IC, let q_I do the wandering ----
.EKF_P0_S = 0.05
.EKF_P0_I = 0.05

# ---- |-unpack the EKF parameter vector into natural-scale params (K seasons) ----
.ekf_unpack = function(par, K){
  list(
    S0  = plogis(par[1:K]),                              # per-season susceptibility (0,1)
    c   = exp(par[(K+1):(2*K)]),                          # per-season reporting fraction / obs scale
    b   = exp(par[[2*K + 1]]),                            # shared additive off-season baseline
    phi = exp(par[[2*K + 2]]),                            # shared observation overdispersion
    qI  = exp(par[[2*K + 3]])                             # shared process-noise sd on I (one per country)
  )
}

# ---- |-weakly-informative log-priors (regularise S0, phi and q_I; c, b left free) ----
# The q_I prior is what keeps the filter honest: centred small so it cannot inflate to absorb the
# wave (which would mask S0). S0 and phi as in the deterministic method.
.ekf_logprior = function(par, K){
  sum(dnorm(par[1:K], qlogis(0.75), 1.0, log = TRUE)) +  # logit_S0 per season
  dnorm(par[[2*K + 2]], log(15),   0.8, log = TRUE) +    # shared log_phi
  dnorm(par[[2*K + 3]], log(1e-4), 1.0, log = TRUE)      # shared log_qI (kept small)
}

# ---- |-EKF negative log-(posterior) over a country's seasons (fixed R0, I0; shared q_I) ----
.ekf_negll = function(par, ylist, R0, gamma, I0, n_sub = 7, use_prior = TRUE, return_fit = FALSE){
  K = length(ylist); beta = R0 * gamma
  p = .ekf_unpack(par, K)

  ll = 0; mu_list = vector("list", K)
  for (s in seq_len(K)){
    f = .ekf_filter(ylist[[s]], p$S0[s], I0, p$c[s], p$b, p$phi, beta, gamma, n_sub, p$qI,
                    p0_S = .EKF_P0_S, p0_I = .EKF_P0_I)
    ll = ll + f$loglik; mu_list[[s]] = f$mu_pred         # filtered one-step-ahead mean (the red curve)
  }

  lp = if (use_prior) .ekf_logprior(par, K) else 0

  if (return_fit) return(list(negloglik = -(ll + lp), loglik = ll, params = p, mu = mu_list))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
}

# ---- |-fit the EKF method for one country via multi-start optim ----
# Same signature as the other methods (so the registry can call it uniformly).
fit_sir_ekf = function(ylist, R0 = 1.5, infectious_period_days = 3, seed_i0 = 1e-5,
                       n_sub = 7, n_starts = 4, seed = 1){
  gamma = 7 / infectious_period_days                     # per-week recovery rate
  K = length(ylist)
  pos    = lapply(ylist, function(y) y[is.finite(y) & y > 0])
  ymax   = max(vapply(pos, function(y) if (length(y)) max(y)                  else NA_real_, numeric(1)), na.rm = TRUE)
  bguess = max(vapply(pos, function(y) if (length(y)) as.numeric(quantile(y, 0.1)) else NA_real_, numeric(1)), na.rm = TRUE)

  # start in the GROWING regime (S0 high) with a small process noise
  base = c(rep(qlogis(0.8), K), rep(log(ymax / 0.02), K), log(max(bguess, 1e-3)), log(15), log(1e-4))
  names(base) = c(paste0("logit_S0_", seq_len(K)), paste0("log_c_", seq_len(K)),
                  "log_b", "log_phi", "log_qI")
  jit_sd = c(rep(0.8, K), rep(0.5, K), 0.5, 0.5, 0.8)    # wider spread on S0; some spread on q_I
  set.seed(seed)
  best = NULL
  for (s in seq_len(n_starts)){
    start = base + if (s == 1) 0 else rnorm(length(base), 0, jit_sd)
    fit = tryCatch(
      optim(start, .ekf_negll, ylist = ylist, R0 = R0, gamma = gamma, I0 = seed_i0,
            n_sub = n_sub, method = "BFGS", control = list(maxit = 200, reltol = 1e-8)),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && (is.null(best) || fit$value < best$value)) best = fit
  }
  if (is.null(best)) stop("fit_sir_ekf: all optim starts failed")

  fitted = .ekf_negll(best$par, ylist = ylist, R0 = R0, gamma = gamma, I0 = seed_i0,
                      n_sub = n_sub, return_fit = TRUE)
  list(method = "ekf", R0 = R0, gamma = gamma, seed_i0 = seed_i0,
       params = fitted$params, mu = fitted$mu, ylist = ylist,
       convergence = best$convergence, negloglik = best$value)
}

# ---- |-deterministic ILI+ trajectory for one fitted season (mechanistic mean / projection) ----
# The EKF's process noise is zero-mean, so the model's expected curve is the deterministic SIR from
# the fitted S0 (used for projection; fit$mu holds the filtered in-sample curve). n_weeks > length(y)
# projects forward.
sir_ekf_trajectory = function(fit, s, n_weeks = length(fit$ylist[[s]]), n_sub = 7){
  inc = .sir_season_incidence(fit$params$S0[s], fit$seed_i0, fit$R0 * fit$gamma, fit$gamma, n_weeks, n_sub)
  fit$params$c[s] * inc + fit$params$b
}
