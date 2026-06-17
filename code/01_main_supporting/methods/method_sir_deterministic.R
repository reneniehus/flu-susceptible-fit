# method_sir_deterministic.R
#
# METHOD: deterministic per-season susceptibility fit (registry name "deterministic").
# Requires sir_core.R (the SIR integrator) to be sourced first.
#
# Generative model: the season's ILI+ wave IS a single deterministic SIR; observations are
# overdispersed noise around it. S0 has to be earned by the dynamics -- the wave's rise rate,
# r = gamma*(R0*S0 - 1), identifies it. (No filter; misfit shows up as honest residuals.)
#
# Fixed by assumption (from settings): R0, gamma (infectious period) and the seed I0 (a small
# constant import, planted at week 1 of every season). The ABSOLUTE S0 is conditional on R0/I0;
# the RELATIVE S0 across seasons is the quantity of interest.
#
# Fitted, per country (K seasons):
#   per season : S0[1..K] (susceptibility) and c[1..K] (reporting fraction)
#   shared     : b (additive off-season baseline) and phi (observation overdispersion)
#
# Returns the common method-fit contract consumed by methods_registry.R:
#   list(method, R0, gamma, seed_i0, params = list(S0, c, b, phi, qI), mu = <per-season curves>,
#        ylist, convergence, negloglik)

# ---- |-unpack the parameter vector into natural-scale params (K seasons) ----
.det_unpack = function(par, K){
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
.det_logprior = function(par, K){
  sum(dnorm(par[1:K], qlogis(0.75), 1.0, log = TRUE)) +  # logit_S0 per season
  dnorm(par[[2*K + 2]], log(15), 0.8, log = TRUE)        # shared log_phi
}

# ---- |-deterministic negative log-(posterior) over a country's seasons (fixed R0, I0) ----
# mu = c[s] * (weekly new infections) + b, with overdispersed obs noise Var = mu + mu^2/phi.
.det_negll = function(par, ylist, R0, gamma, I0, n_sub = 7, use_prior = TRUE, return_fit = FALSE){
  K = length(ylist); beta = R0 * gamma
  p = .det_unpack(par, K)

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

  lp = if (use_prior) .det_logprior(par, K) else 0

  if (return_fit) return(list(negloglik = -(ll + lp), loglik = ll, params = p, mu = mu_list))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
}

# ---- |-fit the deterministic method for one country via multi-start optim ----
# ylist: list of weekly ILI+ vectors (one per season, seeded at week 1). Returns the common
# method-fit contract.
fit_sir_deterministic = function(ylist, R0 = 1.5, infectious_period_days = 3, seed_i0 = 1e-5,
                                 n_sub = 7, n_starts = 4, seed = 1, ...){   # ... ignores cross-method args
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
      optim(start, .det_negll, ylist = ylist, R0 = R0, gamma = gamma, I0 = seed_i0,
            n_sub = n_sub, method = "BFGS", control = list(maxit = 200, reltol = 1e-8)),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && (is.null(best) || fit$value < best$value)) best = fit
  }
  if (is.null(best)) stop("fit_sir_deterministic: all optim starts failed")

  fitted = .det_negll(best$par, ylist = ylist, R0 = R0, gamma = gamma, I0 = seed_i0,
                      n_sub = n_sub, return_fit = TRUE)
  fitted$params$qI = NA_real_                            # no process noise (schema uniformity)
  list(method = "deterministic", R0 = R0, gamma = gamma, seed_i0 = seed_i0,
       params = fitted$params, mu = fitted$mu, ylist = ylist,
       convergence = best$convergence, negloglik = best$value)
}

# ---- |-deterministic ILI+ trajectory for one fitted season (smooth curve / projection) ----
# Identical to the fitted mu within the data window; set n_weeks > length(y) to project forward.
sir_deterministic_trajectory = function(fit, s, n_weeks = length(fit$ylist[[s]]), n_sub = 7){
  inc = .sir_season_incidence(fit$params$S0[s], fit$seed_i0, fit$R0 * fit$gamma, fit$gamma, n_weeks, n_sub)
  fit$params$c[s] * inc + fit$params$b
}
