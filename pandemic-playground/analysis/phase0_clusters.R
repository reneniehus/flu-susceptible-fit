# phase0_clusters.R
#
# PHASE 0 -- before local introduction, when the pathogen is a FEW STUTTERING CHAINS, not one growing
# curve.
# QUESTION: "Is it spreading person-to-person, how fast -- AND how heterogeneous (superspreading)?"
#
# THE METHOD (branching-process inference from cluster final sizes). When all you have is a handful of
# small, mostly self-limiting introduction clusters, exponential-growth fitting is premature and can
# actively mislead (overdispersion makes the early curve stochastic and biases a naive growth-rate R
# DOWNWARD). The right tool is the size DISTRIBUTION of the chains, which -- via the Galton-Watson
# theory -- identifies BOTH the reproduction number R and the offspring dispersion k. Cluster
# inference is the ONLY early method that yields k, and k in turn tells you how much to trust any later
# growth-rate R.
#
#   * Crude R from the mean size. For a sub-critical process the expected chain size is 1/(1-R) for
#     ANY offspring distribution, so R_crude = 1 - 1/mean(size). (Mean pins R; the spread pins k.)
#   * Joint (R, k) by maximum likelihood on the exact chain-size PMF. For negative-binomial offspring
#     (mean R, dispersion k), the total-progeny PMF follows from the Dwass / Otter identity
#         P(N = n) = (1/n) P( X_1 + ... + X_n = n - 1 ),
#     and the sum of n iid NB(size k) offspring is NB(size n*k), giving a closed form (k = Inf reduces
#     to the Borel(R) distribution). Many singletons plus a rare large cluster => small k.
#   * The singleton fraction gives a second equation for k: P(size = 1) = P(0 offspring) = (1 + R/k)^-k.
#   * Extinction probability q solves q = (1 + (R/k)(1 - q))^-k (q = 1 for R <= 1).
#
# THE CAVEAT THAT BITES HARDEST. Extinct minor clusters cannot, alone, tell a sub-critical R from its
# super-critical conjugate: a super-critical process conditioned on extinction is distributionally
# identical to a sub-critical one. Finite clusters therefore identify the SUB-CRITICAL branch; to
# recover a super-critical R you need an external anchor -- whether a major outbreak occurred, or the
# fraction of introductions that took off (Nishiura et al. 2012). And detection matters: singletons
# are the easiest clusters to miss, so a zero-truncated fit is offered for that bias.
#
# References
#   Blumberg S, Lloyd-Smith JO. PLoS Comput Biol. 2013;9(5):e1002993.  (R and k from chain sizes)
#   Endo A, et al. Wellcome Open Res. 2020;5:67.  (COVID R and k from exported-case clusters)
#   Nishiura H, et al. J Theor Biol. 2012;294:48-55.  (super-critical R from minor-outbreak sizes)
#
# Requires: analysis_common.R.

# ---- |-log chain-size PMF for NB(mean R, dispersion k) offspring (Dwass/Otter; Borel if k = Inf) ----
.chain_logpmf <- function(n, R, k) {
  if (!is.finite(k)) return(-R * n + (n - 1) * log(R * n) - lgamma(n + 1))    # Borel(R) limit
  nk <- n * k
  -log(n) + lgamma(n - 1 + nk) - lgamma(nk) - lgamma(n) +
    nk * log(k / (k + R)) + (n - 1) * log(R / (k + R))
}

# ---- |-joint MLE of R and k from a vector of cluster final sizes ----
# sizes : integer chain sizes (>= 1). Returns R, k with 95% CIs (from the observed information), the
# crude mean-based R, and a super-critical flag. `zero_truncate = m` conditions on size >= m (m = 2
# removes the singleton-under-ascertainment bias when singletons are missed).
cluster_size_fit <- function(sizes, level = 0.95, zero_truncate = 1) {
  sizes <- sizes[is.finite(sizes) & sizes >= zero_truncate]
  if (length(sizes) < 5) stop("cluster_size_fit: need >= 5 clusters to fit")
  # truncation normaliser: log P(N >= zero_truncate) so the likelihood is conditioned correctly
  ltrunc <- function(R, k) if (zero_truncate <= 1) 0 else
    log1p(-sum(exp(vapply(seq_len(zero_truncate - 1), .chain_logpmf, numeric(1), R = R, k = k))))
  nll <- function(par) {
    R <- exp(par[1]); k <- exp(par[2])
    -(sum(vapply(sizes, .chain_logpmf, numeric(1), R = R, k = k)) - length(sizes) * ltrunc(R, k))
  }
  R0 <- min(max(1 - 1 / mean(sizes), 0.3), 0.95)
  fit <- stats::optim(c(log(R0), log(0.5)), nll, method = "L-BFGS-B",
                      lower = c(log(0.05), log(0.02)), upper = c(log(4), log(50)), hessian = TRUE)
  est <- exp(fit$par); R <- est[1]; k <- est[2]
  se <- tryCatch(sqrt(diag(solve(fit$hessian))), error = function(e) c(NA, NA))   # SE on the log scale
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(R = R, k = k,
       R_ci = R * exp(c(-1, 1) * z * se[1]), k_ci = k * exp(c(-1, 1) * z * se[2]),
       R_crude = 1 - 1 / mean(sizes), n = length(sizes),
       supercritical_possible = R >= 0.95, convergence = fit$convergence)
}

# ---- |-extinction probability of a chain: q solves q = (1 + (R/k)(1-q))^-k  (q = 1 if R <= 1) ----
cluster_extinction <- function(R, k) {
  if (R <= 1) return(1)
  g <- if (!is.finite(k)) function(q) exp(-R * (1 - q)) else function(q) (1 + (R / k) * (1 - q))^(-k)
  q <- 0.1; for (i in 1:200) q <- g(q); q                        # fixed-point iteration
}

# ---- |-second estimate of k from the singleton fraction, given R (moment method A3) ----
# P(size = 1) = P(0 offspring) = (1 + R/k)^-k ; solve for k. Sensitive to missed singletons -- pair
# with zero-truncation when detection is imperfect.
cluster_singleton_k <- function(sizes, R) {
  f1 <- mean(sizes == 1)
  if (f1 <= 0 || f1 >= 1) return(NA_real_)
  if (exp(-R) >= f1) return(Inf)                                 # more singletons than Poisson -> ~no overdispersion
  stats::uniroot(function(k) (1 + R / k)^(-k) - f1, c(0.01, 100))$root
}

# ---- |-Phase-0 cluster analysis from the observed cluster sizes ----
cluster_analysis <- function(input, zero_truncate = 1, level = 0.95) {
  if (is.null(input$clusters) || !length(input$clusters$sizes))
    stop("cluster_analysis: no observed cluster sizes in the input")
  sizes <- input$clusters$sizes
  fit   <- cluster_size_fit(sizes, level = level, zero_truncate = zero_truncate)
  fit$singleton_fraction <- mean(sizes == 1)
  fit$k_from_singletons  <- cluster_singleton_k(sizes, fit$R)
  fit$extinction_prob    <- cluster_extinction(fit$R, fit$k)
  fit
}

# ---- |-score the cluster fit against the true R and dispersion k ----
score_clusters <- function(sim, cf) {
  tr <- truth_clusters(sim)
  list(R = cf$R, true_R = tr$R, R_in_ci = tr$R >= cf$R_ci[1] && tr$R <= cf$R_ci[2],
       k = cf$k, true_k = tr$k, k_in_ci = is.finite(tr$k) && tr$k >= cf$k_ci[1] && tr$k <= cf$k_ci[2],
       R_crude = cf$R_crude)
}
