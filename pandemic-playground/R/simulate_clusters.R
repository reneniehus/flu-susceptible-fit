# simulate_clusters.R
#
# The EARLY-CLUSTER data structure: a handful of independent introductions, each seeding a small,
# mostly self-limiting transmission chain that we observe through its FINAL SIZE. This is the "few
# stuttering chains" regime -- fundamentally different from a single growing incidence curve -- and it
# is the regime in which branching-process inference (phase0_clusters.R) recovers BOTH the reproduction
# number R and the offspring dispersion k. Cluster inference is the only early method that yields k.
#
# THE PROCESS. Each introduction is a Galton-Watson branching process: the index case, and every case
# after it, produces a random number of secondary cases drawn from a NEGATIVE-BINOMIAL offspring
# distribution with mean R and dispersion k (the Lloyd-Smith 2005 superspreading parameterisation;
# k = Inf -> Poisson). The chain grows generation by generation until it dies out (or hits a size cap
# that guards the R >= 1 case). We record each cluster's total size and generational depth.
#
# OBSERVATION. Each case is detected independently with probability `detection_prob`; a cluster with
# zero detected cases is missed entirely. This is realistic and consequential: SINGLETONS (chains of
# size 1) are the easiest clusters to miss, so the observed size distribution is biased -- the exact
# ascertainment problem the zero-truncated estimator in phase0_clusters.R addresses.
#
# Returns the LATENT truth (true cluster sizes, the true R and k) and the OBSERVED (detected) sizes.
#
# References
#   Lloyd-Smith JO, et al. Superspreading and the effect of individual variation on disease emergence.
#     Nature. 2005;438:355-359.
#   Blumberg S, Lloyd-Smith JO. Inference of R0 and transmission heterogeneity from the size
#     distribution of stuttering chains. PLoS Comput Biol. 2013;9(5):e1002993.
#   Endo A, Abbott S, Kucharski AJ, Funk S. Estimating the overdispersion in COVID-19 transmission
#     using outbreak sizes outside China. Wellcome Open Res. 2020;5:67.

# ---- |-simulate independent introduction clusters as Galton-Watson branching processes ----
# cfg : the config (uses cfg$clusters). Uses the current RNG state (seeded by the caller).
simulate_clusters <- function(cfg) {
  cl   <- cfg$clusters
  R    <- cl$R; k <- cl$dispersion_k; n <- cl$n_introductions; maxN <- cl$max_size %||% 1000

  size <- integer(n); gens <- integer(n); capped <- logical(n)
  for (i in seq_len(n)) {
    total   <- 0L                                       # cases accumulated in this chain
    current <- 1L                                       # the index case starts the chain
    g       <- 0L
    while (current > 0L && total < maxN) {
      total   <- total + current
      # offspring of the `current` active cases (NB per case; Poisson when k = Inf)
      current <- if (is.finite(k)) sum(stats::rnbinom(current, size = k, mu = R))
                 else               sum(stats::rpois(current, R))
      g <- g + 1L
    }
    size[i] <- total; gens[i] <- g; capped[i] <- total >= maxN
  }

  # observation: thin each cluster's cases by detection; a cluster with 0 detected is not seen at all
  detected       <- stats::rbinom(n, size, cl$detection_prob)
  observed_sizes <- detected[detected > 0]

  list(
    true_sizes     = size,                              # complete cluster sizes (LATENT)
    generations    = gens,
    capped         = capped,                            # TRUE where the size cap bound (only if R >= 1)
    observed_sizes = observed_sizes,                    # detected sizes (missed clusters dropped)
    R = R, k = k, n_introductions = n,
    n_observed = length(observed_sizes), detection_prob = cl$detection_prob
  )
}
