# simulate_source.R
#
# Stage 1: the source epidemic at location X. A stochastic renewal process (renewal.R) seeded by a
# handful of index infections on day 0. If the variant is enabled, a second, fitter strain is planted
# on its introduction day and runs alongside the wild type, sharing X's susceptible pool. This stage
# owns the epidemic that everything else is downstream of: its infectious PREVALENCE drives exports,
# and its VARIANT FREQUENCY sets the strain mix that travellers carry abroad.
#
# Returns the LATENT truth for X: infections by day and strain, infectious prevalence (count and
# fraction of X's population), the susceptible trajectory, and the variant frequency series.
#
# Requires: R/renewal.R, R/utils.R sourced first.

# ---- |-simulate the source epidemic at X (one or two strains) ----
# cfg : the config; par : draw_parameters(cfg). Uses the current RNG state (seeded by the caller).
simulate_source <- function(cfg, par) {
  n   <- cfg$n_days
  gi  <- par$pmf$gi
  N   <- cfg$source$population
  K   <- par$n_strains

  # ---- per-strain reproduction numbers: the variant inherits X's R trajectory x (1 + fitness) ----
  rt_list <- list(cfg$rt_source)
  if (K == 2L) {
    variant_rt <- cfg$rt_source
    variant_rt$value <- variant_rt$value * (1 + cfg$variant$fitness)
    rt_list[[2]] <- variant_rt
  }

  # ---- seeding: wild-type index cases on day 0; the variant on its introduction day ----
  seeding <- matrix(0, nrow = n, ncol = K)
  seeding[1, 1] <- cfg$source_seed_infections
  if (K == 2L && !is.na(par$variant_day))
    seeding[par$variant_day + 1, 2] <- cfg$variant$seed_infections

  sim <- simulate_renewal(n, N, gi, rt_list, seeding, stochastic = TRUE)

  total_inc  <- rowSums(sim$incidence)
  variant_freq <- if (K == 2L) {
    denom <- sim$incidence[, 1] + sim$incidence[, 2]
    ifelse(denom > 0, sim$incidence[, 2] / denom, NA_real_)   # NA before any infections exist
  } else rep(0, n)

  list(
    infections        = sim$incidence,                 # n x K (columns: wild type, variant)
    infections_total  = total_inc,
    prevalence        = sim$prevalence,                # currently-infectious COUNT at X
    prevalence_frac   = sim$prevalence / N,            # as a fraction of X's population (drives exports)
    prevalence_by_strain = sim$prevalence_by_strain,
    susceptible       = sim$susceptible,
    variant_freq      = variant_freq                   # true variant frequency among new infections at X
  )
}
