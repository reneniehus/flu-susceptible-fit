# simulate_importations.R
#
# Stage 3: importation. Each day, travellers leaving X carry infection in proportion to X's current
# infectious prevalence, so the expected number of infected arrivals in country c is
#   E[imports_{c,t}] = volume_{X->c, t} * prevalence_fraction_{X, t}
# and the realised imports are a Poisson draw around it. When the variant is active, each imported
# infection is labelled wild-type or variant by a binomial split on X's current variant prevalence --
# so the variant travels abroad exactly in proportion to the connectivity and the source dynamics,
# with no separate per-country seeding. These imported infections are the seeds of the local
# epidemics (Stage 4).
#
# We record the TRUE imports per country per day and strain. What an analyst can actually SEE -- the
# DETECTED imports, thinned by surveillance quality -- is produced later in observe.R; keeping the
# true seeding here separate from its detection there is again the truth/observation split.
#
# Requires: R/renewal.R (upstream), R/utils.R.

# ---- |-simulate true imported infections into every country (by day and strain) ----
# cfg : config; par : draw_parameters(cfg); source : simulate_source() output; flights : simulate_flights().
# Uses the current RNG state.
simulate_importations <- function(cfg, par, source, flights) {
  n  <- cfg$n_days
  cc <- cfg$countries
  K  <- par$n_strains

  prev_frac <- source$prevalence_frac                              # length n: fraction of X infectious
  # expected infected arrivals: true volume (n x c) scaled row-wise by the source prevalence fraction
  expected <- sweep(flights$true, 1, prev_frac, `*`)               # n x c

  # realised imports ~ Poisson(expected), drawn per cell
  imports_total <- matrix(stats::rpois(length(expected), as.vector(expected)),
                          nrow = n, dimnames = list(NULL, cc))

  # split each country-day's imports into strains by X's current variant prevalence fraction
  imports <- array(0, dim = c(n, length(cc), K), dimnames = list(NULL, cc, NULL))
  if (K == 1L) {
    imports[, , 1] <- imports_total
  } else {
    p_variant <- ifelse(source$prevalence > 0,
                        source$prevalence_by_strain[, 2] / source$prevalence, 0)   # length n
    for (t in seq_len(n)) {
      if (sum(imports_total[t, ]) == 0) next                                     # no imports that day -> nothing to split
      n_var <- stats::rbinom(length(cc), imports_total[t, ], p_variant[t])          # variant arrivals
      imports[t, , 2] <- n_var
      imports[t, , 1] <- imports_total[t, ] - n_var
    }
  }

  list(imports = imports,                                           # n x n_countries x K (TRUE seeds)
       imports_total = imports_total,                               # n x n_countries
       expected = expected)                                         # n x n_countries (expected, for reference)
}
