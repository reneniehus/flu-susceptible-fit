# epidist.R
#
# Delay distributions as first-class objects, plus their discretisation to a day-indexed PMF.
# This is the currency the whole playground trades in: the generation interval, incubation period,
# onset->death and onset->report delays all arrive as `epidist` objects and leave as PMF vectors.
#
# WHY A LOCAL CLASS (and how {epiparameter} plugs in). The brief asks for delay distributions "as
# {epiparameter} objects". {epiparameter} is the right home for a curated, citable epidemiological
# distribution, but it is not installable in every environment (CRAN is blocked here; see
# documentation/decisions.md). So we define a deliberately tiny, {epiparameter}-SHAPED class -- a
# name, a family, and the family's natural parameters -- that (a) runs with base R alone and (b)
# accepts a real <epiparameter> object via `as_epidist()` the moment the package IS available. The
# rest of the code only ever calls `discretise()`, so swapping a hand-set distribution for a curated
# {epiparameter} one is a one-line change at the config, nothing downstream.
#
# DISCRETISATION. A continuous delay X with CDF F is turned into P(delay = d) for integer days d by
# interval censoring, P(d) = F(d + 1/2) - F(d - 1/2), renormalised over the support (Cori et al.
# 2013, AJE, appendix). Two conventions on the day-0 bin:
#   - "interval" (default): keep the [-1/2, 1/2) mass at day 0. Right for incubation / reporting /
#     onset->death delays, where a zero-day delay is real and common.
#   - "cori": force P(0) = 0 and start the support at day 1. Right for the generation interval, where
#     an infector cannot infect on the very day they were infected (matches EpiEstim's w_0 = 0).
# Both renormalise so the returned PMF sums to 1 over its support.
#
# References
#   Cori A, et al. A new framework and software to estimate time-varying reproduction numbers during
#     epidemics. Am J Epidemiol. 2013;178(9):1505-1512.  (discretisation of the serial interval)
#   epiparameter R package (Epiverse-TRACE): the curated-distribution object this class mirrors.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Construct an epidist ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-construct a delay distribution from its family + natural parameters ----
# family is one of "gamma", "lognormal", "weibull". `params` are that family's native parameters:
#   gamma     : shape, rate            (or use epidist_gamma(mean, sd) for the epidemiologist's input)
#   lognormal : meanlog, sdlog
#   weibull   : shape, scale
# `name` is a free-text label ("generation interval", "incubation period", ...) carried for printing.
epidist <- function(name, family, params) {
  family <- match.arg(family, c("gamma", "lognormal", "weibull"))
  needed <- switch(family,
    gamma     = c("shape", "rate"),
    lognormal = c("meanlog", "sdlog"),
    weibull   = c("shape", "scale"))
  miss <- setdiff(needed, names(params))
  if (length(miss)) stop(sprintf("epidist('%s', '%s'): missing parameter(s): %s",
                                 name, family, paste(miss, collapse = ", ")))
  structure(list(name = name, family = family, params = params[needed]),
            class = "epidist")
}

# ---- |-convenience constructors taking the epidemiologist's (mean, sd) input ----
# The natural way to state a delay is "mean X days, sd Y days"; these convert to the family's native
# parameters so the caller never has to do the moment-matching by hand.
epidist_gamma <- function(name, mean, sd) {
  shape <- (mean / sd)^2; rate <- mean / sd^2
  epidist(name, "gamma", list(shape = shape, rate = rate))
}
epidist_lognormal <- function(name, mean, sd) {
  # match the natural-scale mean & sd to meanlog / sdlog
  sdlog   <- sqrt(log1p((sd / mean)^2))
  meanlog <- log(mean) - sdlog^2 / 2
  epidist(name, "lognormal", list(meanlog = meanlog, sdlog = sdlog))
}
epidist_weibull <- function(name, shape, scale) {
  epidist(name, "weibull", list(shape = shape, scale = scale))
}

# ---- |-CDF of an epidist at x (family dispatch) ----
.epidist_cdf <- function(d, x) {
  p <- d$params
  switch(d$family,
    gamma     = pgamma(x,  shape = p$shape, rate  = p$rate),
    lognormal = plnorm(x,  meanlog = p$meanlog, sdlog = p$sdlog),
    weibull   = pweibull(x, shape = p$shape, scale = p$scale))
}

# ---- |-mean of an epidist (for reporting / sanity checks) ----
epidist_mean <- function(d) {
  p <- d$params
  switch(d$family,
    gamma     = p$shape / p$rate,
    lognormal = exp(p$meanlog + p$sdlog^2 / 2),
    weibull   = p$scale * gamma(1 + 1 / p$shape))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Discretise to a day-indexed PMF ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-discretise an epidist to a PMF over integer days 0..max_days ----
# Returns a numeric vector `w` of length max_days + 1 with w[d + 1] = P(delay = d), summing to 1.
# `boundary`:
#   "interval" -> keep the day-0 mass (incubation, reporting, onset->death)
#   "cori"     -> force P(0) = 0 (generation / serial interval; an infector needs >= 1 day)
# `max_days` should be chosen well beyond the distribution's mass; the tail beyond it is folded in by
# the renormalisation. A helper `epidist_support()` picks a sensible default from a high quantile.
discretise <- function(d, max_days = epidist_support(d), boundary = c("interval", "cori")) {
  boundary <- match.arg(boundary)
  days  <- 0:max_days
  # interval-censored mass P(d) = F(d + 1/2) - F(d - 1/2); the day-0 lower edge is clamped at 0
  upper <- .epidist_cdf(d, days + 0.5)
  lower <- .epidist_cdf(d, pmax(days - 0.5, 0))
  w <- upper - lower
  if (boundary == "cori") w[1] <- 0            # no same-day transmission for a generation interval
  s <- sum(w)
  if (s <= 0) stop(sprintf("discretise('%s'): PMF has zero mass -- check max_days / parameters", d$name))
  w / s                                        # renormalise (folds the > max_days tail back in)
}

# ---- |-a sensible discretisation horizon: the 0.999 quantile, rounded up ----
# So `discretise()` can be called without the caller having to know each distribution's tail length.
epidist_support <- function(d, q = 0.999) {
  p <- d$params
  hi <- switch(d$family,
    gamma     = qgamma(q,  shape = p$shape, rate  = p$rate),
    lognormal = qlnorm(q,  meanlog = p$meanlog, sdlog = p$sdlog),
    weibull   = qweibull(q, shape = p$shape, scale = p$scale))
  max(1, ceiling(hi))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Interop with the real {epiparameter} package (extension front) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-coerce a real <epiparameter> object into an epidist (used only when the package is present) ----
# Lets a caller drop a curated, citable distribution straight into a config:
#   library(epiparameter)
#   gi <- as_epidist(epiparameter_db(disease = "COVID-19", epi_name = "generation")[[1]])
# We read the family and its parameters through {epiparameter}'s accessors and rebuild the local
# object, so everything downstream (discretise(), the renewal engine) is unchanged. No-op unless the
# package is installed -- the playground never depends on it.
as_epidist <- function(x, name = NULL) {
  if (inherits(x, "epidist")) return(x)
  if (!requireNamespace("epiparameter", quietly = TRUE))
    stop("as_epidist(): {epiparameter} is not installed. Build the distribution with epidist_gamma()/",
         "epidist_lognormal()/epidist_weibull() instead, or install {epiparameter}.")
  fam <- tolower(as.character(epiparameter::family(x)))
  par <- epiparameter::get_parameters(x)
  nm  <- name %||% tryCatch(as.character(x$epi_name), error = function(e) "delay")
  if (grepl("gamma", fam)) {
    epidist(nm, "gamma", list(shape = unname(par["shape"]), rate = unname(1 / par["scale"])))
  } else if (grepl("lnorm|lognormal", fam)) {
    epidist(nm, "lognormal", list(meanlog = unname(par["meanlog"]), sdlog = unname(par["sdlog"])))
  } else if (grepl("weibull", fam)) {
    epidist(nm, "weibull", list(shape = unname(par["shape"]), scale = unname(par["scale"])))
  } else {
    stop(sprintf("as_epidist(): unsupported {epiparameter} family '%s' (add it to epidist.R).", fam))
  }
}

# ---- |-null-coalescing helper (used above; base R has no %||% before 4.4) ----
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- |-print method: one readable line, mean in days ----
print.epidist <- function(x, ...) {
  pars <- paste(sprintf("%s=%.3g", names(x$params), unlist(x$params)), collapse = ", ")
  cat(sprintf("<epidist> %-20s %s(%s)  mean=%.2f d\n", x$name, x$family, pars, epidist_mean(x)))
  invisible(x)
}
