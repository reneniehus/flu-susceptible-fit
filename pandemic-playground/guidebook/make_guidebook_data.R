# Generate the data behind the visual guidebook: run one default pandemic and, at each outbreak
# stage, capture the surveillance data an analyst would see plus what each tool returns, scored
# against truth. Writes a single JSON the HTML guidebook embeds. Run from the project root.
suppressWarnings(suppressMessages(source("setup.R")))
library(jsonlite)

set.seed(1)
sim   <- simulate_pandemic(default_config())
input <- as_analysis_input(sim)
cfg   <- sim$config
country <- "IT"                                        # the worked destination country
src <- cfg$source$code

round2 <- function(x, d = 3) ifelse(is.finite(x), round(x, d), NA)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 0 -- before local introduction (the source X) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

src_cases <- loc_series(input$cases_by_report, src)
src_inf   <- truth_infections(sim, src)
win0 <- 20:45
ga <- growth_analysis(input, src, window = win0)
sg <- score_growth(sim, ga)

# exponential fit line for plotting (over the fit window, extended a little)
gfit_days <- min(win0):(max(win0) + 6)
d0 <- src_cases[src_cases$day %in% win0, ]
gfit <- stats::glm(cases ~ day, data = d0, family = poisson())
gfit_line <- data.frame(day = gfit_days, fit = as.numeric(exp(predict(gfit, newdata = data.frame(day = gfit_days)))))

phase0_curve <- data.frame(day = src_cases$day,
                           reported = src_cases$cases,
                           infections = src_inf$infections[match(src_cases$day, src_inf$day)])
phase0_curve <- phase0_curve[phase0_curve$day <= 60, ]

# importation-risk regression across countries
ir <- importation_risk(input, window = 25:70, min_surveillance = 0.8)
sir <- score_importation_risk(ir)
ir_out <- ir[, c("country","volume","surveillance_quality","detected_imports","expected",
                 "pi_lower","pi_upper","under_detecting")]

# catchment back-calculation, reported as an honest range across the traveller/detection corrections
gi_mean0 <- epidist_mean(cfg$delays$generation_interval)
cb  <- catchment_backcalc(input, 30:60, min_surveillance = 0.75, source_pop = cfg$source$population)
scb <- score_catchment(sim, cb, input)
crg <- catchment_range(input, 30:60, source_pop = cfg$source$population,
                       growth_rate = ga$r, recovery_rate = 1 / gi_mean0, min_surveillance = 0.75)

# cluster / branching-process inference: recover R AND the dispersion k from transmission-chain sizes
cf  <- cluster_analysis(input)
scl <- score_clusters(sim, cf)
csizes <- input$clusters$sizes
ctab   <- table(csizes)
cluster_obs <- data.frame(size = as.integer(names(ctab)), share = as.numeric(ctab) / length(csizes))
cluster_fit <- data.frame(size = 1:max(csizes),
                          pmf  = exp(vapply(1:max(csizes), function(n) .chain_logpmf(n, cf$R, cf$k), numeric(1))))

phase0 <- list(
  curve = phase0_curve,
  growth_fit = gfit_line,
  growth = list(r = round2(ga$r), doubling = round2(ga$doubling_time, 1),
                R = round2(ga$R, 2), R_lower = round2(ga$R_ci[1], 2), R_upper = round2(ga$R_ci[2], 2),
                truth_R = round2(sg$truth_R, 2), window = range(win0)),
  clusters = list(
    R = round2(cf$R, 2), R_lower = round2(cf$R_ci[1], 2), R_upper = round2(cf$R_ci[2], 2),
    k = round2(cf$k, 2), k_lower = round2(cf$k_ci[1], 2), k_upper = round2(cf$k_ci[2], 2),
    truth_R = round2(scl$true_R, 2), truth_k = round2(scl$true_k, 2),
    R_in_ci = scl$R_in_ci, k_in_ci = scl$k_in_ci,
    n = cf$n, singleton_frac = round2(mean(csizes == 1), 2),
    extinction = round2(cf$extinction_prob, 2), max_size = max(csizes),
    obs = cluster_obs, fit = cluster_fit),
  importation = ir_out,
  importation_flagged = sum(ir$under_detecting),
  importation_cor = round2(sir$cor_residual_surveillance, 2),
  catchment = list(est_prevalence = round(cb$est_prevalence), true_prevalence = round(scb$true_prevalence),
                   ratio = round2(scb$ratio, 2), under_ascertainment = round2(scb$true_underascertainment, 1),
                   reported = scb$reported_source_cases,
                   low = round(crg$low), central = round(crg$central), high = round(crg$high))
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 1 -- early local growth (country) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

as_of1 <- 100
by_report <- loc_series(input$cases_by_report, country)
by_onset  <- loc_series(input$cases_by_onset, country)          # right-truncated
truth_onset <- truth_cases_by_onset(sim, country)              # eventual

nc <- nowcast_truncation(input, country, as_of = as_of1)
nc_recent <- nc[nc$onset_day > (as_of1 - 21) & nc$onset_day <= as_of1, ]
nc_recent$truth <- truth_onset$cases[match(nc_recent$onset_day, truth_onset$day)]
sn <- score_nowcast(sim, nc, recent = 14)

# epidemic curve for the country up to as_of (by report date), plus the nowcast overlay on the tail
phase1_curve <- by_report[by_report$day <= as_of1, c("day","cases")]

# Rt over time
rt <- rt_analysis(input, country)
rt_truth <- truth_rt(sim, country)
rt$truth <- rt_truth$Rt_effective[match(rt$day, rt_truth$day)]
sr <- score_rt(sim, rt)

# delay-adjusted CFR over time
cfr_cut <- seq(50, 180, by = 10)
cfr_roll <- do.call(rbind, lapply(cfr_cut, function(t) {
  s <- suppressWarnings(cfr_static(input, country, as_of = t))
  sk <- score_cfr(sim, s, input)
  data.frame(as_of = t, naive = s$cfr_naive, adjusted = s$cfr_adjusted,
             lower = s$cfr_lower, upper = s$cfr_upper, truth = sk$true_confirmed_cfr)
}))

phase1 <- list(
  as_of = as_of1,
  curve = phase1_curve,
  nowcast = nc_recent[, c("onset_day","observed","nowcast","nowcast_lower","nowcast_upper","truth","completeness")],
  nowcast_improvement = round2(sn$rmse_improvement * 100, 0),
  rt = data.frame(day = rt$day, Rt = round2(rt$Rt, 2), lower = round2(rt$Rt_lower, 2),
                  upper = round2(rt$Rt_upper, 2), truth = round2(rt$truth, 2)),
  rt_score = list(cor = round2(sr$cor, 2), mae = round2(sr$mae, 2)),
  cfr = data.frame(as_of = cfr_roll$as_of, naive = round2(cfr_roll$naive, 4),
                   adjusted = round2(cfr_roll$adjusted, 4), lower = round2(cfr_roll$lower, 4),
                   upper = round2(cfr_roll$upper, 4), truth = round2(cfr_roll$truth, 4)),
  true_ifr = cfg$ifr
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 2 -- growth to peak, healthcare demand ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

as_of2 <- 95
cap <- truth_capacity(sim, country)
pop <- sim$par$country_pop[[country]]
fc <- forecast_capacity(input, country, as_of = as_of2, horizon = 30, capacity = cap, population = pop)
sf <- score_forecast(sim, fc, input)

adm_obs <- loc_series(input$admissions, country)
adm_hist <- adm_obs[adm_obs$day <= as_of2, c("day","admissions")]
adm_truth_future <- adm_obs[adm_obs$day > as_of2 & adm_obs$day <= as_of2 + 30, c("day","admissions")]
# reshape the scenario fan wide
proj <- fc$projection
fan <- reshape(proj[, c("day","scenario","admissions")], idvar = "day", timevar = "scenario", direction = "wide")
names(fan) <- gsub("admissions.R x ", "R", names(fan)); names(fan) <- gsub(" ", "", names(fan))

# intervention ITS
iday <- cfg$rt_country[[country]]$t[2]
it <- intervention_its(input, country, intervention_day = iday, window = 21)
sv <- score_intervention(sim, it, window = 21)
its_curve <- by_onset[by_onset$day >= (it$breakpoint - 21) & by_onset$day <= (it$breakpoint + 21), c("day","cases")]

phase2 <- list(
  as_of = as_of2, capacity = round(cap),
  admissions_history = adm_hist,
  admissions_truth = adm_truth_future,
  forecast = fan,
  R_now = round2(fc$R_now, 2),
  breach = list(forecast = sf$forecast_breach, forecast_day = sf$forecast_breach_day,
                truth = sf$true_breach, truth_day = sf$true_breach_day),
  intervention_day = iday,
  its_curve = its_curve,
  its = list(breakpoint = it$breakpoint,
             R_before = round2(it$R_before, 2), R_after = round2(it$R_after, 2),
             true_R_before = round2(sv$true_R_before, 2), true_R_after = round2(sv$true_R_after, 2),
             slowed = it$slowed, depletion_suspected = it$depletion_suspected)
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 3 -- sustained transmission, variant, next season ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

win3 <- 40:230
vs <- variant_selection(input, country, window = win3)
sq <- score_variant_selection(sim, vs, window = win3)
vf_truth <- truth_variant_freq(sim, country)
vlines <- vs$freq_line
vlines$truth <- vf_truth$variant_freq[match(vlines$day, vf_truth$day)]

# SIRS next-season scenarios
immune0 <- 1 - tail(sim$latent$local$susceptible[, country], 1) / pop
gi_mean <- epidist_mean(cfg$delays$generation_interval)
grid <- list(waning = c(slow = 1/365, med = 1/200, fast = 1/120),
             R0 = c(wildtype = 1.3, variant = 1.8), uptake = c(none = 0, campaign = 0.4))
res <- sirs_scenarios(grid, immune0 = immune0, gamma = 1/gi_mean, days = 300)
be <- boosting_effect(res, "variant")
# a few representative trajectories for a small-multiple
traj_of <- function(R0v, uptake) {
  S0 <- (1 - immune0) * (1 - grid$uptake[[uptake]])
  tr <- sirs_integrate(R0 = grid$R0[[R0v]], gamma = 1/gi_mean, omega = grid$waning[["med"]], S0 = S0, I0 = 1e-4, days = 300)
  data.frame(day = tr$day, I = tr$I)
}

phase3 <- list(
  variant = data.frame(day = vlines$day, observed = round2(vlines$observed, 3),
                       fitted = round2(vlines$fitted, 3), truth = round2(vlines$truth, 3)),
  selection = list(s = round2(vs$s, 3), lower = round2(vs$s_ci[1], 3), upper = round2(vs$s_ci[2], 3),
                   realized = round2(sq$s_realized, 3), crossover = round(vs$crossover_day)),
  scenarios = res[, c("scenario","waning","R0","uptake","peak_prevalence","cumulative_incidence","rel_peak")],
  boosting = list(peak_reduction = round2(be$peak_reduction * 100, 0),
                  incidence_reduction = round2(be$incidence_reduction * 100, 0)),
  immune0 = round2(immune0 * 100, 0),
  traj = list(wildtype_none = traj_of("wildtype","none"),
              variant_none  = traj_of("variant","none"),
              variant_campaign = traj_of("variant","campaign"))
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Assemble + write ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
out <- list(
  meta = list(source = src, country = country, n_days = cfg$n_days,
              start_date = format(cfg$start_date),
              variant_intro = format(cfg$variant$intro_date),
              ifr = cfg$ifr, gen_interval = round2(gi_mean, 1),
              total_infections = round(sum(sim$truth$infections$infections)),
              total_cases = round(sum(input$cases_by_report$cases)),
              total_deaths = round(sum(input$deaths$deaths_by_date))),
  phase0 = phase0, phase1 = phase1, phase2 = phase2, phase3 = phase3
)

writeLines(toJSON(out, auto_unbox = TRUE, digits = 6, na = "null", pretty = FALSE),
           "guidebook/guidebook_data.json")
cat("wrote guidebook/guidebook_data.json (", file.size("guidebook/guidebook_data.json"), "bytes )\n")
cat("meta:", out$meta$total_infections, "infections,", out$meta$total_cases, "cases,", out$meta$total_deaths, "deaths\n")
cat("phase0 growth R:", out$phase0$growth$R, "(truth", out$phase0$growth$truth_R, ")\n")
cat("phase1 nowcast improvement:", out$phase1$nowcast_improvement, "%\n")
cat("phase2 breach forecast day:", out$phase2$breach$forecast_day, "truth", out$phase2$breach$truth_day, "\n")
cat("phase3 selection s:", out$phase3$selection$s, "boosting peak reduction:", out$phase3$boosting$peak_reduction, "%\n")
