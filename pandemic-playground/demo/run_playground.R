#!/usr/bin/env Rscript
# run_playground.R
#
# The guided tour: simulate one synthetic pandemic, then walk the analysis toolbox phase by phase,
# scoring each tool against the known truth. Run from the project root:
#   Rscript demo/run_playground.R
#
# This is the playground's whole promise in one script -- every analytical answer is put next to the
# truth it was trying to recover, so you can SEE how well (or badly) each method does, and why.

source("setup.R")

rule <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Simulate a pandemic ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule("SIMULATE  --  latent truth + observed surveillance from one config + seed")
sim   <- simulate_pandemic(default_config(), quiet = FALSE)
input <- as_analysis_input(sim)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 0 -- before local introduction (at the source X) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule("PHASE 0  --  before local introduction (source X still dominates)")

cat("\nQ: Is it spreading person-to-person, and how fast?\n")
ga <- growth_analysis(input, "X", window = 20:45)
sg <- score_growth(sim, ga)
cat(sprintf("   growth rate r = %.3f/day (doubling %.1f d)  ->  R = %.2f [%.2f, %.2f]\n",
            ga$r, ga$doubling_time, ga$R, ga$R_ci[1], ga$R_ci[2]))
cat(sprintf("   TRUTH: early R = %.2f   (on the OBSERVED case curve the scaling-up testing inflates R --\n",
            sg$truth_R))
cat(sprintf("          fit the true infection curve instead and the same tool returns R = %.2f)\n",
            { i <- truth_infections(sim, "X"); d <- i[i$day %in% 20:45, ]
              gp <- discretise(input$delays$generation_interval, boundary = "cori")
              round(r_to_R(estimate_growth_rate(d$infections, d$day)$r, gp), 2) }))

cat("\nQ: How much are we missing at source -- how big is it really?\n")
cb <- catchment_backcalc(input, 30:60, min_surveillance = 0.75, source_pop = sim$config$source$population)
sc <- score_catchment(sim, cb, input)
cat(sprintf("   back-calculated source prevalence = %s infectious   (TRUTH %s; ratio %.2f)\n",
            format(round(cb$est_prevalence), big.mark = ","), format(round(sc$true_prevalence), big.mark = ","), sc$ratio))
cat(sprintf("   meanwhile the source reported only %s cases in that window -- the true epidemic is %.0fx bigger\n",
            format(sc$reported_source_cases, big.mark = ","), sc$true_underascertainment))

cat("\nQ: Which of our countries are under-detecting their imports?\n")
ir <- importation_risk(input, window = 25:70, min_surveillance = 0.8)
si <- score_importation_risk(ir)
cat(sprintf("   flagged %d countries below the imports-vs-flights line;\n", si$n_flagged))
cat(sprintf("   their mean surveillance quality %.2f vs %.2f for the rest (cor(residual, true quality) = %.2f)\n",
            si$mean_surv_flagged, si$mean_surv_unflagged, si$cor_residual_surveillance))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 1 -- early local exponential growth (a worked country) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
country <- "IT"
rule(sprintf("PHASE 1  --  early local growth (%s)", country))

cat("\nQ: How deadly is it?  (delay-adjusted CFR)\n")
for (t in c(70, 100, 160)) {
  s  <- cfr_static(input, country, as_of = t)
  sk <- score_cfr(sim, s, input)
  cat(sprintf("   day %3d: naive CFR %.3f  vs  delay-adjusted %.3f [%.3f, %.3f]   (true confirmed CFR %.3f; true IFR %.3f)\n",
              t, s$cfr_naive, s$cfr_adjusted, s$cfr_lower, s$cfr_upper, sk$true_confirmed_cfr, sk$true_ifr))
}
cat("   -> the adjusted CFR is stable early; the gap to the IFR is the (uncorrected) case under-ascertainment\n")

cat("\nQ: Is transmission growing or shrinking right now?  (Cori Rt)\n")
rt <- rt_analysis(input, country)                     # window derived from the generation interval
sr <- score_rt(sim, rt)
cat(sprintf("   Rt vs the realized truth: cor %.2f, MAE %.2f, %.0f%% CI coverage\n",
            sr$cor, sr$mae, 100 * sr$coverage))

cat("\nQ: How many cases really occurred recently, before reports catch up?  (nowcast)\n")
nc <- nowcast_truncation(input, country, as_of = 100)
sn <- score_nowcast(sim, nc, recent = 14)
cat(sprintf("   last-14-day onset counts: nowcast RMSE %.0f vs %.0f if you just trust the raw truncated data (%.0f%% better)\n",
            sn$nowcast$rmse, sn$observed_naive$rmse, 100 * sn$rmse_improvement))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 2 -- growth to peak, healthcare demand ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule(sprintf("PHASE 2  --  growth to peak, healthcare demand (%s)", country))
cap <- truth_capacity(sim, country)

cat("\nQ: Will we breach hospital capacity in the next few weeks?  (renewal forecast)\n")
pop_c <- sim$par$country_pop[[country]]
fc <- forecast_capacity(input, country, as_of = 95, horizon = 28, capacity = cap, population = pop_c)
sf <- score_forecast(sim, fc, input)
cat(sprintf("   R_now %.2f, capacity %.0f/day. Forecast: breach = %s (day %s);  TRUTH: breach = %s (day %s)\n",
            fc$R_now, cap, sf$forecast_breach, sf$forecast_breach_day, sf$true_breach, sf$true_breach_day))

cat("\nQ: Are our control measures working?  (interrupted time series)\n")
iday <- sim$config$rt_country[[country]]$t[2]
it <- intervention_its(input, country, intervention_day = iday, window = 21)
sv <- score_intervention(sim, it, window = 21)
cat(sprintf("   at the intervention: R %.2f -> %.2f (growth %s);  TRUTH: realized R %.2f -> %.2f  [association, not proof]\n",
            it$R_before, it$R_after, if (it$slowed) "significantly slowed" else "not clearly slowed",
            sv$true_R_before, sv$true_R_after))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PHASE 3 -- variant takeover and next-season scenarios ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule("PHASE 3  --  sustained transmission, a fitter variant, next season")

cat("\nQ: Is the new variant taking over, and how fast?  (selection coefficient)\n")
vs <- variant_selection(input, country, window = 40:200)
sq <- score_variant_selection(sim, vs, window = 40:200)
cat(sprintf("   selection coefficient s = %.3f/day [%.3f, %.3f]  (realized truth %.3f);  variant passes 50%% on day %.0f\n",
            vs$s, vs$s_ci[1], vs$s_ci[2], sq$s_realized, vs$crossover_day))

cat("\nQ: What might next season look like, and does boosting help?  (SIRS scenarios)\n")
immune0 <- 1 - tail(sim$latent$local$susceptible[, country], 1) / sim$par$country_pop[[country]]
gi_mean <- epidist_mean(sim$config$delays$generation_interval)           # match the SIRS speed to the pathogen
grid <- list(waning = c(slow = 1/365, med = 1/200, fast = 1/120),
             R0 = c(wildtype = 1.3, variant = 1.8), uptake = c(none = 0, campaign = 0.4))
res <- sirs_scenarios(grid, immune0 = immune0, gamma = 1 / gi_mean, days = 365)
be  <- boosting_effect(res, "variant")
cat(sprintf("   seeding next season at %.0f%% immune. A fitter variant raises the peak ~%.1fx vs wild type (no booster).\n",
            100 * immune0, mean(res$rel_peak[res$R0 == "variant" & res$uptake == "none"])))
cat(sprintf("   A 40%% booster campaign cuts the variant peak by ~%.0f%% and the cumulative incidence by ~%.0f%% (ensembled over waning).\n",
            100 * be$peak_reduction, 100 * be$incidence_reduction))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Beyond COVID -- the framework represents other respiratory pathogens ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule("BEYOND COVID  --  the same framework, a superspreading (SARS-like) pathogen")
cat("\nSuperspreading: with the offspring dispersion k, most infections transmit to no-one and imported\n")
cat("chains often fizzle -- the defining SARS/MERS behaviour a Poisson process cannot produce.\n")
gi_pmf <- discretise(sim$config$delays$generation_interval, boundary = "cori")
establish <- function(k, reps = 300) {
  mean(vapply(seq_len(reps), function(i) {
    set.seed(i); seed <- matrix(0, 120, 1); seed[1, 1] <- 1
    sum(simulate_renewal(120, 1e7, gi_pmf, list(step_schedule(0, 2.4)), seed, dispersion = k)$incidence) > 50
  }, logical(1)))
}
cat(sprintf("   establishment probability from one imported case (R0 = 2.4):  Poisson %.0f%%   vs   k=0.16 (SARS) %.0f%%\n",
            100 * establish(Inf), 100 * establish(0.16)))

rule("DONE  --  every answer above was scored against a truth the analyst never saw")
