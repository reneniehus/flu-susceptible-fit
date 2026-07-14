# Shared test setup: source the whole playground and build one small, fast simulation the analysis
# tests reuse. No external data -- the simulator IS the fixture.
library(testthat)
PP_ROOT <- Sys.getenv("PP_ROOT", unset = normalizePath(file.path("..", "..")))

engine <- c("epidist", "config", "utils", "renewal", "draw_parameters", "simulate_source",
            "simulate_flights", "simulate_importations", "simulate_local", "simulate_clusters",
            "observe", "assemble")
for (f in engine) source(file.path(PP_ROOT, "R", paste0(f, ".R")))

toolbox <- c("analysis_common", "phase0_growth_R", "phase0_clusters", "phase0_catchment",
             "phase0_importation_risk", "phase1_cfr", "phase1_rt", "phase1_nowcast",
             "phase2_forecast", "phase2_intervention", "phase3_scenarios", "phase3_variant_selection")
for (f in toolbox) source(file.path(PP_ROOT, "analysis", paste0(f, ".R")))

# a small, fast fixture: 10 countries (a deliberate mix of high- and low-surveillance ones, so the
# importation-risk regression has anchor countries AND under-detectors) over 160 days -- enough for
# growth, a peak, deaths, capacity breaches and a rising variant.
test_cfg   <- config_subset(default_config(),
                            c("IT", "BE", "ES", "FR", "CZ", "DE", "NL", "PL", "DK", "RO"), n_days = 160)
test_sim   <- simulate_pandemic(test_cfg)
test_input <- as_analysis_input(test_sim)
