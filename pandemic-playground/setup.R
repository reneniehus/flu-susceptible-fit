# setup.R
#
# Source the whole playground in dependency order, from the project root:
#   source("setup.R")            # simulation engine + analysis toolbox
#   sim <- simulate_pandemic()   # a full synthetic pandemic (truth + observed)
#
# The engine (R/) has NO external dependencies -- base R only -- so it runs anywhere. The analysis
# toolbox (analysis/) uses only base R + the tidyverse for tabulation; a couple of tools optionally
# use {nnet}. Nothing needs the epidemiological packages the methods emulate ({epiparameter},
# EpiEstim, cfr, EpiNow2, ...), which are documented extension fronts (see documentation/decisions.md).

.pp_root <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA)
if (is.na(.pp_root) || !length(.pp_root)) .pp_root <- "."

# ---- |-simulation engine (base R; order matters -- later files use earlier ones) ----
.engine <- c("epidist", "config", "utils", "renewal", "draw_parameters",
             "simulate_source", "simulate_flights", "simulate_importations",
             "simulate_local", "observe", "assemble")
for (f in .engine) source(file.path(.pp_root, "R", paste0(f, ".R")))

# ---- |-analysis toolbox (base R + tidyverse; the phase files + shared scoring) ----
.toolbox <- c("analysis_common",
              "phase0_growth_R", "phase0_catchment", "phase0_importation_risk",
              "phase1_cfr", "phase1_rt", "phase1_nowcast",
              "phase2_forecast", "phase2_intervention",
              "phase3_scenarios", "phase3_variant_selection")
for (f in .toolbox) {
  path <- file.path(.pp_root, "analysis", paste0(f, ".R"))
  if (file.exists(path)) source(path)
}

rm(.engine, .toolbox, .pp_root)
