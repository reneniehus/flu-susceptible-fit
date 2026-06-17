# Shared setup, run once before the test files. Builds the data lists and models_in from the
# LOCAL snapshots (no internet needed), so the tests exercise the real pipeline offline.
library(here)
suppressMessages(source(here("code/01_main_supporting/setup.R")))
source(here("code/02_settings/settings_version0.R")); params <- settings()
for (f in c("flu_functions", "validate", "load_data", "gen_model_input", "eyeballing",
            "sir_core", "methods/method_sir_deterministic", "methods/method_sir_ekf",
            "methods/method_descriptive", "methods_registry"))
  source(here(paste0("code/01_main_supporting/", f, ".R")))

data      <- load_data(params, regenerate = FALSE, new_from_online = FALSE)
data      <- load_data_epi(data, params, regenerate = TRUE, new_from_online = FALSE) # fresh epi from snapshot
models_in <- gen_model_input(params, data)
