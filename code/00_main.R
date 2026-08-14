# 00_main.R -- build the data layer: raw committed snapshots -> model-ready inputs.
#
# RUN ORDER for the whole project (each later stage reads the earlier stage's outputs):
#   1. this script                                  -> output/models_in.rds (+ eyeballing manifest)
#   2. code/04_modelling/build_slim_panel.R         -> data/slim_flu_iliplus.csv (committed panel)
#   3. code/05_analysis/prepare_descriptors.R       -> output/descriptors{,_vax}.csv
#      code/05_analysis/dominant_subtype.R          -> output/dominant_subtype*.csv, descriptors_subtype.csv
#   4. the analyses/figures, each standalone from those tables:
#      code/03_report/{data,driver}_availability.R, code/04_modelling/{fit_methods_demo,
#      descriptive_overview,ekf_overview}.R, code/05_analysis/{analyse_patterns,plot_patterns,
#      plot_vax_scatter,hierarchical_models,bayes_subtype,subtype_8season,bayes_prior_burden,
#      bayes_precovid_ve_subtype,precovid_predict_postcovid}.R
#   Tests: Rscript run_tests.R (offline; rebuilds this pipeline from the committed snapshots).

# ---- |-Clear ----
gc() # free memory before the build (does NOT clear the environment)

# ---- |-Set up ----
source("code/01_main_supporting/setup.R")

# ---- |-load task specific settings ----
source("code/02_settings/settings_version0.R"); params=settings() # settings_version_X.R script to be changed by high-level user

# ---- |-sourcing support scripts ----
source("code/01_main_supporting/flu_functions.R")
source("code/01_main_supporting/validate.R")
source("code/01_main_supporting/load_data.R")
source("code/01_main_supporting/gen_model_input.R")
source("code/01_main_supporting/eyeballing.R")
source("code/01_main_supporting/sir_core.R")                          # shared SIR engine + data loaders
source("code/01_main_supporting/methods/method_sir_deterministic.R")  # method: deterministic SIR fit
source("code/01_main_supporting/methods/method_sir_ekf.R")            # method: EKF SIR (process noise)
source("code/01_main_supporting/methods/method_descriptive.R")        # method: descriptive curve features
source("code/01_main_supporting/methods_registry.R")                  # swappable-method registry + summaries

# ---- |-load flu data ----
data = load_data( params, regenerate = F, new_from_online = F) # loads the data # regenerate=T recreates the data lists, new_from_online=T uses the online versions for recreation
validate_data(data)          # fail loudly here, not deep inside gen_model_input (see validate.R)

# ---- |-generate model inputs ----
models_in = gen_model_input( params, data )
validate_models_in(models_in)
saveRDS(models_in, "output/models_in.rds") # persist the model-ready inputs (gitignored cache); the panel build + analysis read this (code/04_modelling/build_slim_panel.R, code/05_analysis/)

# ---- |-report ----
# render the data-eyeballing report: rmarkdown::render("code/03_report/eyeballing_report.Rmd")
# email it (manual utility, off by default): source("code/01_main_supporting/send_report.R"); send_report(params)

# ---- |-The end
eb = eyeballing(models_in, params, data) # build the quality + dynamics figure manifest for the report (see eyeballing.R)
