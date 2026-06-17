# ---- |-Clear ----
gc() # clear environment & memory

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
source("code/01_main_supporting/run_model.R")
source("code/01_main_supporting/sir_core.R")                          # shared SIR engine + data loaders
source("code/01_main_supporting/methods/method_sir_deterministic.R")  # method: deterministic SIR fit
source("code/01_main_supporting/methods/method_sir_ekf.R")            # method: EKF SIR (process noise)
source("code/01_main_supporting/methods/method_descriptive.R")        # method: descriptive curve features
source("code/01_main_supporting/methods_registry.R")                  # swappable-method registry + summaries
source("code/01_main_supporting/process_and_save.R")
source("code/01_main_supporting/send_report.R")

# ---- |-load flu data ----
data = load_data( params, regenerate = F, new_from_online = F) # loads the data # regenerate=T recreates the data lists, new_from_online=T uses the online versions for recreation

# ---- |-generate model inputs ----
models_in = gen_model_input( params, data )
# ---- |-run flu models----
models_out = run_model( params, data , models_in ) # runs the model scripts

# ---- |-process and save model output ----

# ---- |-report ----
# render the data-eyeballing report: rmarkdown::render("code/03_report/eyeballing_report.Rmd")

# ---- |-The end
# (temporary code for any quick checking)
eb = eyeballing(models_in, params, data) # quality + dynamics figure manifest (see eyeballing.R)
