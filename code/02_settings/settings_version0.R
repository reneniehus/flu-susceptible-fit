settings = function() {
  params = list()
  
  # ---- |-Run modes ----
  params$save_submission = F # T: saves the file ready for respicompass, F; will be faster
  
  # debug/fast modes
  params$rapid_stan_fit = T # T: runs scripts with settings that reduce run-time
  
  # ---- |-Report email settings (used by send_report(); off by default) ----
  params$send_report       = FALSE                       # TRUE: email the rendered report
  params$report_from       = "you@example.org"           # sender address
  params$report_recipients = c("you@example.org")        # one or more recipients
  params$report_subject    = "Model run complete"
  params$report_attachments= c("code/03_report/eyeballing_report.html") # files to attach (missing ones skipped)
  params$smtp_host         = ""                           # your SMTP server, e.g. "smtp.example.org"
  params$smtp_port         = 25
  params$smtp_insecure     = TRUE
  
  # ---- |-Names/identifiers ----
  params$four_age_groups = c("0-4","5-14","15-64","65+") # the order is important
  
  # ---- |-Disease parameters ----
  params$Rnull = 1.5 #

  # immunity parameters
  params$ve_spread = 0.20 # vaccine effect on onward spread when vaccinated individual is infected

  # ---- |-Susceptibility fits (methods in code/01_main_supporting/methods/) ----
  # R0 and the seed are FIXED here; the per-season susceptibility S0 and reporting fraction c are
  # fitted. Only the RELATIVE S0 across seasons is interpreted (absolute S0 is conditional on these).
  # Shared by every method (deterministic, EKF, ...).
  params$susc_R0                    = 1.5   # fixed seasonal-influenza R0 (literature)
  params$susc_infectious_period_days = 3    # mean infectious period -> gamma = 7/this (per week)
  params$susc_seed_i0               = 1e-5  # constant seed (~0.001% of pop, southern-hemisphere import)
  params$susc_countries             = c("DK", "FR", "IE", "HU")  # countries for the SIR / method fits + demos (NOT the slim panel, which spans all 25 countries with usable data; see build_slim_panel.R)
  params$susc_smooth_window         = 4     # centered moving-average window (weeks) for the descriptive method
  
  # ---- |-Data ----
  params$latest_start_year = 2025 # if the last partly/fully observed season is 2024/25, put 2024
  params$season_start_monthday = "-08-01" # initial date of for SIR initiation
  params$season_end_monthday = "-07-31" # end date of SIR process

  # data-loading settings (consumed by code/01_main_supporting/load_data.R)
  params$respicompass_round = "2024-2025_round_1" # RespiCompass hub round folder; bump this for a new season
  params$demography_year = 2024 # ReportYear used when querying the ECDC population database
  params$use_ecdc_db = FALSE # keep FALSE: demography comes from the committed snapshot (the ECDC SQL client is not part of this repo)

  # ---- |-Simulations ----
  params$simulation_seed = 12
  
  # ---- |-Countries ----
  params$run_countries = c("IT", "AT", "BE", "BG", "HR")
  
  # ---- |-Model-specific  settings ----
  
  # ---- |-Fitting and uncertainty ----
  
  # ---- |-Flu scenarios ----
  
  # ---- |-Folder paths ----
  
  return(params)
}

 



