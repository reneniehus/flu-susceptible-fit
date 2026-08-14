# Live settings only: every knob here has a consumer in code/ (grep it before adding back a
# removed one -- dead knobs invite silent drift, e.g. a season stamp once came from a knob that
# had drifted from the pinned data round). Historical knobs live in git history.
settings = function() {
  params = list()

  # ---- |-Report email settings (used by send_report(), a manual utility; off by default) ----
  params$send_report       = FALSE                       # TRUE: email the rendered report
  params$report_from       = "you@example.org"           # sender address
  params$report_recipients = c("you@example.org")        # one or more recipients
  params$report_subject    = "Model run complete"
  params$report_attachments= c("code/03_report/eyeballing_report.html") # files to attach (missing ones skipped)
  params$smtp_host         = ""                           # your SMTP server, e.g. "smtp.example.org"
  params$smtp_port         = 25
  params$smtp_insecure     = TRUE

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
  params$season_start_monthday = "-08-01" # season window start (Aug -> Jul); drives season labels + season_week
  params$season_end_monthday = "-07-31" # season window end

  # data-loading settings (consumed by code/01_main_supporting/load_data.R)
  params$respicompass_round = "2024-2025_round_1" # RespiCompass hub round folder; bump this for a new season (also stamps the vaccination-scenario season)
  params$use_ecdc_db = FALSE # keep FALSE: demography comes from the committed snapshot (the ECDC SQL client is not part of this repo)

  return(params)
}

 



