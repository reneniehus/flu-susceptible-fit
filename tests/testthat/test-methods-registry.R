# Tests the method registry + the common per-season summary schema (methods_registry.R). This is
# method-agnostic: it runs EVERY registered method through run_method() / summarise_method_fit()
# and checks they all honour the contract, so newly added methods are covered automatically.

sl <- load_flu_iliplus_slim("DK", path = here::here("data/slim_flu_iliplus.csv"))

test_that("every registered method runs and produces the standard summary schema", {
  schema <- c("method", "country", "season", "S0", "R_eff", "c", "auc", "peak_height",
              "peak_week", "onset_week", "steepness", "cor", "process_noise", "convergence")
  for (m in names(sir_methods())){
    fit  <- run_method(m, sl, params, n_starts = 1)   # contract/schema only -> cheapest fit
    expect_equal(fit$country, "DK")
    expect_equal(fit$method_name, m)

    smry <- summarise_method_fit(fit)
    expect_named(smry, schema)
    expect_equal(nrow(smry), length(sl$seasons))               # one row per season
    expect_true(all(smry$method == m))
    expect_true(all(is.finite(smry$S0) & smry$S0 > 0))         # not bounded < 1 (descriptive S0 is implied)
    expect_true(all(smry$R_eff > 0))
    expect_true(all(smry$auc > 0 & smry$peak_height > 0))
    expect_true(all(smry$peak_week >= 1, na.rm = TRUE))
    expect_true(all(is.finite(smry$onset_week)))               # threshold onset resolves for real waves
  }
})
