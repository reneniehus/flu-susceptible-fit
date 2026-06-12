# Tests for the canonical tables produced by gen_model_input().

test_that("models_in exposes the three canonical tables and passes its contract", {
  expect_true(all(c("data_timeseries_long", "data_timeseries_wide", "data_season_summary") %in% names(models_in)))
  expect_no_error(validate_models_in(models_in))
})

test_that("ILI+ is constructed for influenza, SARS-CoV-2 and RSV", {
  pathogens <- models_in$data_timeseries_long %>%
    filter(indicator == "ili_plus") %>% pull(pathogen) %>% unique()
  expect_true(all(c("Influenza", "SARS-CoV-2", "RSV") %in% pathogens))
})

test_that("season summary carries both completeness measures, bounded in [0, 1]", {
  s <- models_in$data_season_summary
  expect_true(all(c("completeness", "completeness_active") %in% names(s)))
  expect_true(all(s$completeness        >= 0 & s$completeness        <= 1, na.rm = TRUE))
  expect_true(all(s$completeness_active >= 0 & s$completeness_active <= 1, na.rm = TRUE))
})

test_that("the wide table is keyed one row per country x season x date x agegroup", {
  w <- models_in$data_timeseries_wide
  expect_true(all(c("country_short", "season", "date", "agegroup") %in% names(w)))
  expect_equal(nrow(distinct(w, country_short, season, date, agegroup)), nrow(w))
})
