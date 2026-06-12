# Contract tests: shields against upstream ERVISS column/file changes.

test_that("all expected epi streams load and pass the data contract", {
  expect_no_error(validate_data(data))
})

test_that("detailed typing streams carry the pathogen columns the extractors need", {
  for (nm in c("erviss_typing_sentinel", "erviss_typing_nonsentinel")) {
    expect_no_error(
      check_required_columns(data$epi[[nm]], c("pathogen", "pathogentype", "pathogensubtype",
                                               "indicator", "value", "date", "country_short"), nm)
    )
  }
})

test_that("the three respiratory pathogens are present in the typing data", {
  pathogens <- sort(unique(data$epi$erviss_typing_sentinel$pathogen))
  expect_true(all(c("Influenza", "SARS-CoV-2", "RSV") %in% pathogens))
})
