# Invariants for the ILI+ construction.

test_that("ILI+ equals ILI consultation rate x positivity (exact internal invariant)", {
  L   <- models_in$data_timeseries_long
  ili <- L %>% filter(stream == "ili_ari", indicator == "ILIconsultationrate") %>%
    select(country_short, date, agegroup, ili = value)
  pos <- L %>% filter(stream == "typing_sentinel", indicator == "positivity", pathogen == "Influenza") %>%
    select(country_short, date, pos = value)
  ip  <- L %>% filter(stream == "ili_plus_sentinel", pathogen == "Influenza") %>%
    select(country_short, date, agegroup, iliplus = value)

  recon <- ili %>%
    inner_join(pos, by = c("country_short", "date")) %>%
    inner_join(ip,  by = c("country_short", "date", "agegroup")) %>%
    mutate(diff = abs(iliplus - ili * pos))

  expect_gt(nrow(recon), 1000)
  expect_lt(max(recon$diff, na.rm = TRUE), 1e-9)
})

test_that("influenza ILI+ reproduces the RespiCompass ILI+ (median country correlation ~1)", {
  cmp <- models_in$data_timeseries_long %>%
    filter(indicator == "ili_plus", agegroup == "age_total",
           (stream == "ili_plus_sentinel" & pathogen == "Influenza") | stream == "ili_plus_respicompass") %>%
    select(country_short, date, stream, value) %>%
    pivot_wider(names_from = stream, values_from = value) %>%
    filter(!is.na(ili_plus_sentinel), !is.na(ili_plus_respicompass))

  per_country <- cmp %>%
    group_by(country_short) %>% filter(n() > 20) %>%
    summarise(cor = cor(ili_plus_sentinel, ili_plus_respicompass), .groups = "drop")

  expect_gt(nrow(per_country), 10)
  expect_gt(median(per_country$cor), 0.99)   # robust to the few countries using different positivity
})
