# Contract tests for the committed slim ILI+ panel (data/slim_flu_iliplus.csv) -- the offline
# input every susceptibility method fits. The panel encodes the data decisions recorded in
# documentation/decisions.md (COVID-season exclusion, >=15-observed-week inclusion rule,
# single-source countries); these tests pin them directly on the committed file, so a rebuild that
# silently drops a rule fails here even if the modelling code still runs. The last test re-runs the
# build pipeline in memory and requires it to reproduce the committed CSV exactly, guarding against
# the cached models_in, the build script and the committed file drifting apart.

slim_csv <- read.csv(here::here("data/slim_flu_iliplus.csv"), stringsAsFactors = FALSE)

test_that("the four COVID-disrupted seasons are excluded, and both eras are represented", {
  # decisions.md: 2019/2020-2022/2023 are treated as disrupted by the pandemic (NPIs, collapsed
  # influenza circulation, changed care-seeking/testing) and must never enter the fits.
  covid <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")
  expect_false(any(slim_csv$season %in% covid))
  # the whole point of the two-source stitch is a long record spanning BOTH eras
  expect_true(any(slim_csv$season <= "2018/2019"))   # pre-COVID (RespiCompass era)
  expect_true(any(slim_csv$season >= "2023/2024"))   # post-COVID (ERVISS era)
})

test_that("every country-season has >= 15 weeks with finite, strictly positive ILI+", {
  # The inclusion rule counts STRICTLY positive weeks (a season of zeros carries no wave to fit);
  # country-seasons below 15 such weeks are dropped at build time.
  wk <- slim_csv %>%
    group_by(country_short, season) %>%
    summarise(n_pos = sum(is.finite(value) & value > 0), .groups = "drop")
  expect_gte(nrow(wk), 100)                          # the panel spans many country-seasons (25 countries)
  expect_true(all(wk$n_pos >= 15))
})

test_that("single-source overrides hold: NO/ES are RespiCompass-only, SK/LV ERVISS-only", {
  # decisions.md: these countries have no clean 2023/24 overlap (or an anomalous series, LV), so
  # each is kept on ONE source and the other is dismissed -- mixing them back in would reintroduce
  # exactly the source/era artefact the stitch is designed to avoid.
  for (cc in c("NO", "ES"))
    expect_equal(unique(slim_csv$source[slim_csv$country_short == cc]), "RespiCompass")
  for (cc in c("SK", "LV"))
    expect_equal(unique(slim_csv$source[slim_csv$country_short == cc]), "ERVISS")
})

test_that("the panel rebuilds identically from the cached models_in", {
  # Replicates code/04_modelling/build_slim_panel.R in memory (no file is written) from the same
  # input the script reads (output/models_in.rds), and requires the result to equal the committed
  # CSV. If the build script's rules and this replica ever diverge, or the committed file goes
  # stale against the cache, this fails. KEEP IN SYNC with build_slim_panel.R.
  skip_if_not(file.exists(here::here("output/models_in.rds")),
              "output/models_in.rds not built (run code/00_main.R first)")
  mi <- readRDS(here::here("output/models_in.rds"))

  covid          <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")
  min_wk         <- 15
  overlap_season <- "2023/2024"                     # season where RespiCompass and ERVISS overlap
  nonsentinel    <- c("MT", "IS", "HR", "RO", "LV", "FI")  # ERVISS ILI+ from NON-sentinel positivity
  per_1000       <- c("CY", "LU", "MT")             # per-100-consultations -> x1000 onto per-100 000
  resp_only      <- c("NO", "ES")
  erviss_only    <- c("SK", "LV")

  ili_plus <- mi$data_timeseries_long %>%
    filter(indicator == "ili_plus", pathogen == "Influenza", agegroup == "age_total") %>%
    select(stream, country_short, season, date, season_week, value)

  erviss <- ili_plus %>% filter(stream %in% c("ili_plus_sentinel", "ili_plus_nonsentinel")) %>%
    mutate(chosen_stream = ifelse(country_short %in% nonsentinel, "ili_plus_nonsentinel", "ili_plus_sentinel")) %>%
    filter(stream == chosen_stream) %>%
    mutate(value = value * ifelse(country_short %in% per_1000, 1000, 1)) %>%
    transmute(country_short, season, date, season_week, erviss = value)
  respicompass <- ili_plus %>% filter(stream == "ili_plus_respicompass") %>%
    transmute(country_short, season, date, season_week, respicompass = value)

  combined <- full_join(erviss, respicompass, by = c("country_short", "season", "date", "season_week"))

  align_factors <- combined %>%
    filter(season == overlap_season, is.finite(erviss), erviss > 0, is.finite(respicompass), respicompass > 0) %>%
    group_by(country_short) %>% summarise(align_factor = median(respicompass / erviss), .groups = "drop")

  combined <- combined %>% left_join(align_factors, by = "country_short") %>%
    mutate(align_factor = ifelse(is.na(align_factor), 1, align_factor),
           value  = case_when(country_short %in% resp_only   ~ respicompass,
                              country_short %in% erviss_only ~ erviss,            # native ERVISS scale
                              !is.na(respicompass) ~ respicompass,
                              TRUE ~ erviss * align_factor),
           source = case_when(country_short %in% resp_only   ~ "RespiCompass",
                              country_short %in% erviss_only ~ "ERVISS",
                              !is.na(respicompass) ~ "RespiCompass", TRUE ~ "ERVISS")) %>%
    filter(!season %in% covid, is.finite(value))

  rebuilt <- combined %>%
    group_by(country_short, season) %>%
    filter(sum(is.finite(value) & value > 0) >= min_wk) %>%
    group_modify(function(df, key){
      season_source <- df$source[which(!is.na(df$source))][1]
      last_week     <- max(df$season_week[is.finite(df$value)])
      tibble(season_week = 1:last_week) %>%
        left_join(df %>% select(season_week, date, value), by = "season_week") %>%
        arrange(season_week) %>%
        transmute(week = season_week, season_week, date, value, source = season_source)
    }) %>% ungroup() %>%
    arrange(country_short, season, week)

  # normalise both sides to plain data frames (CSV round-trip turns Date into character)
  norm <- function(d){
    d <- as.data.frame(d)
    d$date <- as.character(d$date)
    rownames(d) <- NULL
    d
  }
  expect_equal(norm(rebuilt), norm(slim_csv), tolerance = 1e-8)  # 1e-8 >> write.csv round-trip error
})
