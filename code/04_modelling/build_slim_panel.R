# build_slim_panel.R
#
# Build the committed combined ILI+ panel (data/slim_flu_iliplus.csv) from models_in, encoding every
# data decision (see documentation/decisions.md):
#   - ERVISS ILI+ = ILI rate x positivity, SENTINEL except NON-sentinel for MT/IS/HR/RO/LV/FI;
#     per-100-consultations countries CY/LU/MT scaled x1000 onto the per-100 000 basis (units);
#   - source per country-season: RespiCompass <=2023/24, ERVISS 2024/25+ (default stitch), with the
#     ERVISS era aligned to the RespiCompass scale by a per-country factor estimated from 2023/24;
#   - single-source overrides: NO,ES = RespiCompass only; SK,LV = ERVISS only (native scale);
#   - exclude the four COVID seasons; require >=15 observed weeks; contiguous weekly grid from the
#     season start (week 1).
#
# Run from the repo root:  Rscript code/04_modelling/build_slim_panel.R

source("code/01_main_supporting/setup.R")
models_in <- readRDS("output/models_in.rds")

covid          <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")
min_wk         <- 15
overlap_season <- "2023/2024"                     # the season RespiCompass and ERVISS overlap, used to align their scales
nonsentinel    <- c("MT", "IS", "HR", "RO", "LV", "FI")  # countries whose ERVISS ILI+ uses NON-sentinel positivity
per_1000       <- c("CY", "LU", "MT")             # ERVISS reports these per-100-consultations -> x1000 onto the per-100 000 basis
resp_only      <- c("NO", "ES")                   # single source = RespiCompass
erviss_only    <- c("SK", "LV")                   # single source = ERVISS

ili_plus <- models_in$data_timeseries_long %>%
  filter(indicator=="ili_plus", pathogen=="Influenza", agegroup=="age_total") %>%
  select(stream, country_short, season, date, season_week, value)

erviss <- ili_plus %>% filter(stream %in% c("ili_plus_sentinel","ili_plus_nonsentinel")) %>%
  mutate(chosen_stream = ifelse(country_short %in% nonsentinel, "ili_plus_nonsentinel", "ili_plus_sentinel")) %>%
  filter(stream==chosen_stream) %>%
  mutate(value = value * ifelse(country_short %in% per_1000, 1000, 1)) %>%
  transmute(country_short, season, date, season_week, erviss = value)
respicompass <- ili_plus %>% filter(stream=="ili_plus_respicompass") %>%
  transmute(country_short, season, date, season_week, respicompass = value)

combined <- full_join(erviss, respicompass, by=c("country_short","season","date","season_week"))

# per-country alignment factor from the overlap season (median RespiCompass / ERVISS over weeks where BOTH streams are finite and > 0)
align_factors <- combined %>% filter(season==overlap_season, is.finite(erviss), erviss>0, is.finite(respicompass), respicompass>0) %>%
  group_by(country_short) %>% summarise(align_factor = median(respicompass/erviss), .groups="drop")

combined <- combined %>% left_join(align_factors, by="country_short") %>%
  mutate(align_factor = ifelse(is.na(align_factor), 1, align_factor),
         value  = case_when(country_short %in% resp_only   ~ respicompass,
                            country_short %in% erviss_only ~ erviss,               # native ERVISS scale
                            !is.na(respicompass) ~ respicompass,                   # default: RespiCompass where present
                            TRUE        ~ erviss * align_factor),                  # default ERVISS era, aligned to RespiCompass
         source = case_when(country_short %in% resp_only   ~ "RespiCompass",
                            country_short %in% erviss_only ~ "ERVISS",
                            !is.na(respicompass) ~ "RespiCompass", TRUE ~ "ERVISS")) %>%
  filter(!season %in% covid, is.finite(value))

# keep country-seasons with >=15 observed weeks; lay each on a contiguous weekly grid from week 1
slim <- combined %>%
  group_by(country_short, season) %>%
  filter(sum(is.finite(value) & value>0) >= min_wk) %>%
  group_modify(function(df, key){
    season_source <- df$source[which(!is.na(df$source))][1]   # source is constant within a country-season; take the first non-NA
    last_week     <- max(df$season_week[is.finite(df$value)])
    tibble(season_week = 1:last_week) %>%
      left_join(df %>% select(season_week, date, value), by="season_week") %>%
      arrange(season_week) %>%
      transmute(week = season_week, season_week, date, value, source = season_source)
  }) %>% ungroup() %>%
  arrange(country_short, season, week)

write.csv(slim, "data/slim_flu_iliplus.csv", row.names=FALSE)
cat(sprintf("wrote data/slim_flu_iliplus.csv: %d rows | %d countries | %d seasons | %d country-seasons\n",
            nrow(slim), n_distinct(slim$country_short), n_distinct(slim$season),
            nrow(distinct(slim, country_short, season))))
slim %>% distinct(country_short, season, source) %>% count(season, source) %>%
  tidyr::pivot_wider(names_from=source, values_from=n, values_fill=0) %>% as.data.frame() %>% print(row.names=FALSE)
