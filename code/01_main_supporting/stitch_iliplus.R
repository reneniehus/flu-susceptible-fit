# stitch_iliplus.R
#
# THE single implementation of the combined ILI+ panel-assembly rules (see
# documentation/decisions.md). Before this file existed the same stitch was copy-pasted in
# build_slim_panel.R, precovid_predict_postcovid.R and (approximately) data_availability.R,
# and the copies had started to diverge -- every consumer now calls stitch_iliplus_panel().
#
# The rules it encodes:
#   - ERVISS ILI+ = ILI rate x positivity, SENTINEL except NON-sentinel for MT/IS/HR/RO/LV/FI;
#     per-100-consultations countries CY/LU/MT scaled x1000 onto the per-100 000 basis;
#   - stitch PER WEEK: RespiCompass where present, else ERVISS aligned to the RespiCompass scale
#     by a per-country factor (median RespiCompass/ERVISS over the 2023/24 overlap weeks). In
#     practice RespiCompass covers <=2023/24 and ERVISS 2024/25+, but within the 2023/24 overlap
#     season aligned-ERVISS weeks DO fill RespiCompass gaps (14 country-seasons mix sources);
#     the per-season `source` label is the FIRST contributing source = the dominant early-season
#     one (RespiCompass for all current mixed seasons);
#   - single-source overrides: NO,ES = RespiCompass only; SK,LV = ERVISS only (native scale);
#   - optionally exclude the four COVID seasons; require >= min_wk strictly-POSITIVE observed
#     weeks (a season of zeros carries no wave); contiguous weekly grid from season week 1.
#
# Requires the tidyverse (source setup.R first). models_in comes from output/models_in.rds.

stitch_covid_seasons <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")

stitch_iliplus_panel <- function(models_in, exclude_covid = TRUE, min_wk = 15,
                                 overlap_season = "2023/2024"){
  nonsentinel <- c("MT", "IS", "HR", "RO", "LV", "FI")  # ERVISS ILI+ uses NON-sentinel positivity
  per_1000    <- c("CY", "LU", "MT")                    # per-100-consultations -> x1000 to per-100 000
  resp_only   <- c("NO", "ES")                          # single source = RespiCompass
  erviss_only <- c("SK", "LV")                          # single source = ERVISS

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
                              !is.na(respicompass) ~ "RespiCompass", TRUE ~ "ERVISS"))
  if (exclude_covid) combined <- combined %>% filter(!season %in% stitch_covid_seasons)
  combined <- combined %>% filter(is.finite(value))

  # keep country-seasons with >= min_wk positive observed weeks; contiguous weekly grid from week 1
  combined %>%
    group_by(country_short, season) %>%
    filter(sum(is.finite(value) & value>0) >= min_wk) %>%
    group_modify(function(df, key){
      season_source <- df$source[which(!is.na(df$source))][1]   # first contributing source = the dominant early-season one (overlap seasons can mix sources; see header)
      last_week     <- max(df$season_week[is.finite(df$value)])
      tibble(season_week = 1:last_week) %>%
        left_join(df %>% select(season_week, date, value), by="season_week") %>%
        arrange(season_week) %>%
        transmute(week = season_week, season_week, date, value, source = season_source)
    }) %>% ungroup() %>%
    arrange(country_short, season, week)
}
