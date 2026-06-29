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

covid       <- c("2019/2020","2020/2021","2021/2022","2022/2023")
min_wk      <- 15
nonsent     <- c("MT","IS","HR","RO","LV","FI")
x1000       <- c("CY","LU","MT")
resp_only   <- c("NO","ES")
erviss_only <- c("SK","LV")

ip <- models_in$data_timeseries_long %>%
  filter(indicator=="ili_plus", pathogen=="Influenza", agegroup=="age_total") %>%
  select(stream, country_short, season, date, season_week, value)

erv <- ip %>% filter(stream %in% c("ili_plus_sentinel","ili_plus_nonsentinel")) %>%
  mutate(pick = ifelse(country_short %in% nonsent, "ili_plus_nonsentinel", "ili_plus_sentinel")) %>%
  filter(stream==pick) %>%
  mutate(value = value * ifelse(country_short %in% x1000, 1000, 1)) %>%
  transmute(country_short, season, date, season_week, erv = value)
rsp <- ip %>% filter(stream=="ili_plus_respicompass") %>%
  transmute(country_short, season, date, season_week, rsp = value)

m <- full_join(erv, rsp, by=c("country_short","season","date","season_week"))

# per-country alignment factor from the 2023/24 overlap (median RespiCompass / ERVISS over shared weeks)
fac <- m %>% filter(season=="2023/2024", is.finite(erv), erv>0, is.finite(rsp), rsp>0) %>%
  group_by(country_short) %>% summarise(factor = median(rsp/erv), .groups="drop")

m <- m %>% left_join(fac, by="country_short") %>%
  mutate(factor = ifelse(is.na(factor), 1, factor),
         value  = case_when(country_short %in% resp_only   ~ rsp,
                            country_short %in% erviss_only ~ erv,            # native ERVISS scale
                            !is.na(rsp) ~ rsp,                               # default: RespiCompass where present
                            TRUE        ~ erv * factor),                     # default ERVISS era, aligned
         source = case_when(country_short %in% resp_only   ~ "RespiCompass",
                            country_short %in% erviss_only ~ "ERVISS",
                            !is.na(rsp) ~ "RespiCompass", TRUE ~ "ERVISS")) %>%
  filter(!season %in% covid, is.finite(value))

# keep country-seasons with >=15 observed weeks; lay each on a contiguous weekly grid from week 1
slim <- m %>%
  group_by(country_short, season) %>%
  filter(sum(is.finite(value) & value>0) >= min_wk) %>%
  group_modify(function(df, key){
    src  <- df$source[which(!is.na(df$source))][1]
    w_hi <- max(df$season_week[is.finite(df$value)])
    tibble(season_week = 1:w_hi) %>%
      left_join(df %>% select(season_week, date, value), by="season_week") %>%
      arrange(season_week) %>%
      transmute(week = season_week, season_week, date, value, source = src)
  }) %>% ungroup() %>%
  arrange(country_short, season, week)

write.csv(slim, "data/slim_flu_iliplus.csv", row.names=FALSE)
cat(sprintf("wrote data/slim_flu_iliplus.csv: %d rows | %d countries | %d seasons | %d country-seasons\n",
            nrow(slim), n_distinct(slim$country_short), n_distinct(slim$season),
            nrow(distinct(slim, country_short, season))))
slim %>% distinct(country_short, season, source) %>% count(season, source) %>%
  tidyr::pivot_wider(names_from=source, values_from=n, values_fill=0) %>% as.data.frame() %>% print(row.names=FALSE)
