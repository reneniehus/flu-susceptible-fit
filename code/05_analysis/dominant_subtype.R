# dominant_subtype.R
#
# Determine the DOMINANT influenza (sub)type per country-season from ERVISS typing, HIERARCHICALLY:
#   1. TYPE level: plurality of ALL characterised A vs ALL characterised B detections -- including
#      the type-known-but-unsubtyped rows ('A (unknown)', 'B (unknown)'), which are the majority of
#      what labs report. (An earlier rule counted 'B (unknown)' toward B but DISCARDED 'A (unknown)';
#      in 2024/25 unsubtyped A was ~206k of ~260k A detections EU/EEA-wide, so that asymmetry called
#      'B' in 22/30 countries while type-level A actually exceeded B in 28/30 -- an artifact of the
#      counting rule, not the data. See documentation/decisions.md.)
#   2. SUBTYPE level: if type A wins, the dominant subtype is the plurality among SUBTYPED A
#      (A(H1)pdm09 vs A(H3)); a country-season with type A dominant but zero subtyped A keeps
#      dominant_type = "A" and dominant = NA (unassignable at subtype level, left out of subtype
#      analyses). If type B wins, dominant = "B" (lineages recorded separately by ERVISS).
# We combine sentinel + non-sentinel detections and require >= MIN_TYPED_DETECTIONS characterised
# (A+B) detections for a call; top_share (type level) and subtype_share (within subtyped A) keep
# "clear" vs "mixed/co-circulating" transparent.
#
# Subtype data is ERVISS typing, which starts 2021 -> only the post-COVID analysis seasons get a
# subtype; the pre-COVID (RespiCompass) seasons have none. The 'EU/EEA' aggregate rows are NOT a
# country: they are split off into a separate continental table (printed for cross-checking
# data/external/dominant_subtype_by_season.csv). Writes output/dominant_subtype.csv and
# output/descriptors_subtype.csv (descriptors + dominant subtype). Run from the repo root.

source("code/01_main_supporting/setup.R")

MIN_TYPED_DETECTIONS <- 20   # a country-season needs at least this many characterised (A+B) detections

typing <- read.csv("data/erviss_flu_type_subtype.csv", stringsAsFactors=FALSE) %>%
  filter(pathogen=="Influenza", indicator=="detections") %>%
  mutate(date        = ISOweek::ISOweek2date(paste0(yearweek, "-3")),
         season_year = year(date) - (month(date) < 8), season = paste0(season_year, "/", season_year+1),
         type_group  = case_when(pathogentype == "Influenza A" ~ "A",     # includes 'A (unknown)'
                                 pathogentype == "Influenza B" ~ "B",     # includes 'B (unknown)'
                                 TRUE                          ~ NA_character_),   # 'Influenza untyped' -> dropped
         subtype_group = case_when(pathogensubtype == "A(H1)pdm09" ~ "A(H1N1)",
                                   pathogensubtype == "A(H3)"      ~ "A(H3N2)",
                                   TRUE                            ~ NA_character_))

# ---- |-hierarchical dominance for one set of typing rows (already grouped by country/season) ----
dominance_call <- function(df){
  n_A <- sum(df$value[df$type_group=="A"], na.rm=TRUE)
  n_B <- sum(df$value[df$type_group=="B"], na.rm=TRUE)
  n_h1 <- sum(df$value[!is.na(df$subtype_group) & df$subtype_group=="A(H1N1)"], na.rm=TRUE)
  n_h3 <- sum(df$value[!is.na(df$subtype_group) & df$subtype_group=="A(H3N2)"], na.rm=TRUE)
  dominant_type <- if (n_A >= n_B) "A" else "B"
  dominant <- if (dominant_type == "B") "B"
              else if (n_h1 + n_h3 == 0) NA_character_                 # type A clear, subtype unknowable
              else if (n_h1 >= n_h3) "A(H1N1)" else "A(H3N2)"
  tibble(typed_n = n_A + n_B, n_A = n_A, n_B = n_B,
         dominant_type = dominant_type, dominant = dominant,
         top_share     = max(n_A, n_B) / (n_A + n_B),                  # type-level clarity
         subtype_share = ifelse(dominant_type=="A" & n_h1+n_h3 > 0,    # within-subtyped-A clarity
                                max(n_h1, n_h3) / (n_h1 + n_h3), NA_real_))
}

typed <- typing %>% filter(!is.na(type_group))
dominant_by_season <- typed %>%
  filter(countryname != "EU/EEA") %>%
  mutate(country = EU_short(countryname)) %>% filter(!is.na(country)) %>%
  group_by(country, season) %>% group_modify(~ dominance_call(.x)) %>% ungroup() %>%
  filter(typed_n >= MIN_TYPED_DETECTIONS)
write.csv(dominant_by_season, "output/dominant_subtype.csv", row.names=FALSE)

# continental call from the ERVISS 'EU/EEA' aggregate rows (the cross-check for
# data/external/dominant_subtype_by_season.csv, which drives the season-level analyses)
dominant_continental <- typed %>% filter(countryname == "EU/EEA") %>%
  group_by(season) %>% group_modify(~ dominance_call(.x)) %>% ungroup()
write.csv(dominant_continental, "output/dominant_subtype_continental.csv", row.names=FALSE)

descriptors <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE)
descriptors_subtype <- descriptors %>%
  left_join(dominant_by_season %>% select(country, season, dominant, typed_n, top_share), by=c("country","season"))
write.csv(descriptors_subtype, "output/descriptors_subtype.csv", row.names=FALSE)

cat("Dominant subtype by season (countries, typed_n>=20; NA = type A dominant but no subtyped A):\n")
print(dominant_by_season %>% count(season, dominant) %>%
        tidyr::pivot_wider(names_from=dominant, values_from=n, values_fill=0) %>% as.data.frame(), row.names=FALSE)
cat("\nContinental (ERVISS 'EU/EEA' aggregate):\n")
print(dominant_continental %>%
        select(season, typed_n, n_A, n_B, dominant_type, dominant, top_share, subtype_share) %>%
        as.data.frame(), row.names=FALSE, digits=3)
cat(sprintf("\n'clear' at type level (top share >=0.5 of A+B): %.0f%%  | median type-level top share: %.2f\n",
            100*mean(dominant_by_season$top_share>=0.5), median(dominant_by_season$top_share)))
cat(sprintf("Descriptor country-seasons with a determined dominant subtype: %d / %d (%.0f%%)\n",
            sum(!is.na(descriptors_subtype$dominant)), nrow(descriptors_subtype), 100*mean(!is.na(descriptors_subtype$dominant))))
