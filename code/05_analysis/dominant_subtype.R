# dominant_subtype.R
#
# Determine the DOMINANT influenza (sub)type per country-season, the WHO/ECDC way: the (sub)type
# accounting for the largest share of CHARACTERISED detections over the season, among A(H1N1)pdm09,
# A(H3N2), and B. We combine sentinel + non-sentinel detections (all characterised viruses), and
# require >=20 typed detections in the country-season for a reliable call (sparse country-seasons are
# left unassigned). This mirrors ECDC/WHO "predominant (sub)type" reporting (plurality of typed
# detections; we also record the leading-type share so "clear" vs "mixed/co-circulating" is transparent).
#
# Subtype data is ERVISS typing, which starts 2021 -> only the post-COVID analysis seasons get a
# subtype; the pre-COVID (RespiCompass) seasons have none. Writes output/dominant_subtype.csv and
# output/descriptors_subtype.csv (descriptors + dominant subtype). Run from the repo root.

source("code/01_main_supporting/setup.R")

MIN_TYPED_DETECTIONS <- 20   # a country-season needs at least this many typed detections for a reliable call

typing <- read.csv("data/erviss_flu_type_subtype.csv", stringsAsFactors=FALSE) %>%
  filter(pathogen=="Influenza", indicator=="detections") %>%
  mutate(date        = ISOweek::ISOweek2date(paste0(yearweek, "-3")),
         country     = EU_short(countryname),
         season_year = year(date) - (month(date) < 8), season = paste0(season_year, "/", season_year+1),
         subtype_group = case_when(pathogensubtype == "A(H1)pdm09"                       ~ "A(H1N1)",
                                   pathogensubtype == "A(H3)"                            ~ "A(H3N2)",
                                   pathogensubtype %in% c("B/Vic","B/Yam","B (unknown)") ~ "B",
                                   TRUE                                                  ~ "unsubtyped"))

dominant_by_season <- typing %>% filter(subtype_group %in% c("A(H1N1)","A(H3N2)","B")) %>%
  group_by(country, season, subtype_group) %>% summarise(detections = sum(value, na.rm=TRUE), .groups="drop") %>%
  group_by(country, season) %>%
  summarise(typed_n = sum(detections), dominant = subtype_group[which.max(detections)],
            top_share = max(detections)/sum(detections), .groups="drop") %>%
  filter(typed_n >= MIN_TYPED_DETECTIONS)
write.csv(dominant_by_season, "output/dominant_subtype.csv", row.names=FALSE)

descriptors <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE)
descriptors_subtype <- descriptors %>%
  left_join(dominant_by_season %>% select(country, season, dominant, typed_n, top_share), by=c("country","season"))
write.csv(descriptors_subtype, "output/descriptors_subtype.csv", row.names=FALSE)

cat("Dominant subtype by season (countries, typed_n>=20):\n")
print(dominant_by_season %>% count(season, dominant) %>%
        tidyr::pivot_wider(names_from=dominant, values_from=n, values_fill=0) %>% as.data.frame(), row.names=FALSE)
cat(sprintf("\n'clear' (top share >=0.5): %.0f%%  | median top share: %.2f\n",
            100*mean(dominant_by_season$top_share>=0.5), median(dominant_by_season$top_share)))
cat(sprintf("Descriptor country-seasons with a determined dominant subtype: %d / %d (%.0f%%)\n",
            sum(!is.na(descriptors_subtype$dominant)), nrow(descriptors_subtype), 100*mean(!is.na(descriptors_subtype$dominant))))
