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

ts <- read.csv("data/erviss_flu_type_subtype.csv", stringsAsFactors=FALSE) %>%
  filter(pathogen=="Influenza", indicator=="detections") %>%
  mutate(date    = ISOweek::ISOweek2date(paste0(yearweek, "-3")),
         country = EU_short(countryname),
         syr     = year(date) - (month(date) < 8), season = paste0(syr, "/", syr+1),
         cat = case_when(pathogensubtype == "A(H1)pdm09"                       ~ "A(H1N1)",
                         pathogensubtype == "A(H3)"                            ~ "A(H3N2)",
                         pathogensubtype %in% c("B/Vic","B/Yam","B (unknown)") ~ "B",
                         TRUE                                                  ~ "unsubtyped"))

dom <- ts %>% filter(cat %in% c("A(H1N1)","A(H3N2)","B")) %>%
  group_by(country, season, cat) %>% summarise(det = sum(value, na.rm=TRUE), .groups="drop") %>%
  group_by(country, season) %>%
  summarise(typed_n = sum(det), dominant = cat[which.max(det)],
            top_share = max(det)/sum(det), .groups="drop") %>%
  filter(typed_n >= 20)
write.csv(dom, "output/dominant_subtype.csv", row.names=FALSE)

desc <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE)
mrg  <- desc %>% left_join(dom %>% select(country, season, dominant, typed_n, top_share),
                           by=c("country","season"))
write.csv(mrg, "output/descriptors_subtype.csv", row.names=FALSE)

cat("Dominant subtype by season (countries, typed_n>=20):\n")
print(dom %>% count(season, dominant) %>%
        tidyr::pivot_wider(names_from=dominant, values_from=n, values_fill=0) %>% as.data.frame(), row.names=FALSE)
cat(sprintf("\n'clear' (top share >=0.5): %.0f%%  | median top share: %.2f\n",
            100*mean(dom$top_share>=0.5), median(dom$top_share)))
cat(sprintf("Descriptor country-seasons with a determined dominant subtype: %d / %d (%.0f%%)\n",
            sum(!is.na(mrg$dominant)), nrow(mrg), 100*mean(!is.na(mrg$dominant))))
