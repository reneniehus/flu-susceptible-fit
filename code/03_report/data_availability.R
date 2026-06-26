# data_availability.R
#
# The combined influenza ILI+ analysis dataset at a glance: for each country x season, the source we
# would use (RespiCompass where available, i.e. up to 2023/24; ERVISS reconstruction for 2024/25+),
# the number of observed weeks, and which country-seasons are EXCLUDED -- the four COVID-impacted
# seasons (greyed) and any cell with < 15 observed weeks (red X). The vertical divider after 2023/24
# marks the RespiCompass -> ERVISS hand-over (2023/24 is the overlap used to align the two streams).
#
# ERVISS ILI+ is reconstructed RespiCompass-style: ILI consultation rate x influenza SENTINEL test
# positivity, except NON-sentinel positivity for MT, IS, HR, RO, LV, FI (per the RespiCompass note).
#
# Run from the repo root:  Rscript code/03_report/data_availability.R

source("code/01_main_supporting/setup.R")
models_in <- readRDS("output/models_in.rds")

covid    <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")   # excluded (acute COVID phase)
min_wk   <- 15                                                       # min observed weeks to include a season
nonsent  <- c("MT", "IS", "HR", "RO", "LV", "FI")                   # non-sentinel positivity countries

ip <- models_in$data_timeseries_long %>%
  filter(indicator=="ili_plus", pathogen=="Influenza", agegroup=="age_total")

erv <- ip %>% filter(stream %in% c("ili_plus_sentinel","ili_plus_nonsentinel")) %>%
  mutate(pick = ifelse(country_short %in% nonsent, "ili_plus_nonsentinel", "ili_plus_sentinel")) %>%
  filter(stream==pick) %>%
  group_by(country_short, season) %>% summarise(erv_wk=sum(is.finite(value) & value>0), .groups="drop")
rsp <- ip %>% filter(stream=="ili_plus_respicompass") %>%
  group_by(country_short, season) %>% summarise(rsp_wk=sum(is.finite(value) & value>0), .groups="drop")

av <- full_join(erv, rsp, by=c("country_short","season")) %>%
  mutate(erv_wk=coalesce(erv_wk,0L), rsp_wk=coalesce(rsp_wk,0L),
         source = ifelse(rsp_wk>0, "RespiCompass", "ERVISS"),       # RespiCompass where it exists, else ERVISS
         weeks  = ifelse(rsp_wk>0, rsp_wk, erv_wk)) %>%
  filter(weeks>0) %>%
  mutate(status = case_when(season %in% covid ~ "excluded: COVID",
                            weeks < min_wk      ~ "excluded: <15 wks",
                            TRUE                ~ "included"))

lvl <- sort(unique(av$season)); av$season <- factor(av$season, levels=lvl)
av$country_short <- factor(av$country_short, levels=rev(sort(unique(av$country_short))))
xcov <- match(intersect(covid, lvl), lvl)
xdiv <- match("2023/2024", lvl) + 0.5                                # RespiCompass | ERVISS hand-over

p <- ggplot(av, aes(season, country_short)) +
  annotate("rect", xmin=xcov-0.5, xmax=xcov+0.5, ymin=-Inf, ymax=Inf, fill="grey60", alpha=0.35) +
  geom_tile(aes(fill=weeks, color=source), linewidth=0.6, width=0.92, height=0.92) +
  geom_text(data=subset(av, status=="excluded: <15 wks"), aes(label="x"), color="red", size=3) +
  geom_vline(xintercept=xdiv, linetype="dashed", color="grey20") +
  scale_fill_viridis_c(name="observed\nweeks", option="D", limits=c(0,40), oob=scales::squish) +
  scale_color_manual(name="source", values=c(RespiCompass="#1b9e77", ERVISS="#d95f02")) +
  labs(title="Combined influenza ILI+ analysis dataset: coverage and exclusions",
       subtitle="tile = observed weeks | border = source (RespiCompass <=2023/24, ERVISS 2024/25+) | grey = COVID excluded | red x = <15 weeks | dashed = stream hand-over",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())

dir.create("output", showWarnings=FALSE)
ggsave("output/data_availability.png", p, width=12, height=8.5, dpi=110)
cat("figure -> output/data_availability.png\n\n")

cat("Included country-seasons per season (>=15 weeks, non-COVID):\n")
av %>% filter(status=="included") %>% count(season, source, name="n") %>%
  tidyr::pivot_wider(names_from=source, values_from=n, values_fill=0) %>% as.data.frame() %>% print(row.names=FALSE)
