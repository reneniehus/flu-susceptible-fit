# data_availability.R
#
# The FINAL combined influenza ILI+ analysis dataset at a glance, encoding all the source decisions:
#   - default stitch: RespiCompass up to 2023/24, ERVISS reconstruction for 2024/25+ (2023/24 is the
#     overlap used to align the two streams);
#   - single-source overrides (most non-COVID seasons / data quality): NO, ES = RespiCompass only;
#     SK, LV = ERVISS only (the other source is dismissed for those countries);
#   - four COVID-impacted seasons excluded (greyed); cells with < 15 observed weeks excluded (red X).
# Border colour = source; the dashed line marks the RespiCompass -> ERVISS hand-over.
#
# ERVISS ILI+ is reconstructed RespiCompass-style: ILI rate x SENTINEL positivity, except NON-sentinel
# for MT, IS, HR, RO, LV, FI; and the per-100-consultations rates of CY, LU, MT are x1000 (units).
# (The x1000 affects values, not the week counts shown here.)
#
# Run from the repo root:  Rscript code/03_report/data_availability.R

source("code/01_main_supporting/setup.R")
models_in <- readRDS("output/models_in.rds")

covid       <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")
min_wk      <- 15
nonsent     <- c("MT", "IS", "HR", "RO", "LV", "FI")
resp_only   <- c("NO", "ES")    # single source = RespiCompass
erviss_only <- c("SK", "LV")    # single source = ERVISS

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
         source = case_when(country_short %in% resp_only   ~ ifelse(rsp_wk>0,"RespiCompass",NA_character_),
                            country_short %in% erviss_only ~ ifelse(erv_wk>0,"ERVISS",NA_character_),
                            rsp_wk>0 ~ "RespiCompass", erv_wk>0 ~ "ERVISS", TRUE ~ NA_character_),
         weeks  = case_when(country_short %in% resp_only   ~ rsp_wk,
                            country_short %in% erviss_only ~ erv_wk,
                            rsp_wk>0 ~ rsp_wk, TRUE ~ erv_wk)) %>%
  filter(!is.na(source), weeks>0) %>%
  mutate(status = case_when(season %in% covid ~ "excluded: COVID",
                            weeks < min_wk      ~ "excluded: <15 wks",
                            TRUE                ~ "included"))

lvl <- sort(unique(av$season)); av$season <- factor(av$season, levels=lvl)
av$country_short <- factor(av$country_short, levels=rev(sort(unique(av$country_short))))
xcov <- match(intersect(covid, lvl), lvl)
xdiv <- match("2023/2024", lvl) + 0.5

p <- ggplot(av, aes(season, country_short)) +
  annotate("rect", xmin=xcov-0.5, xmax=xcov+0.5, ymin=-Inf, ymax=Inf, fill="grey60", alpha=0.35) +
  geom_tile(aes(fill=weeks, color=source), linewidth=0.6, width=0.92, height=0.92) +
  geom_text(data=subset(av, status=="excluded: <15 wks"), aes(label="x"), color="red", size=3) +
  geom_vline(xintercept=xdiv, linetype="dashed", color="grey20") +
  scale_fill_viridis_c(name="observed\nweeks", option="D", limits=c(0,40), oob=scales::squish) +
  scale_color_manual(name="source", values=c(RespiCompass="#1b9e77", ERVISS="#d95f02")) +
  labs(title="Final influenza ILI+ analysis dataset: coverage and exclusions",
       subtitle="border = source | grey = COVID excluded | red x = <15 weeks | dashed = stream hand-over | single-source: NO,ES=RespiCompass, SK,LV=ERVISS",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())

dir.create("output", showWarnings=FALSE)
ggsave("output/data_availability.png", p, width=12, height=8.5, dpi=110)
cat("figure -> output/data_availability.png\n\n")

inc <- av %>% filter(status=="included")
cat(sprintf("INCLUDED: %d country-seasons | %d distinct countries | %d seasons\n\n",
            nrow(inc), n_distinct(inc$country_short), n_distinct(inc$season)))
cat("Countries per season (by source):\n")
inc %>% count(season, source, name="n") %>%
  tidyr::pivot_wider(names_from=source, values_from=n, values_fill=0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>% as.data.frame() %>% print(row.names=FALSE)
cat("\nSeasons per country (n):\n")
print(table(inc$country_short))
