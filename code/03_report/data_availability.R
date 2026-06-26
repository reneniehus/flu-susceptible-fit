# data_availability.R
#
# A picture of ILI+ data availability: observed weeks per country x season, for each source
# (ERVISS sentinel reconstruction vs RespiCompass), so the country/season coverage -- and the gap
# between the two sources -- is visible at a glance. COVID-impacted seasons are greyed; the current
# target seasons are boxed. Run from the repo root:  Rscript code/03_report/data_availability.R

source("code/01_main_supporting/setup.R")
models_in <- readRDS("output/models_in.rds")

covid  <- c("2019/2020", "2020/2021", "2022/2023")
target <- c("2017/2018", "2018/2019", "2023/2024", "2024/2025", "2025/2026")

avail <- models_in$data_timeseries_long %>%
  filter(indicator=="ili_plus", pathogen=="Influenza", agegroup=="age_total",
         stream %in% c("ili_plus_respicompass","ili_plus_sentinel")) %>%
  group_by(source = recode(stream, ili_plus_respicompass="RespiCompass ILI+",
                                    ili_plus_sentinel="ERVISS sentinel ILI+ (reconstructed)"),
           country_short, season) %>%
  summarise(obs = sum(is.finite(value) & value>0), .groups="drop") %>%
  filter(obs > 0)

lvl  <- sort(unique(avail$season)); avail$season <- factor(avail$season, levels=lvl)
avail$country_short <- factor(avail$country_short, levels=rev(sort(unique(avail$country_short))))
xc <- match(intersect(covid, lvl), lvl); xt <- match(intersect(target, lvl), lvl)

p <- ggplot(avail, aes(season, country_short, fill=obs)) +
  annotate("rect", xmin=xc-0.5, xmax=xc+0.5, ymin=-Inf, ymax=Inf, fill="grey55", alpha=0.30) +
  geom_tile(color="white", linewidth=0.25) +
  annotate("rect", xmin=xt-0.5, xmax=xt+0.5, ymin=-Inf, ymax=Inf, fill=NA, color="red", linewidth=0.7) +
  facet_wrap(~source, ncol=1) +
  scale_fill_viridis_c(name="observed\nweeks", option="D", limits=c(0,40), oob=scales::squish) +
  labs(title="Influenza ILI+ data availability by country, season and source",
       subtitle="red box = target seasons | grey = COVID-impacted (excluded) | tile colour = observed weeks",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(),
        strip.text=element_text(face="bold"))

dir.create("output", showWarnings=FALSE)
ggsave("output/data_availability.png", p, width=11, height=10, dpi=110)
cat("figure -> output/data_availability.png\n\n")

cat("Target seasons: countries with >=15 observed weeks, by source:\n")
avail %>% filter(season %in% target, obs>=15) %>%
  count(source, season, name="n_countries") %>%
  tidyr::pivot_wider(names_from=source, values_from=n_countries, values_fill=0) %>%
  as.data.frame() %>% print(row.names=FALSE)
