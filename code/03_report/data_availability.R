# data_availability.R
#
# The FINAL combined influenza ILI+ analysis dataset at a glance, encoding all the source decisions:
#   - default stitch: RespiCompass up to 2023/24, ERVISS reconstruction for 2024/25+ (2023/24 is the
#     overlap used to align the two streams);
#   - single-source overrides (most non-COVID seasons / data quality): NO, ES = RespiCompass only;
#     SK, LV = ERVISS only (the other source is dismissed for those countries);
#   - four COVID-impacted seasons excluded (greyed); cells with < 15 observed weeks excluded (red x).
# Border colour = source; the dashed line marks the RespiCompass -> ERVISS hand-over.
#
# ERVISS ILI+ is reconstructed RespiCompass-style: ILI rate x SENTINEL positivity, except NON-sentinel
# for MT, IS, HR, RO, LV, FI; and the per-100-consultations rates of CY, LU, MT are x1000 (units).
# (The x1000 affects values, not the week counts shown here.)
#
# Run from the repo root:  Rscript code/03_report/data_availability.R

source("code/01_main_supporting/setup.R")
source("code/01_main_supporting/stitch_iliplus.R")
models_in <- readRDS("output/models_in.rds")

covid           <- stitch_covid_seasons  # shared with the panel build
min_wk          <- 15
handover_season <- "2023/2024"  # last RespiCompass season / overlap -> RespiCompass->ERVISS divider

# Week counts come from the SAME stitched series the analysis panel uses (min_wk = 0 so excluded
# cells still appear; COVID seasons kept so they can be greyed out). An earlier version counted a
# single stream per cell, so 11 overlap-season cells displayed fewer weeks than the final dataset
# actually contains (RespiCompass gaps that aligned ERVISS fills) -- the union is what is analysed.
stitched <- stitch_iliplus_panel(models_in, exclude_covid=FALSE, min_wk=0)
av <- stitched %>%
  group_by(country_short, season) %>%
  summarise(weeks = sum(is.finite(value) & value>0), source = source[1], .groups="drop") %>%
  filter(weeks > 0) %>%
  mutate(status = case_when(season %in% covid ~ "excluded: COVID",
                            weeks < min_wk      ~ "excluded: <15 wks",
                            TRUE                ~ "included"))

season_levels <- sort(unique(av$season)); av$season <- factor(av$season, levels=season_levels)
av$country_short <- factor(av$country_short, levels=rev(sort(unique(av$country_short))))
covid_x <- match(intersect(covid, season_levels), season_levels)
handover_x <- match(handover_season, season_levels) + 0.5

p <- ggplot(av, aes(season, country_short)) +
  annotate("rect", xmin=covid_x-0.5, xmax=covid_x+0.5, ymin=-Inf, ymax=Inf, fill="grey60", alpha=0.35) +
  geom_tile(aes(fill=weeks, color=source), linewidth=0.6, width=0.92, height=0.92) +
  geom_text(data=subset(av, status=="excluded: <15 wks"), aes(label="x"), color="red", size=3) +
  geom_vline(xintercept=handover_x, linetype="dashed", color="grey20") +
  scale_fill_viridis_c(name="observed\nweeks", option="D", limits=c(0,40), oob=scales::squish) +   # 40 ~ a full season of weekly reports; higher values squish to the top colour
  scale_color_manual(name="source", values=c(RespiCompass="#1b9e77", ERVISS="#d95f02")) +
  labs(title="Final influenza ILI+ analysis dataset: coverage and exclusions",
       subtitle="border = source | grey = COVID excluded | red x = <15 weeks | dashed = stream hand-over | single-source: NO,ES=RespiCompass, SK,LV=ERVISS",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())

dir.create("output", showWarnings=FALSE)
ggsave("output/data_availability.png", p, width=12, height=8.5, dpi=110)
cat("figure -> output/data_availability.png\n\n")

included <- av %>% filter(status=="included")
cat(sprintf("INCLUDED: %d country-seasons | %d distinct countries | %d seasons\n\n",
            nrow(included), n_distinct(included$country_short), n_distinct(included$season)))
cat("Countries per season (by source):\n")
included %>% count(season, source, name="n") %>%
  tidyr::pivot_wider(names_from=source, values_from=n, values_fill=0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>% as.data.frame() %>% print(row.names=FALSE)
cat("\nSeasons per country (n):\n")
print(table(included$country_short))
