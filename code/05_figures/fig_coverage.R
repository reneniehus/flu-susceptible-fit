# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Static figure companions for the forecast-coverage analysis ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The HTML dashboard is the interactive home of these views; this renders the two
# headline coverage figures as PNGs for the repo / a paper / a slide:
#   (1) the COVID-19 hospitalisation continuity track (the "since 2021" claim), and
#   (2) the full five-indicator coverage grid, 2021-2026.
# Run (after code/00_main.R):  Rscript code/05_figures/fig_coverage.R

source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()
dir.create(params$figure_dir, showWarnings = FALSE, recursive = TRUE)

ind_levels <- c("COVID-19 cases", "COVID-19 hospitalisations", "COVID-19 deaths",
                "ILI incidence", "ARI incidence")
weekly <- read_csv(file.path(params$output_dir, "hub_coverage_weekly.csv"), show_col_types = FALSE) %>%
  mutate(week = as.Date(week), indicator = factor(indicator, levels = ind_levels))

# ---- |-Fig 1: COVID-19 hospitalisation continuity, archive -> RespiCast ----
hosp <- weekly %>% filter(indicator == "COVID-19 hospitalisations")
p_cont <- ggplot(hosp, aes(week, n_models, fill = era)) +
  geom_col(width = 6) +
  scale_fill_manual(values = c(archive = "#7a8699", current = "#2a78d6"),
                    labels = c(archive = "EU COVID-19 Forecast Hub", current = "RespiCast-Covid19"),
                    name = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %y") +
  labs(title = "COVID-19 hospitalisation forecasts: almost every week since mid-2021",
       subtitle = "Contributing models per ISO week, stitched across the archive -> RespiCast handover (dashed)",
       x = NULL, y = "models") +
  geom_vline(xintercept = as.Date("2024-10-21"), linetype = "dashed", colour = "#888781") +
  theme_project()

# ---- |-Fig 2: the full coverage grid (indicator x week, fill = #models) ----
p_heat <- ggplot(weekly, aes(week, indicator, fill = n_models)) +
  geom_tile(width = 6, height = 0.8) +
  scale_fill_gradientn(colours = c("#e6f0fb", "#9ec5f4", "#3f8ae2", "#245fb0", "#123f7d"),
                       name = "models", limits = c(0, max(weekly$n_models))) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_discrete(limits = rev(ind_levels)) +
  geom_vline(xintercept = as.Date("2024-10-21"), linetype = "dashed", colour = "#888781") +
  labs(title = "Five indicators, five years: European respiratory forecast coverage",
       subtitle = "Fill = contributing models. COVID cases/deaths end at the 2024 handover; hospitalisations continue; ILI/ARI run winter-only",
       x = NULL, y = NULL) +
  theme_project() + theme(legend.position = "right")

fig <- p_cont / p_heat + patchwork::plot_layout(heights = c(1, 1.15))
ggsave(file.path(params$figure_dir, "coverage.png"), fig, width = 12, height = 8.5, dpi = 120)
cat("figure -> output/figures/coverage.png\n")
