# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Static figure companions for the forecast-coverage analysis ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The HTML dashboard is the interactive home of these views; this script renders the
# two headline coverage figures as PNGs for the repo / a paper / a slide.
# Run (after code/00_main.R):  Rscript code/05_figures/fig_coverage.R

source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()
dir.create(params$figure_dir, showWarnings = FALSE, recursive = TRUE)

weekly <- read_csv(file.path(params$output_dir, "hub_coverage_weekly.csv"), show_col_types = FALSE) %>%
  mutate(origin_date = as.Date(origin_date),
         indicator   = factor(indicator, levels = c("COVID-19 hospitalisations", "ILI incidence", "ARI incidence")))

# ---- |-Fig 1: the coverage grid (indicator x week, fill = #models) ----
p_heat <- ggplot(weekly, aes(origin_date, indicator, fill = n_models)) +
  geom_tile(width = 6, height = 0.82) +
  # mark weeks that carried an ensemble with a thin tick along the bottom of each row
  geom_point(data = filter(weekly, has_ensemble), aes(y = as.integer(indicator) - 0.44),
             shape = 15, size = 0.7, colour = "#0b6f8c") +
  scale_fill_gradientn(colours = c("#e6f0fb", "#9ec5f4", "#3f8ae2", "#245fb0", "#123f7d"),
                       name = "models", limits = c(0, max(weekly$n_models))) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
  labs(title = "RespiCast forecast coverage: models per indicator per week",
       subtitle = "Fill = contributing models (ensemble/baseline excluded); teal ticks = weeks an ensemble was published",
       x = NULL, y = NULL) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

# ---- |-Fig 2: models per week, by indicator ----
p_line <- ggplot(weekly, aes(origin_date, n_models, colour = indicator)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = palette_indicator, name = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %y") +
  labs(title = "Contributing models per week",
       subtitle = "ILI draws ~2x the models of ARI in the same hub; winter peaks, summer troughs",
       x = NULL, y = "models") +
  theme_project()

fig <- p_heat / p_line + patchwork::plot_layout(heights = c(1, 1.1))
ggsave(file.path(params$figure_dir, "coverage.png"), fig, width = 11, height = 8.5, dpi = 120)
cat("figure -> output/figures/coverage.png\n")
