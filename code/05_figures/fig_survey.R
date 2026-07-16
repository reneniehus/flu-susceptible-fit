# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Static figure companions for the NFP survey analysis ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PNG companions to the interactive dashboard: the engagement ranking, the capacity
# cross-tab, and the Q7 decision-relevance Likert.
# Run (after code/00_main.R):  Rscript code/05_figures/fig_survey.R

source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()
dir.create(params$figure_dir, showWarnings = FALSE, recursive = TRUE)

od <- params$output_dir
engagement <- read_csv(file.path(od, "survey_engagement.csv"),          show_col_types = FALSE)
capacity   <- read_csv(file.path(od, "survey_engagement_capacity.csv"), show_col_types = FALSE)
decisions  <- read_csv(file.path(od, "survey_decisions_dist.csv"),      show_col_types = FALSE)

class_col <- c("Short-term forecasting" = "#2a78d6", "Scenario modelling" = "#4a3aa7", "Guidance" = "#8b92a2")

# ---- |-Fig A: engagement by output, coloured by class ----
p_eng <- engagement %>%
  mutate(output = reorder(output, mean_score)) %>%
  ggplot(aes(mean_score, output, fill = class)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", mean_score)), hjust = -0.2, size = 3, colour = "#52514e") +
  scale_fill_manual(values = class_col, name = NULL) +
  scale_x_continuous(limits = c(0, 5.2), expand = c(0, 0)) +
  labs(title = "Engagement with ECDC outputs (Q5)",
       subtitle = "Mean 0-5; ordinal, read the ranking. Newer products (RespiCast, RespiCompass) sit lower.",
       x = "mean engagement", y = NULL) +
  theme_project() + theme(legend.position = "top")

# ---- |-Fig B: engagement by class x capacity band ----
p_cap <- capacity %>%
  mutate(staff = factor(staff, levels = params$levels_staff),
         class = factor(class, levels = names(class_col))) %>%
  ggplot(aes(staff, mean_score, fill = class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_fill_manual(values = class_col, name = NULL) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  labs(title = "Engagement by in-house capacity (Q5 x Q6)",
       subtitle = "Short-term forecasting is roughly flat across capacity - the need is universal, the uptake is not",
       x = "in-house modelling staff", y = "mean engagement") +
  theme_project() + theme(legend.position = "top")

# ---- |-Fig C: Q7 decision-relevance, stacked Likert ----
lik_lv  <- params$levels_q7
lik_col <- c("Very unlikely" = "#184f95", "Unlikely" = "#7fb0ec", "Unsure" = "#c9cdd7",
             "Likely" = "#eb8c6a", "Very likely" = "#c0392b")
p_dec <- decisions %>%
  filter(action != "Other activities") %>%
  mutate(level = factor(level, levels = lik_lv)) %>%
  group_by(action) %>% mutate(order_key = sum(pct[level %in% c("Likely", "Very likely")])) %>% ungroup() %>%
  mutate(action = reorder(action, order_key)) %>%
  ggplot(aes(pct, action, fill = level)) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = lik_col, name = NULL) +
  labs(title = "Which decisions would ECDC modelling inform? (Q7)",
       subtitle = "Surveillance leads; healthcare-capacity planning - the canonical forecast use - trails",
       x = "% of respondents", y = NULL) +
  theme_project() + theme(legend.position = "top")

fig <- (p_eng / p_cap / p_dec) + patchwork::plot_layout(heights = c(1, 1, 1))
ggsave(file.path(params$figure_dir, "survey.png"), fig, width = 10, height = 12, dpi = 120)
cat("figure -> output/figures/survey.png\n")
