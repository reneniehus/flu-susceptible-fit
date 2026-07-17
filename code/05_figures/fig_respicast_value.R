# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Figures for the RespiCast-value supplement ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Five stand-alone PNGs, each with a plain-language title so they read without knowing
# the survey questions. They back the "what value does RespiCast have and could have"
# summary. Written to output/figures/respicast/.
# Run (after code/00_main.R):  Rscript code/05_figures/fig_respicast_value.R

source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()
fig_dir <- file.path(params$figure_dir, "respicast")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
od <- params$output_dir

BLUE <- "#2a78d6"; GREY <- "#b9bec9"; INK <- "#52514e"

# a local theme: title left-aligned across the whole plot (more room, no clipping), plain hyphens only
theme_supp <- function() {
  theme_project() +
    theme(plot.title.position = "plot",
          plot.title    = element_text(size = 12.5, face = "bold"),
          plot.subtitle = element_text(size = 9.3, colour = INK),
          plot.margin   = margin(10, 16, 8, 8))
}
save_fig <- function(p, name, h = 3.4) ggsave(file.path(fig_dir, name), p, width = 7.4, height = h, dpi = 150)

# ---- |-Fig 1: the need -- most institutes have no modelling team of their own ----
staff <- read_csv(file.path(od, "survey_staff.csv"), show_col_types = FALSE) %>%
  mutate(label = recode(staff, "0 staff" = "No in-house\nmodelling team",
                        "1-5 staff" = "1-5 staff", ">10 staff" = ">10 staff"),
         label = factor(label, levels = c("No in-house\nmodelling team", "1-5 staff", ">10 staff")))
p1 <- ggplot(staff, aes(label, n, fill = staff == "0 staff")) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d  (%d%%)", n, pct)), vjust = -0.4, size = 3.6, colour = INK) +
  scale_fill_manual(values = c("TRUE" = BLUE, "FALSE" = GREY), guide = "none") +
  scale_y_continuous(limits = c(0, 16), expand = c(0, 0)) +
  labs(title = "Most national institutes cannot forecast in-house",
       subtitle = "The 19 responding countries, by size of their in-house mathematical-modelling team",
       x = NULL, y = "countries") +
  theme_supp()
save_fig(p1, "fig1_capacity.png")

# ---- |-Fig 2: the preference -- forecasts chosen over scenarios ----
vc <- read_csv(file.path(od, "survey_value_choice.csv"), show_col_types = FALSE) %>%
  mutate(label = recode(value_choice,
                        "Both will be useful" = "Both will be useful",
                        "RespiCast" = "Short-term forecasts\n(RespiCast) only",
                        "Neither will be useful" = "Neither",
                        "RespiCompass" = "Scenarios\n(RespiCompass) only"),
         hl = value_choice %in% c("RespiCast", "Both will be useful")) %>%
  arrange(n) %>% mutate(label = factor(label, levels = label))
p2 <- ggplot(vc, aes(n, label, fill = hl)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d  (%d%%)", n, pct)), hjust = -0.15, size = 3.5, colour = INK) +
  scale_fill_manual(values = c("TRUE" = BLUE, "FALSE" = GREY), guide = "none") +
  scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) +
  labs(title = "Asked which is most useful, forecasts win by a wide margin",
       subtitle = "Choice between ECDC short-term forecasts (RespiCast) and seasonal scenarios (RespiCompass)",
       x = "countries", y = NULL) +
  theme_supp()
save_fig(p2, "fig2_preference.png")

# ---- |-Fig 3: current standing -- engagement with each ECDC output ----
eng <- read_csv(file.path(od, "survey_engagement.csv"), show_col_types = FALSE) %>%
  mutate(nice = recode(output,
                       "COVID-19 risk assessments" = "COVID-19 guidance documents (2020)",
                       "COVID-19 Forecast Hub" = "COVID-19 Forecast Hub (2021)",
                       "COVID-19 Scenario Hub" = "COVID-19 Scenario Hub (2022)",
                       "RespiCast" = "RespiCast forecasts (2023)",
                       "RespiCompass" = "RespiCompass scenarios (2024)"),
         hl = output == "RespiCast") %>%
  arrange(mean_score) %>% mutate(nice = factor(nice, levels = nice))
p3 <- ggplot(eng, aes(mean_score, nice, fill = hl)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", mean_score)), hjust = -0.35, size = 3.6, colour = INK) +
  scale_fill_manual(values = c("TRUE" = BLUE, "FALSE" = GREY), guide = "none") +
  scale_x_continuous(limits = c(0, 4.6), expand = c(0, 0)) +
  labs(title = "The newer forecasting products are less used so far",
       subtitle = "Average self-rated engagement (0 = not aware, 5 = actively used); launch year in brackets; RespiCast in blue",
       x = "average engagement (0-5)", y = NULL) +
  theme_supp()
save_fig(p3, "fig3_engagement.png", h = 3.7)

# ---- |-Fig 4: where value lands -- decisions modelling would inform ----
dec <- read_csv(file.path(od, "survey_decisions_headline.csv"), show_col_types = FALSE) %>%
  filter(action != "Other activities") %>%
  mutate(nice = recode(action,
                       "Inform surveillance activities" = "Surveillance / situational awareness",
                       "Planning vaccination campaigns" = "Vaccination campaigns",
                       "Procuring medical countermeasures" = "Medical countermeasures",
                       "Healthcare capacity planning" = "Healthcare capacity planning")) %>%
  arrange(pct_likely) %>% mutate(nice = factor(nice, levels = nice))
p4 <- ggplot(dec, aes(pct_likely, nice)) +
  geom_col(width = 0.6, fill = BLUE) +
  geom_text(aes(label = sprintf("%d%%", pct_likely)), hjust = -0.25, size = 3.6, colour = INK) +
  scale_x_continuous(limits = c(0, 66), expand = c(0, 0)) +
  labs(title = "Forecasts help awareness most, resource planning least",
       subtitle = "Share of countries saying ECDC modelling is 'likely' or 'very likely' to inform each national decision",
       x = "% of countries", y = NULL) +
  theme_supp()
save_fig(p4, "fig4_decisions.png")

# ---- |-Fig 5: the ceiling -- is there a way to use the modelling? ----
lev <- c("Strongly disagree", "Disagree", "Neither agree nor disagree", "Agree", "Strongly agree")
lcol <- c("Strongly disagree" = "#184f95", "Disagree" = "#7fb0ec",
          "Neither agree nor disagree" = "#c9cdd7", "Agree" = "#eb8c6a", "Strongly agree" = "#c0392b")
integ <- read_csv(file.path(od, "survey_integration.csv"), show_col_types = FALSE) %>%
  mutate(level = factor(level, levels = lev)) %>% arrange(level) %>%
  mutate(short = recode(as.character(level), "Neither agree nor disagree" = "Neither"))
p5 <- ggplot(integ, aes(x = 1, y = n, fill = level)) +
  geom_col(width = 0.5, colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%s\n%d", short, n)), position = position_stack(vjust = 0.5), size = 2.9, colour = "white") +
  scale_fill_manual(values = lcol, guide = "none") +
  coord_flip() +
  labs(title = "Only about half can readily use modelling in decisions",
       subtitle = "'There is a clear mechanism to integrate modelling into decisions' (19 countries)\nAgree: 9 (48%)   |   Disagree: 6 (32%)   |   Neither: 4 (21%)",
       x = NULL, y = "countries") +
  theme_supp() + theme(axis.text.y = element_blank(), panel.grid = element_blank())
save_fig(p5, "fig5_integration.png", h = 2.9)

cat("5 figures -> output/figures/respicast/\n")
