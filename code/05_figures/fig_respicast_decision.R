# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Two neutral figures for the RespiCast decision-annex (v2) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# For an internal decision paper the figures are deliberately plain -- no baked-in
# titles (the caption in the document carries the neutral description and source).
# Only the two figures the review said to keep: in-house capacity, and the decision
# pathway. Written to output/figures/respicast/.
# Run (after code/00_main.R):  Rscript code/05_figures/fig_respicast_decision.R

source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()
fig_dir <- file.path(params$figure_dir, "respicast")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
od <- params$output_dir

INK <- "#333844"; SLATE <- "#4a5a70"
theme_plain <- function() theme_project() + theme(plot.margin = margin(6, 14, 4, 4),
                                                  panel.grid.major.y = element_blank())

# ---- |-Fig A: in-house modelling capacity (the need) ----
staff <- read_csv(file.path(od, "survey_staff.csv"), show_col_types = FALSE) %>%
  mutate(label = recode(staff, "0 staff" = "No in-house\nmodelling team",
                        "1-5 staff" = "1-5 staff", ">10 staff" = ">10 staff"),
         label = factor(label, levels = c("No in-house\nmodelling team", "1-5 staff", ">10 staff")))
pA <- ggplot(staff, aes(label, n)) +
  geom_col(width = 0.6, fill = SLATE) +
  geom_text(aes(label = sprintf("%d  (%d%%)", n, pct)), vjust = -0.4, size = 3.5, colour = INK) +
  scale_y_continuous(limits = c(0, 15.5), expand = c(0, 0)) +
  labs(x = NULL, y = "countries (of 19)") +
  theme_plain()
ggsave(file.path(fig_dir, "dec_fig1_capacity.png"), pA, width = 6.6, height = 2.7, dpi = 150)

# ---- |-Fig B: a clear mechanism to integrate modelling into decisions (the ceiling) ----
lev  <- c("Strongly disagree", "Disagree", "Neither agree nor disagree", "Agree", "Strongly agree")
lcol <- c("Strongly disagree" = "#2f5c93", "Disagree" = "#7fa8d4",
          "Neither agree nor disagree" = "#c9cdd7", "Agree" = "#e6a06a", "Strongly agree" = "#c96f43")
integ <- read_csv(file.path(od, "survey_integration.csv"), show_col_types = FALSE) %>%
  mutate(level = factor(level, levels = lev)) %>% arrange(level) %>%
  mutate(short = recode(as.character(level), "Neither agree nor disagree" = "Neither"))
pB <- ggplot(integ, aes(x = 1, y = n, fill = level)) +
  geom_col(width = 0.42, colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%s\n%d", short, n)), position = position_stack(vjust = 0.5),
            size = 2.8, colour = "white") +
  scale_fill_manual(values = lcol, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "countries (of 19)") +
  theme_plain() + theme(axis.text.y = element_blank(), panel.grid = element_blank())
ggsave(file.path(fig_dir, "dec_fig2_pathway.png"), pB, width = 6.6, height = 1.9, dpi = 150)

cat("2 neutral figures -> output/figures/respicast/dec_fig{1,2}_*.png\n")
