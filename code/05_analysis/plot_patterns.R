# plot_patterns.R -- analysis figures for the descriptor patterns + vaccination linkage.
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
descriptors <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>%
  mutate(log_auc=log(auc),   # era by SEASON (2023/24 is post-COVID even where RespiCompass-sourced); source kept separately
         era=ifelse(covid_era=="post", "post-COVID (2023/24 on)", "pre-COVID"))
within_country <- descriptors %>% group_by(country) %>% mutate(steep_c=steepness-mean(steepness), log_auc_c=log_auc-mean(log_auc),
                                         onset_c=onset_week-mean(onset_week)) %>% ungroup()
country_means <- descriptors %>% group_by(country) %>% summarise(onset=mean(onset_week), steep=mean(steepness), peak_wk=mean(peak_week), .groups="drop")
plot_theme <- theme_minimal(base_size=10)
cor_burden_speed <- round(cor(within_country$log_auc_c, within_country$steep_c, use="pairwise"), 2)    # within-country burden vs speed (computed live so titles never drift)
cor_onset_steep <- round(cor(country_means$onset, country_means$steep), 2)                        # across-country mean onset vs mean steepness

p1 <- ggplot(within_country, aes(log_auc_c, steep_c)) + geom_point(alpha=0.5, color="#1b9e77") + geom_smooth(method="lm",se=FALSE,color="grey30") +
  labs(title=sprintf("Within-country: steepness vs burden (decoupled, r=%.2f)", cor_burden_speed), x="log AUC (country-demeaned)", y="steepness (demeaned)") + plot_theme
p2 <- ggplot(country_means, aes(onset, steep)) + geom_point(size=2, color="#d95f02") + geom_smooth(method="lm",se=FALSE,color="grey30") +
  geom_text(aes(label=country), size=2.5, vjust=-0.6) +
  labs(title=sprintf("Across countries: later onset -> steeper (r=%.2f)", cor_onset_steep), x="mean onset week", y="mean steepness") + plot_theme
p3 <- descriptors %>% group_by(season, era) %>%
  summarise(onset=mean(onset_week), source_mix=paste(sort(unique(source)), collapse="+"), .groups="drop") %>%
  ggplot(aes(season, onset, color=era, shape=source_mix, group=1)) + geom_point(size=3) +
  labs(title="Season onset over time (post-COVID earlier)", shape="source(s)",
       subtitle="era by season; point shape shows the measurement source(s), which shift one season later (2024/25)") +
  labs(x=NULL, y="mean onset week") + plot_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top", plot.subtitle=element_text(size=7.5))
vax_country_means <- descriptors %>% filter(is.finite(vax_cov_65)) %>% group_by(country) %>%
  summarise(vax=mean(vax_cov_65), peak_wk=mean(peak_week), steep=mean(steepness), .groups="drop")
cor_vax_peak <- round(cor(vax_country_means$vax, vax_country_means$peak_wk), 2)                        # across-country 65+ coverage vs peak timing
p4 <- ggplot(vax_country_means, aes(vax, peak_wk)) + geom_point(size=2, color="#7570b3") + geom_smooth(method="lm",se=FALSE,color="grey30") +
  labs(title=sprintf("Vaccination 65+ vs peak timing (across countries, r=%.2f)", cor_vax_peak), x="mean 65+ coverage", y="mean peak week") + plot_theme
ggsave("output/analysis_patterns.png", (p1|p2)/(p3|p4), width=13, height=9, dpi=110)
cat("figure -> output/analysis_patterns.png\n")
