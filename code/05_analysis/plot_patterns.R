# plot_patterns.R -- analysis figures for the descriptor patterns + vaccination linkage.
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
d <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>%
  mutate(lauc=log(auc), era=ifelse(source=="RespiCompass","pre-COVID (RespiCompass)","post-COVID (ERVISS)"))
dw <- d %>% group_by(country) %>% mutate(steep_c=steepness-mean(steepness), lauc_c=lauc-mean(lauc),
                                         onset_c=onset_week-mean(onset_week)) %>% ungroup()
cm <- d %>% group_by(country) %>% summarise(onset=mean(onset_week), steep=mean(steepness), peak_wk=mean(peak_week), .groups="drop")
th <- theme_minimal(base_size=10)
r1 <- round(cor(dw$lauc_c, dw$steep_c, use="pairwise"), 2)    # within-country burden vs speed (computed live so titles never drift)
r2 <- round(cor(cm$onset, cm$steep), 2)                        # across-country mean onset vs mean steepness

p1 <- ggplot(dw, aes(lauc_c, steep_c)) + geom_point(alpha=0.5, color="#1b9e77") + geom_smooth(method="lm",se=FALSE,color="grey30") +
  labs(title=sprintf("Within-country: steepness vs burden (decoupled, r=%.2f)", r1), x="log AUC (country-demeaned)", y="steepness (demeaned)") + th
p2 <- ggplot(cm, aes(onset, steep)) + geom_point(size=2, color="#d95f02") + geom_smooth(method="lm",se=FALSE,color="grey30") +
  geom_text(aes(label=country), size=2.5, vjust=-0.6) +
  labs(title=sprintf("Across countries: later onset -> steeper (r=%.2f)", r2), x="mean onset week", y="mean steepness") + th
p3 <- d %>% group_by(season, era) %>% summarise(onset=mean(onset_week), .groups="drop") %>%
  ggplot(aes(season, onset, color=era, group=1)) + geom_point(size=3) +
  labs(title="Season onset over time (post-COVID earlier; confounded w/ source)", x=NULL, y="mean onset week") + th +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top")
dv <- d %>% filter(is.finite(vax_cov_65)) %>% group_by(country) %>%
  summarise(vax=mean(vax_cov_65), peak_wk=mean(peak_week), steep=mean(steepness), .groups="drop")
r4 <- round(cor(dv$vax, dv$peak_wk), 2)                        # across-country 65+ coverage vs peak timing
p4 <- ggplot(dv, aes(vax, peak_wk)) + geom_point(size=2, color="#7570b3") + geom_smooth(method="lm",se=FALSE,color="grey30") +
  labs(title=sprintf("Vaccination 65+ vs peak timing (across countries, r=%.2f)", r4), x="mean 65+ coverage", y="mean peak week") + th
ggsave("output/analysis_patterns.png", (p1|p2)/(p3|p4), width=13, height=9, dpi=110)
cat("figure -> output/analysis_patterns.png\n")
