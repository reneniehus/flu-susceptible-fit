# plot_vax_scatter.R -- per-country-normalised 65+ coverage vs each descriptor, all country-seasons
# pooled into one panel per outcome, with the within-country correlation line. (AUC/peak log first.)
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
zc <- function(x){ if (sum(is.finite(x))<2 || sd(x,na.rm=TRUE)==0) return(rep(NA_real_,length(x)))
                   (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE) }
d <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>%
  filter(is.finite(vax_cov_65)) %>%
  mutate(auc=log(auc), peak_height=log(peak_height)) %>%
  group_by(country) %>%
  mutate(across(c(vax_cov_65, steepness, peak_week, peak_height, onset_week, auc), zc)) %>% ungroup()

long <- d %>% select(country, season, vax_cov_65, steepness, peak_week, peak_height, onset_week, auc) %>%
  pivot_longer(c(steepness, peak_week, peak_height, onset_week, auc), names_to="outcome", values_to="y") %>%
  filter(is.finite(vax_cov_65), is.finite(y))
labs_o <- c(auc="AUC (log, z)", peak_height="peak height (log, z)", peak_week="peak week (z)",
            onset_week="onset week (z)", steepness="steepness (z)")
long$outcome <- factor(long$outcome, levels=names(labs_o), labels=labs_o)
rr <- long %>% group_by(outcome) %>% summarise(r=cor(vax_cov_65,y), n=n(), .groups="drop") %>%
  mutate(lab=sprintf("r = %.2f (n=%d)", r, n))

p <- ggplot(long, aes(vax_cov_65, y)) +
  geom_hline(yintercept=0, color="grey85") + geom_vline(xintercept=0, color="grey85") +
  geom_point(alpha=0.45, size=1.1, color="#1b9e77") +
  geom_smooth(method="lm", se=TRUE, color="#d95f02", fill="#d95f02", alpha=0.15) +
  geom_text(data=rr, aes(x=-Inf, y=Inf, label=lab), hjust=-0.1, vjust=1.4, size=3.4, inherit.aes=FALSE) +
  facet_wrap(~outcome, nrow=2) +
  labs(title="Within-country: 65+ vaccination coverage vs each descriptor (each normalised per country)",
       subtitle="all country-seasons with coverage (pre-COVID) pooled; line = within-country correlation",
       x="65+ coverage (within-country z-score)", y="descriptor (within-country z-score)") +
  theme_minimal(base_size=11)
ggsave("output/vax_scatter.png", p, width=12, height=7.5, dpi=110)
cat("figure -> output/vax_scatter.png\n"); print(as.data.frame(rr[,c("outcome","r","n")]), row.names=FALSE)
