# bayes_subtype.R
#
# Bayesian hierarchical model: does the dominant influenza subtype predict each descriptor, WITHIN
# country? Model:  y_i = X_i beta + u_country[i] + eps_i ,  u ~ N(0, tau^2), eps ~ N(0, sigma^2);
# X = intercept + dominant-subtype dummies (population-level / shared fixed effects), country random
# intercept = partial pooling. Weak priors (beta ~ N(0,100), sigma^2/tau^2 ~ InvGamma(0.01,0.01)).
# Fit by a Gibbs sampler (3 chains); reported as posterior subtype CONTRASTS with 95% credible
# intervals, cross-validated against lme4 fixed effects.
#
# HEAVY CAVEAT: subtype is still largely season-determined here (under the hierarchical type-first
# dominance rule: 2023/24 and 2024/25 mostly A(H1N1), 2025/26 mostly A(H3N2), with a minority of
# countries departing -- see code/05_analysis/dominant_subtype.R), so a "subtype effect" remains
# hard to separate from a SEASON effect (only 3 post-COVID seasons). Read these as among-season
# differences labelled by subtype, NOT causal subtype effects; subtype_8season.R is the wider,
# season-RE version. Outcomes standardised (z; log for AUC/peak) -> contrasts in SD units.
# Run from the repo root.
suppressMessages({library(dplyr); library(lme4)})
set.seed(1)
source("code/05_analysis/analysis_helpers.R")   # gibbs_ri, rhat_cols, q95 (shared, identical draws)

d <- read.csv("output/descriptors_subtype.csv", stringsAsFactors=FALSE) %>%
  filter(!is.na(dominant)) %>%
  mutate(auc=log(auc), peak_height=log(peak_height),
         dominant=factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")))
cat(sprintf("n = %d country-seasons, %d countries, subtypes: %s\n",
            nrow(d), n_distinct(d$country), paste(table(d$dominant), names(table(d$dominant)), collapse="  ")))

g <- as.integer(factor(d$country))
outcome_labels <- c(auc="AUC (log)", peak_height="peak height (log)", peak_week="peak week",
          onset_week="onset week", steepness="steepness")
res <- list()
for (outcome_key in names(outcome_labels)){
  y <- as.numeric(scale(d[[outcome_key]])); X <- model.matrix(~ dominant, d)   # ref = A(H1N1)
  ch <- gibbs_ri(y, X, g); beta_draws <- do.call(rbind, ch); rh <- rhat_cols(ch)
  # contrasts: cols 2=H3N2-H1N1, 3=B-H1N1 ; derive B-H3N2
  draws <- cbind(`H3N2 - H1N1`=beta_draws[,2], `B - H1N1`=beta_draws[,3], `B - H3N2`=beta_draws[,3]-beta_draws[,2])
  for (cn in colnames(draws)){
    q <- quantile(draws[,cn], c(.5,.025,.975))
    res[[length(res)+1]] <- data.frame(outcome=outcome_labels[outcome_key], contrast=cn, est=round(q[1],2),
      lo=round(q[2],2), hi=round(q[3],2), excl0=ifelse(q[2]>0|q[3]<0,"*",""),
      rhat=round(max(rh),3))
  }
  # lme4 cross-check
  m <- suppressWarnings(lmer(y ~ dominant + (1|country), d, REML=TRUE))
  cat(sprintf("  [%s] lme4 fixef (H3N2-H1N1, B-H1N1): %.2f, %.2f | gibbs: %.2f, %.2f\n", outcome_key,
      fixef(m)[2], fixef(m)[3], mean(beta_draws[,2]), mean(beta_draws[,3])))
}
out <- do.call(rbind, res); rownames(out)<-NULL
cat("\n=== Bayesian within-country subtype contrasts (SD units, 95% CrI) ===\n")
print(out, row.names=FALSE)
write.csv(out, "output/bayes_subtype_contrasts.csv", row.names=FALSE)

# forest plot of the subtype contrasts
suppressMessages(library(ggplot2))
forest <- read.csv("output/bayes_subtype_contrasts.csv", stringsAsFactors=FALSE)
forest$contrast <- factor(forest$contrast, levels=rev(c("H3N2 - H1N1","B - H1N1","B - H3N2")))
p <- ggplot(forest, aes(est, contrast, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) + facet_wrap(~outcome, nrow=1) +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Bayesian within-country subtype contrasts (SD units, 95% credible intervals)",
       subtitle="CAVEAT: subtype is ~season-determined (3 post-COVID seasons) -> these are among-SEASON differences labelled by subtype",
       x="contrast (SD units)", y=NULL) + theme_minimal(base_size=10)
ggsave("output/bayes_subtype.png", p, width=13, height=3.6, dpi=110)
cat("figure -> output/bayes_subtype.png\n")
