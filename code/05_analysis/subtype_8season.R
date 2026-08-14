# subtype_8season.R
#
# Extends the dominant-subtype -> descriptor analysis to ALL EIGHT panel seasons. Dominance labels
# are HYBRID: the country's OWN typed call where ERVISS typing supports one (post-COVID seasons;
# code/05_analysis/dominant_subtype.R, hierarchical type-first rule), and the CONTINENTAL
# season label otherwise (pre-COVID seasons, from literature: data/external/
# dominant_subtype_by_season.csv). Using real country-level calls where they exist both respects
# the data and injects some WITHIN-season subtype variation (a minority of countries depart from
# the continental label), which softens the subtype~season collinearity.
#
# WHAT 8 SEASONS BUY -- STATED HONESTLY (this changed when the dominance rule was corrected):
# A(H1N1) and A(H3N2) each recur in BOTH covid eras (H1N1: 2015/16, 2018/19 pre + 2023/24,
# 2024/25 post; H3N2: 2014/15, 2016/17 pre + 2025/26 post), so their contrast is genuinely
# de-confounded from the era. B does NOT recur: it dominates only 2017/18 (pre) plus a handful of
# post-COVID country-level calls, so B contrasts remain nearly a single-season reading -- treat
# them as such. (An earlier header claimed all three subtypes recur in both eras; that was an
# artifact of the miscounted 2024/25 'B' call and of era being coded from the data SOURCE.)
#
# MODEL (the honest-replication upgrade): y = X beta + u_country + w_season + eps. The season
# random intercept matters: subtype/era are season-level, so with a country intercept only the
# ~20+ countries sharing a season would be treated as independent replicates and the CrIs would
# be anti-conservative -- w_season restores effective replication = #seasons (8). era = covid_era
# (BY SEASON, 2023/24 on = post; the measurement SOURCE shifts one season later and its remaining
# confounding is absorbed by the season intercepts). Gibbs (3 chains, shared helpers),
# cross-checked against lme4 with crossed REs. Outcomes standardised (z; log for AUC/peak).
# Run from the repo root:  Rscript code/05_analysis/subtype_8season.R
suppressMessages({library(dplyr); library(lme4)}); set.seed(1)
source("code/05_analysis/analysis_helpers.R")

desc <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE)
sub  <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE) %>%
  transmute(season, dominant_cont = dominant)
cty  <- read.csv("output/dominant_subtype.csv", stringsAsFactors=FALSE) %>%
  filter(!is.na(dominant)) %>% select(country, season, dominant_cty = dominant)

d <- desc %>% inner_join(sub, by="season") %>% left_join(cty, by=c("country","season")) %>%
  mutate(dominant = coalesce(dominant_cty, dominant_cont),      # country call where typed, continental otherwise
         auc      = log(auc), peak_height = log(peak_height),
         dominant = factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")),
         era      = factor(covid_era, levels=c("pre","post")))
cat(sprintf("n = %d country-seasons | %d countries | %d seasons | country-level dominance for %d rows\n",
            nrow(d), n_distinct(d$country), n_distinct(d$season), sum(!is.na(d$dominant_cty))))
cat("subtype x covid_era table (country-seasons; note B is essentially pre-only):\n"); print(table(d$dominant, d$era))

g_country <- as.integer(factor(d$country))
g_season  <- as.integer(factor(d$season))
outcome_labels <- c(auc="AUC (log)", peak_height="peak height (log)", peak_week="peak week",
          onset_week="onset week", steepness="steepness")
res <- list()
for (outcome_key in names(outcome_labels)){
  y <- as.numeric(scale(d[[outcome_key]])); X <- model.matrix(~ dominant + era, d)   # ref = A(H1N1); era = covid era (by season)
  ch <- gibbs_ri2(y, X, g_country, g_season); beta_draws <- do.call(rbind, ch); rh <- rhat_cols(ch)
  dominant_cols <- grep("^dominant", colnames(X))           # columns 2 = H3N2-H1N1, 3 = B-H1N1
  draws <- cbind(`H3N2 - H1N1`=beta_draws[,dominant_cols[1]], `B - H1N1`=beta_draws[,dominant_cols[2]], `B - H3N2`=beta_draws[,dominant_cols[2]]-beta_draws[,dominant_cols[1]])
  rhs   <- c(rh[dominant_cols[1]], rh[dominant_cols[2]],
             rhat_derived(ch, function(M) M[,dominant_cols[2]] - M[,dominant_cols[1]]))
  for (k in seq_along(colnames(draws))){
    cn <- colnames(draws)[k]; q <- q95(draws[,cn])
    res[[length(res)+1]] <- data.frame(outcome=outcome_labels[outcome_key], contrast=cn, est=round(q[1],2),
      lo=round(q[2],2), hi=round(q[3],2), excl0=ifelse(q[2]>0|q[3]<0,"*",""), rhat=round(rhs[k],3))
  }
  m <- suppressWarnings(lmer(y ~ dominant + era + (1|country) + (1|season), d, REML=TRUE))
  cat(sprintf("  [%-12s] lme4 (H3N2-H1N1, B-H1N1): %.2f, %.2f | gibbs: %.2f, %.2f\n",
              outcome_key, fixef(m)[2], fixef(m)[3], mean(beta_draws[,dominant_cols[1]]), mean(beta_draws[,dominant_cols[2]])))
}
out <- do.call(rbind, res); rownames(out)<-NULL
cat("\n=== Subtype contrasts across 8 seasons (country + SEASON random intercepts, net of covid era; SD units, 95% CrI) ===\n")
print(out, row.names=FALSE)
write.csv(out, "output/subtype_8season_contrasts.csv", row.names=FALSE)

# forest plot
suppressMessages(library(ggplot2))
out$contrast <- factor(out$contrast, levels=rev(c("H3N2 - H1N1","B - H1N1","B - H3N2")))
p <- ggplot(out, aes(est, contrast, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) + facet_wrap(~outcome, nrow=1) +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Dominant-subtype contrasts across 8 seasons (country + season random intercepts, net of covid era)",
       subtitle="Hybrid dominance (country-level where typed, continental otherwise). H1N1/H3N2 recur in both eras; B is essentially single-season (2017/18) -- read B contrasts with that in mind.",
       x="contrast (SD units, 95% CrI)", y=NULL) + theme_minimal(base_size=10)
ggsave("output/subtype_8season.png", p, width=13, height=3.6, dpi=110)
cat("figure -> output/subtype_8season.png\n")
