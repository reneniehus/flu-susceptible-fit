# bayes_precovid_ve_subtype.R
#
# PRE-COVID ONLY Bayesian analysis (2014/15-2018/19; one clean era, no RespiCompass->ERVISS or
# behaviour/reporting confound -- see documentation/reflections.md). Within-country partial pooling
# (country random intercept), Gibbs sampler, lme4-cross-checked. Outcomes: log(AUC), onset week,
# log(peak height) -- standardised (z); log for AUC/peak so a MULTIPLICATIVE burden effect is an
# additive shift (consistent with the vaccine mechanism, see reflections.md).
#
# MODEL 1 (as requested): dominant subtype (categorical) + VE (season-level) -> descriptors.
#   Identifiability caveat: subtype and VE are BOTH season-level over only 5 seasons and are correlated
#   (the two A(H3N2) seasons are the two lowest-VE seasons), so they partly compete -- read with care.
#
# MODEL 2 (mechanistic refinement, motivated by reflections.md): dominant subtype + PROTECTION, where
#   protection = VE x coverage(65+)/100 = effectively-protected fraction of the elderly. Because coverage
#   varies by COUNTRY, protection varies by country x season -> far better identified than VE alone.
#   Falsifiable prediction: EU flu vaccines reduce BURDEN not transmission, so protection should push AUC
#   and peak height DOWN but leave onset/peak TIMING unmoved.
#
# Season-level VE against the DOMINANT subtype is DERIVED from the provenance CSV
# (data/external/vaccine_effectiveness.csv) by the explicit preference rule in
# analysis_helpers.R::ve_vs_dominant() -- primary-care all-ages/target-group rows, end-of-season
# point > interim point > study-range midpoint -- so every value used here is traceable to a
# sourced row (2014/15 14.4, 2015/16 32.9, 2016/17 25.7, 2017/18 45 = midpoint of 36-54,
# 2018/19 37.5 = midpoint of the influenza-A 32-43 range).
# Run from the repo root:  Rscript code/05_analysis/bayes_precovid_ve_subtype.R
suppressMessages({library(dplyr); library(lme4)}); set.seed(1)
source("code/05_analysis/analysis_helpers.R")   # gibbs_ri, rhat helpers, q95, ve_vs_dominant

pre <- c("2014/2015","2015/2016","2016/2017","2017/2018","2018/2019")

sub <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE) %>% transmute(season, dominant)
ve_by_season <- ve_vs_dominant(sub)
d <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>%
  filter(season %in% pre) %>% inner_join(sub, by="season") %>%
  mutate(ve = ve_by_season[season],
         protection = ve * vax_cov_65 / 100,                 # effectively-protected fraction of 65+
         dominant = factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")),
         lauc = log(auc), lpk = log(peak_height))
cat(sprintf("pre-COVID: n=%d country-seasons, %d countries, %d seasons | with coverage: %d\n",
            nrow(d), n_distinct(d$country), n_distinct(d$season), sum(is.finite(d$protection))))
cat("\nseason-level subtype and VE (the confounding to keep in mind):\n")
print(distinct(d, season, dominant, ve) %>% arrange(season), row.names=FALSE)

# ---- |-fit subtype + one continuous predictor, report contrasts + slope per outcome ----
# The continuous predictor is GROUP-MEAN-CENTRED within country before scaling, per the standing
# rule (documentation/analysis_strategy.md: 'predictors group-mean-centred; country random
# intercepts') -- a globally-scaled predictor here would carry mostly BETWEEN-country variation
# (coverage differs far more across countries than across years), quietly turning the 'within-
# country' slope into a blended within+between estimate.
fit_model <- function(dat, xvar, xlabel){
  dd <- dat[is.finite(dat[[xvar]]), ]; g <- as.integer(factor(dd$country))
  x_within <- dd[[xvar]] - ave(dd[[xvar]], dd$country)              # within-country deviations
  xz <- x_within / sd(x_within)
  outcome_labels <- c(lauc="AUC (log)", onset_week="onset week", lpk="peak height (log)")
  res <- list()
  for (outcome_key in names(outcome_labels)){
    y <- as.numeric(scale(dd[[outcome_key]]))
    X <- cbind(model.matrix(~ dominant, dd), xz)                    # [1, H3N2, B, predictor]
    ch <- gibbs_ri(y, X, g); beta_draws <- do.call(rbind, ch); rh <- rhat_cols(ch)
    term_draws <- list(`H3N2 - H1N1`=beta_draws[,2], `B - H1N1`=beta_draws[,3], `B - H3N2`=beta_draws[,3]-beta_draws[,2], predictor=beta_draws[,4])
    term_labels <- c(`H3N2 - H1N1`="subtype: H3N2-H1N1", `B - H1N1`="subtype: B-H1N1",
              `B - H3N2`="subtype: B-H3N2", predictor=paste0("slope: ", xlabel))
    term_rhats <- c(`H3N2 - H1N1`=rh[2], `B - H1N1`=rh[3],
                    `B - H3N2`=rhat_derived(ch, function(M) M[,3]-M[,2]),   # rhat of the DERIVED contrast itself
                    predictor=rh[4])
    for (nm in names(term_draws)){ qq <- q95(term_draws[[nm]])
      res[[length(res)+1]] <- data.frame(outcome=outcome_labels[outcome_key], term=term_labels[[nm]], est=round(qq[1],2),
        lo=round(qq[2],2), hi=round(qq[3],2), excl0=ifelse(qq[2]>0|qq[3]<0,"*",""),
        rhat=round(term_rhats[[nm]],3)) }
    m <- suppressWarnings(lmer(y ~ dominant + xz + (1|country), dd, REML=TRUE))
    cat(sprintf("  [%-16s | %-14s] lme4 slope(%s)=%.2f | gibbs=%.2f\n", xlabel, outcome_labels[outcome_key], xvar, fixef(m)["xz"], mean(beta_draws[,4])))
  }
  out <- do.call(rbind, res); rownames(out)<-NULL; out
}

cat("\n=== MODEL 1 (requested): subtype + VE (both season-level; SD units, 95% CrI) ===\n")
m1 <- fit_model(d, "ve", "VE (season)")
print(m1, row.names=FALSE); write.csv(m1, "output/precovid_ve_subtype_model1.csv", row.names=FALSE)

cat("\n=== MODEL 2 (mechanistic): subtype + PROTECTION = VE x coverage (country x season; SD units, 95% CrI) ===\n")
m2 <- fit_model(d, "protection", "protection")
print(m2, row.names=FALSE); write.csv(m2, "output/precovid_ve_subtype_model2.csv", row.names=FALSE)

# ---- |-mechanistic read-out: protection should hit burden (AUC/peak) but not timing (onset) ----
cat("\nMechanistic check (Model 2 protection slope): expect NEGATIVE on AUC/peak, ~0 on onset\n")
protection_slopes <- m2[grepl("slope:", m2$term), c("outcome","est","lo","hi","excl0")]
print(protection_slopes, row.names=FALSE)

# ---- |-forest plot of both models ----
suppressMessages(library(ggplot2))
m1$model <- "Model 1: subtype + VE"; m2$model <- "Model 2: subtype + protection (VE x coverage)"
pd <- rbind(m1, m2); pd$term <- factor(pd$term, levels=rev(unique(pd$term)))
p <- ggplot(pd, aes(est, term, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) + facet_grid(model ~ outcome, scales="free_y", space="free_y") +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Pre-COVID within-country: subtype + VE (M1) and subtype + protection (M2)",
       subtitle="Protection = VE x 65+ coverage varies by country x season -> better identified than season-level VE. Predictors group-mean-centred within country. SD units, 95% CrI.",
       x="effect (SD units)", y=NULL) + theme_minimal(base_size=9)
ggsave("output/precovid_ve_subtype.png", p, width=12, height=6, dpi=110)
cat("\nfigures/tables -> output/precovid_ve_subtype.png ; output/precovid_ve_subtype_model{1,2}.csv\n")
