# hierarchical_models.R -- WITHIN-country partial-pooling models (lme4) for the descriptor questions.
# Each predictor is group-mean-centered within country, so the fixed slope is the WITHIN-country effect;
# country random intercepts + random slopes give partial pooling (shared EU/EEA slope learns from all
# countries) and estimate across-country heterogeneity. Variables are standardised (log for AUC/peak),
# so slopes are in SD units (~ within-country correlation). Run from the repo root.
suppressMessages({library(dplyr); library(lme4); library(ggplot2)})
descriptors <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE)
LOG_VARS <- c("auc","peak_height")     # reporting-scale-dependent outcomes -> modelled on the log scale

fit_pp <- function(df, outcome, predictor){
  df <- df[is.finite(df[[outcome]]) & is.finite(df[[predictor]]), ]
  y_trans <- if (outcome   %in% LOG_VARS) log(df[[outcome]])   else df[[outcome]]
  x_trans <- if (predictor %in% LOG_VARS) log(df[[predictor]]) else df[[predictor]]
  df$Y <- as.numeric(scale(y_trans)); df$pred_z <- as.numeric(scale(x_trans))
  df$pred_within <- df$pred_z - ave(df$pred_z, df$country)      # group-mean-centre the predictor within country
  model <- suppressWarnings(tryCatch(lmer(Y ~ pred_within + (pred_within|country), df, REML=TRUE), error=function(e) NULL))
  slope_sd <- NA_real_; model_type <- "ri+rs"
  if (is.null(model) || isSingular(model)){ model <- lmer(Y ~ pred_within + (1|country), df, REML=TRUE); model_type <- "ri" }
  else slope_sd <- sqrt(VarCorr(model)$country["pred_within","pred_within"])
  coefs <- summary(model)$coefficients; est <- coefs["pred_within","Estimate"]; se <- coefs["pred_within","Std. Error"]
  data.frame(outcome=outcome, predictor=predictor, beta=round(est,3), lo=round(est-1.96*se,3), hi=round(est+1.96*se,3),
             het_sd=round(slope_sd,2), model=model_type, n=nrow(df), n_ctry=length(unique(df$country)),
             sig=ifelse(abs(est)>1.96*se,"*",""))
}

cat("=== Q1: which 2 of {steepness, peak_height, AUC} share most info (within-country std slope) ===\n")
q1 <- rbind(fit_pp(descriptors,"peak_height","auc"), fit_pp(descriptors,"peak_height","steepness"), fit_pp(descriptors,"auc","steepness"))
print(q1, row.names=FALSE)

cat("\n=== Q2: what does onset_week predict (within-country) ===\n")
q2 <- rbind(fit_pp(descriptors,"peak_height","onset_week"), fit_pp(descriptors,"auc","onset_week"),
            fit_pp(descriptors,"peak_week","onset_week"), fit_pp(descriptors,"steepness","onset_week"))
print(q2, row.names=FALSE)

cat("\n=== Q3: what does 65+ vaccination coverage predict (within-country, pre-COVID) ===\n")
vax_rows <- descriptors[is.finite(descriptors$vax_cov_65), ]
q3 <- rbind(fit_pp(vax_rows,"peak_height","vax_cov_65"), fit_pp(vax_rows,"auc","vax_cov_65"),
            fit_pp(vax_rows,"onset_week","vax_cov_65"), fit_pp(vax_rows,"peak_week","vax_cov_65"),
            fit_pp(vax_rows,"steepness","vax_cov_65"))
print(q3, row.names=FALSE)

# forest plot of Q2 + Q3
forest <- rbind(cbind(q2, Q="Q2: onset_week predicts ->"), cbind(q3, Q="Q3: 65+ coverage predicts ->")) %>%
  mutate(lab=paste0(outcome))
p <- ggplot(forest, aes(beta, lab)) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi), color="#1b9e77") +
  facet_wrap(~Q, scales="free_y") + labs(x="within-country standardised slope (95% CI)", y=NULL,
  title="Partial-pooling within-country effects (lme4)") + theme_minimal(base_size=11)
ggsave("output/hierarchical_effects.png", p, width=11, height=4.5, dpi=110)
cat("\nfigure -> output/hierarchical_effects.png\n")
