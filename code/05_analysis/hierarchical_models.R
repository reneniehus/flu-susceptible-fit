# hierarchical_models.R -- WITHIN-country partial-pooling models (lme4) for the descriptor questions.
# Each predictor is group-mean-centered within country, so the fixed slope is the WITHIN-country effect;
# country random intercepts + random slopes give partial pooling (shared EU/EEA slope learns from all
# countries) and estimate across-country heterogeneity. Variables are standardised (log for AUC/peak),
# so slopes are in SD units (~ within-country correlation). Run from the repo root.
suppressMessages({library(dplyr); library(lme4); library(ggplot2)})
d <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE)
LOG <- c("auc","peak_height")

fit_pp <- function(df, yv, xv){
  df <- df[is.finite(df[[yv]]) & is.finite(df[[xv]]), ]
  ty <- if (yv %in% LOG) log(df[[yv]]) else df[[yv]]
  tx <- if (xv %in% LOG) log(df[[xv]]) else df[[xv]]
  df$Y  <- as.numeric(scale(ty)); df$Xz <- as.numeric(scale(tx))
  df$xw <- df$Xz - ave(df$Xz, df$country)
  m <- suppressWarnings(tryCatch(lmer(Y ~ xw + (xw|country), df, REML=TRUE), error=function(e) NULL))
  het <- NA_real_; rs <- "ri+rs"
  if (is.null(m) || isSingular(m)){ m <- lmer(Y ~ xw + (1|country), df, REML=TRUE); rs <- "ri" }
  else het <- sqrt(VarCorr(m)$country["xw","xw"])
  co <- summary(m)$coefficients; e <- co["xw","Estimate"]; s <- co["xw","Std. Error"]
  data.frame(outcome=yv, predictor=xv, beta=round(e,3), lo=round(e-1.96*s,3), hi=round(e+1.96*s,3),
             het_sd=round(het,2), model=rs, n=nrow(df), n_ctry=length(unique(df$country)),
             sig=ifelse(abs(e)>1.96*s,"*",""))
}

cat("=== Q1: which 2 of {steepness, peak_height, AUC} share most info (within-country std slope) ===\n")
q1 <- rbind(fit_pp(d,"peak_height","auc"), fit_pp(d,"peak_height","steepness"), fit_pp(d,"auc","steepness"))
print(q1, row.names=FALSE)

cat("\n=== Q2: what does onset_week predict (within-country) ===\n")
q2 <- rbind(fit_pp(d,"peak_height","onset_week"), fit_pp(d,"auc","onset_week"),
            fit_pp(d,"peak_week","onset_week"), fit_pp(d,"steepness","onset_week"))
print(q2, row.names=FALSE)

cat("\n=== Q3: what does 65+ vaccination coverage predict (within-country, pre-COVID) ===\n")
dv <- d[is.finite(d$vax_cov_65), ]
q3 <- rbind(fit_pp(dv,"peak_height","vax_cov_65"), fit_pp(dv,"auc","vax_cov_65"),
            fit_pp(dv,"onset_week","vax_cov_65"), fit_pp(dv,"peak_week","vax_cov_65"),
            fit_pp(dv,"steepness","vax_cov_65"))
print(q3, row.names=FALSE)

# forest plot of Q2 + Q3
ff <- rbind(cbind(q2, Q="Q2: onset_week predicts ->"), cbind(q3, Q="Q3: 65+ coverage predicts ->")) %>%
  mutate(lab=paste0(outcome))
p <- ggplot(ff, aes(beta, lab)) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi), color="#1b9e77") +
  facet_wrap(~Q, scales="free_y") + labs(x="within-country standardised slope (95% CI)", y=NULL,
  title="Partial-pooling within-country effects (lme4)") + theme_minimal(base_size=11)
ggsave("output/hierarchical_effects.png", p, width=11, height=4.5, dpi=110)
cat("\nfigure -> output/hierarchical_effects.png\n")
