# precovid_predict_postcovid.R
#
# Fit a PRE-COVID within-country Bayesian model of burden and use it to PREDICT the post-COVID
# seasons 2023/24 and 2024/25, out of sample. Extends bayes_precovid_ve_subtype.R:
#   predictors = dominant SUBTYPE + PROTECTION (VE x 65+ coverage) + PRIOR-SEASON AUC (last season's
#   burden, a country x season predictor). Outcome for the cross-validation: log(AUC).
#
# Prior-season AUC needs the COVID seasons that the committed panel excludes (2023/24's prior is the
# 2022/23 season). We therefore rebuild the ILI+ panel WITHOUT the COVID exclusion via the shared
# stitch (code/01_main_supporting/stitch_iliplus.R -- identical rules to build_slim_panel.R) and
# compute descriptive-method AUC for every season; an assertion below CHECKS (not just asserts in
# prose) that the rebuild reproduces the committed panel's AUCs exactly. Deliberate tension, on the
# record: 2022/23 is excluded as an OUTCOME everywhere (COVID-disrupted wave shape) yet used here as
# a PREDICTOR (its realised burden is what the autoregressive term needs; predictor availability !=
# outcome validity -- see decisions.md).
#
# Cross-validation design: fit on pre-COVID only; predict 2023/24 using 2022/23 AUC as the prior,
# and 2024/25 using 2023/24 AUC as the prior; the country random intercept (estimated pre-COVID)
# carries each country's reporting scale. The TEST set is every train country with a valid prior --
# NOT additionally conditioned on post-COVID coverage availability (protection is no predictor out
# of sample, so requiring coverage would only shrink the test set arbitrarily). VE against the
# dominant subtype is derived from the provenance CSV via analysis_helpers.R::ve_vs_dominant()
# (pre-COVID 14.4/32.9/25.7/45/37.5; test 2023/24 (H1N1) = 52, 2024/25 (H1N1) = 30 end-of-season).
# Run from repo root:  Rscript code/05_analysis/precovid_predict_postcovid.R
suppressMessages({library(dplyr); library(lme4); library(ggplot2)}); set.seed(1)
source("code/02_settings/settings_version0.R"); params <- settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_descriptive.R")
source("code/01_main_supporting/methods_registry.R")
source("code/05_analysis/analysis_helpers.R")

# ---- |-rebuild the ILI+ panel for ALL seasons (shared stitch, minus the COVID exclusion) ----
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R")
models_in <- readRDS("output/models_in.rds")
full <- stitch_iliplus_panel(models_in, exclude_covid=FALSE, min_wk=15)
dir.create("output", showWarnings=FALSE); write.csv(full, "output/slim_flu_iliplus_full.csv", row.names=FALSE)

# ---- |-descriptive AUC etc. for every country-season, then prior (t-1) AUC ----
auc_by_season <- do.call(rbind, lapply(sort(unique(full$country_short)), function(cc){
  sl <- load_flu_iliplus_slim(cc, path="output/slim_flu_iliplus_full.csv")
  s  <- summarise_method_fit(run_method("descriptive", sl, params))
  s[, c("country","season","auc","peak_height","onset_week")] })) %>%
  mutate(syr=as.integer(substr(season,1,4)))
# the rebuild must agree with the committed panel wherever both cover a season (same stitch rules)
if (file.exists("output/descriptors.csv")){
  chk <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE) %>%
    inner_join(auc_by_season, by=c("country","season"), suffix=c("_slim","_full"))
  stopifnot(nrow(chk) > 100, max(abs(log(chk$auc_slim) - log(chk$auc_full))) < 1e-8)
  cat(sprintf("rebuild check: %d country-seasons match the committed panel's AUC exactly\n", nrow(chk)))
}
prior <- auc_by_season %>% transmute(country, syr_next=syr+1, prior_lauc=log(auc))
auc_by_season <- auc_by_season %>% left_join(prior, by=c("country"="country","syr"="syr_next"))

# ---- |-attach subtype, VE-against-dominant, 65+ coverage -> protection ----
sub <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE) %>% transmute(season, dominant)
ve_by_season <- ve_vs_dominant(sub)                                  # traceable to the provenance CSV rows
covpre <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>% transmute(country, season, coverage=vax_cov_65)
covpost<- read.csv("data/external/vaccination_coverage_65plus_postcovid.csv", stringsAsFactors=FALSE) %>%
  filter(panel_country %in% c(TRUE,"TRUE"), !grepl("of invited", age_band)) %>%   # drop NL (60+ of-invited: not comparable)
  transmute(country=country_short, season, coverage=coverage_pct/100)             # percent -> FRACTION (pre-COVID scale)
coverage <- bind_rows(covpre, covpost) %>% filter(is.finite(coverage)) %>% distinct(country, season, .keep_all=TRUE)
stopifnot(max(coverage$coverage) <= 1)   # one scale: fractions (a pct/fraction mix once put test protection ~100x train)

d <- auc_by_season %>% inner_join(sub, by="season") %>% left_join(coverage, by=c("country","season")) %>%
  mutate(ve=ve_by_season[season], protection=ve*coverage/100,        # effectively-protected fraction (VE% x coverage-fraction)
         dominant=factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")), lauc=log(auc))

pre_seasons  <- c("2015/2016","2016/2017","2017/2018","2018/2019")            # 2014/15 has no prior -> excluded
test_seasons <- c("2023/2024","2024/2025")
train <- d %>% filter(season %in% pre_seasons, is.finite(protection), is.finite(prior_lauc))
test  <- d %>% filter(season %in% test_seasons, is.finite(prior_lauc), country %in% unique(train$country))
cat(sprintf("TRAIN pre-COVID: %d country-seasons, %d countries, %d seasons\n", nrow(train), n_distinct(train$country), n_distinct(train$season)))
cat(sprintf("TEST post-COVID: %d country-seasons (%s)\n", nrow(test), paste(table(test$season), names(table(test$season)), collapse=" ")))

# ---- |-within-country deviations, standardised on TRAIN (apply same to TEST) ----
# Demean prior-AUC and protection by the country's TRAINING mean: this breaks the prior-AUC / country-
# intercept ridge (both otherwise encode the country's reporting scale), so the intercept carries the
# LEVEL and the deviation carries the within-country persistence -> stable out-of-sample predictions.
train_means <- train %>% group_by(country) %>% summarise(mean_prior_lauc=mean(prior_lauc), mean_protection=mean(protection), .groups="drop")
dev <- function(df){ df<-left_join(df, train_means, by="country"); df$prior_d<-df$prior_lauc-df$mean_prior_lauc; df$prot_d<-df$protection-df$mean_protection; df }
train<-dev(train); test<-dev(test)
sd_prot<-sd(train$prot_d); sd_prior<-sd(train$prior_d)
train$protection_z<-train$prot_d/sd_prot; train$prior_z<-train$prior_d/sd_prior
test$prior_z <- test$prior_d/sd_prior                                # protection is not used out of sample (see header)

# ---- |-Gibbs varying-intercept model: y = X beta + alpha_country + eps, alpha ~ N(mu_a, t2) ----
# X has NO intercept column -- the country intercept carries each country's full (reporting-scale) level,
# so there is no fixed-intercept/random-intercept ridge. Stores beta, alpha (u), s2 for prediction.
gibbs <- function(y, X, g, n_iter=9000, n_burn=4000, chains=3){
  n<-length(y); p<-ncol(X); G<-max(g); XtX<-crossprod(X); Bm<-Um<-Sm<-NULL
  for (ch in 1:chains){ beta<-rep(0,p); u<-rep(mean(y),G); mua<-mean(y); s2<-var(y); t2<-var(y)/2
    nb<-n_iter-n_burn; B<-matrix(NA,nb,p); U<-matrix(NA,nb,G); S<-numeric(nb)
    for (it in 1:n_iter){
      V<-chol2inv(chol(XtX/s2+diag(1/100,p))); m<-V%*%(crossprod(X,y-u[g])/s2)
      beta<-as.numeric(m+t(chol(V))%*%rnorm(p)); e<-as.numeric(y-X%*%beta)
      for (gi in 1:G){ idx<-which(g==gi); vc<-1/(length(idx)/s2+1/t2); u[gi]<-rnorm(1, vc*(sum(e[idx])/s2+mua/t2), sqrt(vc)) }
      vmu<-1/(G/t2+1/100); mua<-rnorm(1, vmu*sum(u)/t2, sqrt(vmu))
      r<-as.numeric(y-X%*%beta-u[g]); s2<-1/rgamma(1,0.01+n/2,0.01+sum(r^2)/2); t2<-1/rgamma(1,0.01+G/2,0.01+sum((u-mua)^2)/2)
      if (it>n_burn){ k<-it-n_burn; B[k,]<-beta; U[k,]<-u; S[k]<-s2 } }
    Bm<-rbind(Bm,B); Um<-rbind(Um,U); Sm<-c(Sm,S) }
  list(beta=Bm, u=Um, s2=Sm)
}

country_levels <- levels(factor(train$country)); g <- as.integer(factor(train$country, levels=country_levels))
X_train  <- cbind(model.matrix(~ dominant, train)[,-1,drop=FALSE], prot=train$protection_z, prior=train$prior_z)  # [H3N2,B,prot,prior] (no intercept; country RE carries level)

# ---- |-whisker plot: subtype + protection + prior-AUC across the three descriptors (SD units) ----
outcome_labels <- c(lauc="AUC (log)", onset="onset week", peak="peak height (log)")
outcome_data <- list(lauc=log(train$auc), onset=train$onset_week, peak=log(train$peak_height))
whisker <- do.call(rbind, lapply(names(outcome_labels), function(k){
  y <- as.numeric(scale(outcome_data[[k]])); fit <- gibbs(y, X_train, g); A <- fit$beta
  draws <- list(`subtype: H3N2-H1N1`=A[,1], `subtype: B-H1N1`=A[,2], `subtype: B-H3N2`=A[,2]-A[,1],
                `slope: protection`=A[,3], `slope: prior-AUC`=A[,4])
  do.call(rbind, lapply(names(draws), function(nm){ q<-q95(draws[[nm]])
    data.frame(outcome=outcome_labels[k], term=nm, est=round(q[1],2), lo=round(q[2],2), hi=round(q[3],2),
               excl0=ifelse(q[2]>0|q[3]<0,"*","")) })) }))
write.csv(whisker, "output/precovid3_whisker.csv", row.names=FALSE)
cat("\n=== pre-COVID: subtype + protection + prior-AUC (within-country, SD units, 95% CrI) ===\n"); print(whisker, row.names=FALSE)
whisker$term<-factor(whisker$term, levels=rev(c("subtype: H3N2-H1N1","subtype: B-H1N1","subtype: B-H3N2","slope: protection","slope: prior-AUC")))
ggsave("output/precovid3_whisker.png",
  ggplot(whisker, aes(est,term,color=excl0=="*"))+geom_vline(xintercept=0,color="grey60")+
    geom_pointrange(aes(xmin=lo,xmax=hi))+facet_wrap(~outcome,nrow=1)+
    scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"),guide="none")+
    labs(title="Pre-COVID within-country: subtype + protection (VE x coverage) + prior-season AUC",
         subtitle="Country random intercepts; SD units, 95% CrI. Prior-AUC is a country x season predictor (unlike season-level subtype/VE).",
         x="effect (SD units)", y=NULL)+theme_minimal(base_size=10), width=12, height=4, dpi=110)

# ---- |-hierarchical location-scale prediction models: COUNTRY-SPECIFIC residual variance ----
# The season-to-season variability UNEXPLAINED by the predictors is country-specific: sigma_c^2 ~
# InvGamma(a0, b0), the scale b0 estimated across countries (partial pooling), so a volatile country gets
# a WIDER predictive band and a steady one a narrower band -- not a single pooled residual. Two Gibbs models:
#   BASELINE = country intercept alpha_c + sigma_c (knows only the country);
#   FULL     = pooled mean (dominant subtype + prior-season AUC, which carries the reporting scale) + sigma_c.
# Each test season is predicted with ITS country's sigma_c. a0=2 sets moderate variance pooling.
# (PROTECTION is excluded out of sample: post-COVID coverage exists for only ~15 country-seasons, and the
#  VE regime differs; prior-AUC is on the same log scale train/test, so no extrapolation.)
gibbs_ls <- function(y, g, X=NULL, n_iter=9000, n_burn=4000, chains=3, a0=2){
  n<-length(y); G<-max(g); ints<-is.null(X); p<-if (ints) 1L else ncol(X); A<-Sg<-B<-NULL
  for (ch in 1:chains){
    alpha<-rep(mean(y),G); mua<-mean(y); ta2<-var(y); s2<-rep(var(y),G); b0<-var(y); beta<-rep(0,p)
    nb<-n_iter-n_burn; Ac<-matrix(NA,nb,G); Sc<-matrix(NA,nb,G); Bc<-matrix(NA,nb,p)
    for (it in 1:n_iter){
      if (ints){
        for (gi in 1:G){ idx<-which(g==gi); vc<-1/(length(idx)/s2[gi]+1/ta2)
          alpha[gi]<-rnorm(1, vc*(sum(y[idx])/s2[gi]+mua/ta2), sqrt(vc)) }
        vmu<-1/(G/ta2+1/100); mua<-rnorm(1, vmu*sum(alpha)/ta2, sqrt(vmu))
        ta2<-1/rgamma(1, 0.01+G/2, 0.01+sum((alpha-mua)^2)/2); r<-y-alpha[g]
      } else {
        w<-1/s2[g]; V<-chol2inv(chol(crossprod(X*sqrt(w))+diag(1/100,p))); m<-V%*%crossprod(X, w*y)
        beta<-as.numeric(m+t(chol(V))%*%rnorm(p)); r<-y-as.numeric(X%*%beta)
      }
      for (gi in 1:G){ idx<-which(g==gi); s2[gi]<-1/rgamma(1, a0+length(idx)/2, b0+sum(r[idx]^2)/2) }
      b0<-rgamma(1, 1+G*a0, 1+sum(1/s2))
      if (it>n_burn){ k<-it-n_burn; Ac[k,]<-alpha; Sc[k,]<-sqrt(s2); Bc[k,]<-beta } }
    A<-rbind(A,Ac); Sg<-rbind(Sg,Sc); B<-rbind(B,Bc) }
  list(alpha=A, sigma=Sg, beta=B)
}
X_full <- cbind(1, model.matrix(~ dominant, train)[,-1,drop=FALSE], train$prior_lauc)   # [1, H3N2, B, prior_lauc]
fit_full <- gibbs_ls(train$lauc, g, X=X_full)          # FULL: pooled mean + country-specific sigma
fit_base <- gibbs_ls(train$lauc, g)                # BASELINE: country intercept + country-specific sigma
cat(sprintf("\nfull-model coefficients (intercept, H3N2, B, prior-AUC): %s\n", paste(round(colMeans(fit_full$beta),3), collapse=", ")))
cat(sprintf("country-specific baseline residual SD (season-to-season, log): %.2f - %.2f (median %.2f)\n",
            min(colMeans(fit_base$sigma)), max(colMeans(fit_base$sigma)), median(colMeans(fit_base$sigma))))
g_test <- match(test$country, country_levels)
X_test <- cbind(1, as.numeric(test$dominant=="A(H3N2)"), as.numeric(test$dominant=="B"), test$prior_lauc)
predict_ls <- function(mu, sig){ obs<-mu+rnorm(length(mu),0,sig); c(fit=mean(mu), lo=quantile(obs,.025), hi=quantile(obs,.975)) }
pred_full <- sapply(seq_len(nrow(test)), function(i) predict_ls(as.numeric(X_test[i,]%*%t(fit_full$beta)), fit_full$sigma[,g_test[i]]))
test$pred_lauc<-pred_full[1,]; test$pi_lo<-pred_full[2,]; test$pi_hi<-pred_full[3,]
pred_base <- sapply(seq_len(nrow(test)), function(i) predict_ls(fit_base$alpha[,g_test[i]], fit_base$sigma[,g_test[i]]))
test$base_pred<-pred_base[1,]; test$base_lo<-pred_base[2,]; test$base_hi<-pred_base[3,]
test$base_sd <- colMeans(fit_base$sigma)[g_test]
test$in_pi <- test$lauc>=test$pi_lo & test$lauc<=test$pi_hi
rmse  <- sqrt(mean((test$lauc-test$pred_lauc)^2)); r_all <- cor(test$lauc, test$pred_lauc)
rmse0 <- sqrt(mean((test$lauc-test$base_pred)^2)); r0 <- cor(test$lauc, test$base_pred)
rmse_pers <- sqrt(mean((test$lauc-test$prior_lauc)^2))   # PERSISTENCE comparator: predict this season = last season's log AUC
# within-country signal: correlate the two-season deviations for countries with BOTH test seasons only
# (single-season countries would enter as exact (0,0) points and mechanically inflate the correlation)
paired <- test %>% add_count(country, name="n_test") %>% filter(n_test==2) %>% group_by(country) %>%
  mutate(obs_dev=lauc-mean(lauc), pred_dev=pred_lauc-mean(pred_lauc)) %>% ungroup()
r_within <- if (nrow(paired) >= 4) cor(paired$obs_dev, paired$pred_dev) else NA_real_
cat(sprintf("\nCROSS-VAL (log AUC): FULL RMSE=%.2f cor=%.2f (within-dev %.2f over %d paired countries, %.0f%% in PI)\n                     BASELINE RMSE=%.2f cor=%.2f | PERSISTENCE (carry last season forward) RMSE=%.2f\n",
            rmse, r_all, r_within, n_distinct(paired$country), 100*mean(test$in_pi), rmse0, r0, rmse_pers))
print(test %>% transmute(country, season, real=round(lauc,2), full=round(pred_lauc,2),
                        full_pi=sprintf("[%.2f,%.2f]",pi_lo,pi_hi), base_sd=round(base_sd,2), in_pi) %>% as.data.frame(), row.names=FALSE)
write.csv(test %>% select(country,season,lauc,prior_lauc,base_pred,base_lo,base_hi,base_sd,pred_lauc,pi_lo,pi_hi,in_pi), "output/precovid_crossval.csv", row.names=FALSE)

lims <- range(c(test$lauc, test$pred_lauc))
ggsave("output/precovid_crossval.png",
  ggplot(test, aes(pred_lauc, lauc, color=season))+
    geom_abline(slope=1,intercept=0,color="grey60",linetype="dashed")+
    geom_errorbar(aes(xmin=pi_lo,xmax=pi_hi),width=0,alpha=0.4,orientation="y")+ geom_point(size=2.4)+
    geom_text(aes(label=country),size=2.6,vjust=-0.7,show.legend=FALSE)+
    scale_color_manual(values=c("2023/2024"="#1b9e77","2024/2025"="#d95f02"))+
    coord_equal(xlim=lims,ylim=lims)+
    labs(title="Out-of-sample: pre-COVID model predicting post-COVID burden (log AUC)",
         subtitle=sprintf("Predict 2023/24 from 2022/23 prior, 2024/25 from 2023/24 prior. RMSE=%.2f, cor=%.2f, %.0f%% in 95%% PI. Bars = 95%% predictive interval.",
                          rmse, r_all, 100*mean(test$in_pi)),
         x="predicted log(AUC)", y="observed log(AUC)")+theme_minimal(base_size=11), width=8, height=8, dpi=110)

# ---- |-relative view: AUC vs each country's pre-COVID mean, with COUNTRY-SPECIFIC baseline bands ----
# Center each test season on its country's pre-COVID mean (= the baseline point). The baseline predicts 0
# (a "typical season"); its blue 95% band width is now COUNTRY-SPECIFIC (each country's own residual SD),
# so volatile countries get taller bands. The full-model and observed DEVIATIONS are what the eye reads.
# y is on the log scale but tick-labelled as fold-change (1x = pre-COVID average; higher = a bigger season).
# Countries with BOTH test seasons are grouped first and joined by a line; single-season countries follow
# after a dotted gap and stand alone.
rel <- test %>%
  mutate(obs_r=lauc-base_pred, full_r=pred_lauc-base_pred, full_lo=pi_lo-base_pred, full_hi=pi_hi-base_pred,
         seas=ifelse(season=="2023/2024","'23/24","'24/25")) %>%
  add_count(country, name="n_seasons") %>% mutate(paired=n_seasons==2)
country_order <- rel %>% distinct(country, paired) %>% arrange(desc(paired), country) %>% pull(country)   # pairs first, then singles
rel <- rel %>% mutate(country=factor(country, levels=country_order)) %>% arrange(country, season) %>%
  mutate(pos = row_number() + ifelse(paired, 0, 0.8))                    # gap before the singles block
gapx <- if (any(!rel$paired)) min(rel$pos[!rel$paired]) - 0.4 else NA_real_
baseline_bands <- rel %>% group_by(country) %>%                                       # per-country baseline band (its own SD)
  summarise(xmin=min(pos)-0.5, xmax=max(pos)+0.5, blo=(base_lo-base_pred)[1], bhi=(base_hi-base_pred)[1], .groups="drop") %>%
  mutate(shade=rep(c("a","b"), length.out=n()))
y_breaks <- log(c(0.125,0.25,0.5,1,1.5))
ylim <- c(min(rel$obs_r, rel$full_lo, baseline_bands$blo)-0.1, max(rel$obs_r, rel$full_hi, baseline_bands$bhi, 0.15)+0.1)
p <- ggplot(rel)+
  geom_rect(data=baseline_bands, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=shade), alpha=0.6, inherit.aes=FALSE)+
  scale_fill_manual(values=c(a="grey96", b="white"), guide="none")+
  geom_rect(data=baseline_bands, aes(xmin=xmin,xmax=xmax,ymin=blo,ymax=bhi), fill="#8da0cb", alpha=0.32, inherit.aes=FALSE)+  # country-specific baseline 95% PI
  geom_hline(yintercept=0, color="grey55", linewidth=0.4)+                            # baseline point prediction = 'typical season'
  geom_line(aes(pos, obs_r, group=country), color="grey45", linewidth=0.5)+          # joins a country's 2 seasons
  geom_linerange(aes(pos, ymin=full_lo, ymax=full_hi), color="#d95f02", linewidth=0.8)+
  geom_point(aes(pos, full_r), color="#d95f02", size=2.4)+
  geom_point(aes(pos, obs_r), shape=18, size=3.4, color="black")+
  scale_x_continuous(breaks=rel$pos, labels=paste0(rel$country," ",rel$seas))+
  scale_y_continuous(breaks=y_breaks, labels=c("0.12x","0.25x","0.5x","1x","1.5x"))+
  coord_cartesian(ylim=ylim)+
  labs(title="Post-COVID burden relative to each country's pre-COVID average (country scale removed)",
       subtitle=sprintf("diamond = observed | orange = full model (95%% PI) | blue band = baseline 95%% PI, COUNTRY-SPECIFIC width (its prediction = the 1x line)\ngrey line joins a country's two seasons; single-season countries after the gap. RMSE %.2f baseline / %.2f persistence / %.2f full (n=%d test country-seasons).", rmse0, rmse_pers, rmse, nrow(test)),
       x=NULL, y="AUC relative to the country's pre-COVID mean (fold-change)")+
  theme_minimal(base_size=11)+
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        plot.subtitle=element_text(size=9))
if (!is.na(gapx)) p <- p + geom_vline(xintercept=gapx, color="grey70", linetype="dotted")
ggsave("output/precovid_crossval_compare.png", p, width=11, height=6, dpi=110)
write.csv(rel %>% select(country,season,obs_rel=obs_r,full_pred_rel=full_r,full_lo,full_hi), "output/precovid_crossval_compare.csv", row.names=FALSE)
cat("\nfigures -> output/precovid3_whisker.png ; output/precovid_crossval.png ; output/precovid_crossval_compare.png\n")
