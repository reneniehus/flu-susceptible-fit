# precovid_predict_postcovid.R
#
# Fit a PRE-COVID within-country Bayesian model of burden and use it to PREDICT the two post-COVID
# seasons that have 65+ coverage (2023/24, 2024/25), out of sample. Extends bayes_precovid_ve_subtype.R:
#   predictors = dominant SUBTYPE + PROTECTION (VE x 65+ coverage) + PRIOR-SEASON AUC (last season's
#   burden, a country x season predictor). Outcome for the cross-validation: log(AUC).
#
# Prior-season AUC needs the COVID seasons that the committed panel excludes (2023/24's prior is the
# 2022/23 season). We therefore rebuild the ILI+ panel here WITHOUT the COVID exclusion (identical
# reconstruction to build_slim_panel.R otherwise; verified to reproduce the committed 2023/24 AUC
# exactly) and compute descriptive-method AUC for every season.
#
# Cross-validation design (as requested): fit on pre-COVID only; predict 2023/24 using 2022/23 AUC as the
# prior, and 2024/25 using 2023/24 AUC as the prior; the country random intercept (estimated pre-COVID)
# carries each country's reporting scale. VE against the dominant subtype per season (see
# data/external/vaccine_effectiveness.csv): pre-COVID 14.4/32.9/25.7/45/38; test 2023/24 (H1N1)=52,
# 2024/25 (B)=58 (interim). Run from repo root:  Rscript code/05_analysis/precovid_predict_postcovid.R
suppressMessages({library(dplyr); library(lme4); library(ggplot2)}); set.seed(1)
source("code/02_settings/settings_version0.R"); params <- settings()
source("code/01_main_supporting/sir_core.R")
for (f in c("method_sir_deterministic","method_sir_ekf","method_descriptive")) source(paste0("code/01_main_supporting/methods/",f,".R"))
source("code/01_main_supporting/methods_registry.R")

# ---- |-rebuild the ILI+ panel for ALL seasons (mirrors build_slim_panel.R minus the COVID exclusion) ----
mi <- readRDS("output/models_in.rds")
nonsent<-c("MT","IS","HR","RO","LV","FI"); x1000<-c("CY","LU","MT"); resp_only<-c("NO","ES"); erviss_only<-c("SK","LV")
ip <- mi$data_timeseries_long %>% filter(indicator=="ili_plus", pathogen=="Influenza", agegroup=="age_total") %>%
  select(stream, country_short, season, date, season_week, value)
erv <- ip %>% filter(stream %in% c("ili_plus_sentinel","ili_plus_nonsentinel")) %>%
  mutate(pick=ifelse(country_short %in% nonsent,"ili_plus_nonsentinel","ili_plus_sentinel")) %>% filter(stream==pick) %>%
  mutate(value=value*ifelse(country_short %in% x1000,1000,1)) %>% transmute(country_short,season,date,season_week,erv=value)
rsp <- ip %>% filter(stream=="ili_plus_respicompass") %>% transmute(country_short,season,date,season_week,rsp=value)
m <- full_join(erv,rsp,by=c("country_short","season","date","season_week"))
fac <- m %>% filter(season=="2023/2024",is.finite(erv),erv>0,is.finite(rsp),rsp>0) %>% group_by(country_short) %>% summarise(factor=median(rsp/erv),.groups="drop")
m <- m %>% left_join(fac,by="country_short") %>% mutate(factor=ifelse(is.na(factor),1,factor),
  value=case_when(country_short %in% resp_only~rsp, country_short %in% erviss_only~erv, !is.na(rsp)~rsp, TRUE~erv*factor),
  source=case_when(country_short %in% resp_only~"RespiCompass", country_short %in% erviss_only~"ERVISS", !is.na(rsp)~"RespiCompass", TRUE~"ERVISS")) %>%
  filter(is.finite(value))
full <- m %>% group_by(country_short,season) %>% filter(sum(is.finite(value)&value>0)>=15) %>%
  group_modify(function(df,key){ w_hi<-max(df$season_week[is.finite(df$value)])
    tibble(season_week=1:w_hi) %>% left_join(df %>% select(season_week,date,value),by="season_week") %>%
      transmute(week=season_week,season_week,date,value,source=df$source[which(!is.na(df$source))][1]) }) %>% ungroup()
dir.create("output", showWarnings=FALSE); write.csv(full, "output/slim_flu_iliplus_full.csv", row.names=FALSE)

# ---- |-descriptive AUC etc. for every country-season, then prior (t-1) AUC ----
au <- do.call(rbind, lapply(sort(unique(full$country_short)), function(cc){
  sl <- load_flu_iliplus_slim(cc, path="output/slim_flu_iliplus_full.csv")
  s  <- summarise_method_fit(run_method("descriptive", sl, params))
  s[, c("country","season","auc","peak_height","onset_week")] })) %>%
  mutate(syr=as.integer(substr(season,1,4)))
prior <- au %>% transmute(country, syr_next=syr+1, prior_lauc=log(auc))
au <- au %>% left_join(prior, by=c("country"="country","syr"="syr_next"))

# ---- |-attach subtype, VE-against-dominant, 65+ coverage -> protection ----
sub <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE) %>% transmute(season, dominant)
ve_by_season <- c("2014/2015"=14.4,"2015/2016"=32.9,"2016/2017"=25.7,"2017/2018"=45,"2018/2019"=38,
                  "2023/2024"=52,"2024/2025"=58)                     # VE vs dominant subtype (see header)
covpre <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>% transmute(country, season, cov=vax_cov_65)
covpost<- read.csv("data/external/vaccination_coverage_65plus_postcovid.csv", stringsAsFactors=FALSE) %>%
  filter(panel_country %in% c(TRUE,"TRUE"), !grepl("of invited", age_band)) %>%   # drop NL (60+ of-invited: not comparable)
  transmute(country=country_short, season, cov=coverage_pct)
cov <- bind_rows(covpre, covpost) %>% filter(is.finite(cov)) %>% distinct(country, season, .keep_all=TRUE)

d <- au %>% inner_join(sub, by="season") %>% left_join(cov, by=c("country","season")) %>%
  mutate(ve=ve_by_season[season], protection=ve*cov/100,
         dominant=factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")), lauc=log(auc))

pre  <- c("2015/2016","2016/2017","2017/2018","2018/2019")            # 2014/15 has no prior -> excluded
test <- c("2023/2024","2024/2025")
train <- d %>% filter(season %in% pre, is.finite(protection), is.finite(prior_lauc))
tst   <- d %>% filter(season %in% test, is.finite(protection), is.finite(prior_lauc), country %in% unique(train$country))
cat(sprintf("TRAIN pre-COVID: %d country-seasons, %d countries, %d seasons\n", nrow(train), n_distinct(train$country), n_distinct(train$season)))
cat(sprintf("TEST post-COVID: %d country-seasons (%s)\n", nrow(tst), paste(table(tst$season), names(table(tst$season)), collapse=" ")))

# ---- |-within-country deviations, standardised on TRAIN (apply same to TEST) ----
# Demean prior-AUC and protection by the country's TRAINING mean: this breaks the prior-AUC / country-
# intercept ridge (both otherwise encode the country's reporting scale), so the intercept carries the
# LEVEL and the deviation carries the within-country persistence -> stable out-of-sample predictions.
cm <- train %>% group_by(country) %>% summarise(mpc=mean(prior_lauc), mpp=mean(protection), .groups="drop")
dev <- function(df){ df<-left_join(df, cm, by="country"); df$prior_d<-df$prior_lauc-df$mpc; df$prot_d<-df$protection-df$mpp; df }
train<-dev(train); tst<-dev(tst)
sp<-sd(train$prot_d); sl_<-sd(train$prior_d)
train$protection_z<-train$prot_d/sp; train$prior_z<-train$prior_d/sl_
tst$protection_z  <-tst$prot_d/sp;   tst$prior_z  <-tst$prior_d/sl_

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

clev <- levels(factor(train$country)); g <- as.integer(factor(train$country, levels=clev))
Xtr  <- cbind(model.matrix(~ dominant, train)[,-1,drop=FALSE], prot=train$protection_z, prior=train$prior_z)  # [H3N2,B,prot,prior] (no intercept; country RE carries level)

# ---- |-whisker plot: subtype + protection + prior-AUC across the three descriptors (SD units) ----
q95<-function(x) quantile(x,c(.5,.025,.975))
mk <- c(lauc="AUC (log)", onset="onset week", peak="peak height (log)")
src <- list(lauc=log(train$auc), onset=train$onset_week, peak=log(train$peak_height))
wres <- do.call(rbind, lapply(names(mk), function(k){
  y <- as.numeric(scale(src[[k]])); fit <- gibbs(y, Xtr, g); A <- fit$beta
  draws <- list(`subtype: H3N2-H1N1`=A[,1], `subtype: B-H1N1`=A[,2], `subtype: B-H3N2`=A[,2]-A[,1],
                `slope: protection`=A[,3], `slope: prior-AUC`=A[,4])
  do.call(rbind, lapply(names(draws), function(nm){ q<-q95(draws[[nm]])
    data.frame(outcome=mk[k], term=nm, est=round(q[1],2), lo=round(q[2],2), hi=round(q[3],2),
               excl0=ifelse(q[2]>0|q[3]<0,"*","")) })) }))
write.csv(wres, "output/precovid3_whisker.csv", row.names=FALSE)
cat("\n=== pre-COVID: subtype + protection + prior-AUC (within-country, SD units, 95% CrI) ===\n"); print(wres, row.names=FALSE)
wres$term<-factor(wres$term, levels=rev(c("subtype: H3N2-H1N1","subtype: B-H1N1","subtype: B-H3N2","slope: protection","slope: prior-AUC")))
ggsave("output/precovid3_whisker.png",
  ggplot(wres, aes(est,term,color=excl0=="*"))+geom_vline(xintercept=0,color="grey60")+
    geom_pointrange(aes(xmin=lo,xmax=hi))+facet_wrap(~outcome,nrow=1)+
    scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"),guide="none")+
    labs(title="Pre-COVID within-country: subtype + protection (VE x coverage) + prior-season AUC",
         subtitle="Country random intercepts; SD units, 95% CrI. Prior-AUC is a country x season predictor (unlike season-level subtype/VE).",
         x="effect (SD units)", y=NULL)+theme_minimal(base_size=10), width=12, height=4, dpi=110)

# ---- |-cross-validation: pooled AR (prior log-AUC carries the country scale); Gaussian predictive ----
# For PREDICTION we let prior-season log(AUC) carry each country's reporting scale directly (AR-style)
# rather than a country random intercept -- the two are collinear (both encode scale) and the RE model
# cannot extrapolate out of sample. Fit a pooled linear model; the 95% prediction interval is the Gaussian
# posterior predictive (= Bayesian with weak priors). prior_lauc is on the natural log-AUC scale, identical
# for train and test, so there is no extrapolation.
# NOTE: PROTECTION is dropped from the OUT-OF-SAMPLE predictor set. VE jumped pre->post COVID (dominant-
# subtype VE 14-45% pre vs 52-58% test), so post-COVID protection lies outside the training range and does
# not transfer (its within-country effect is ~0 anyway; whisker above). Predict from the transferable
# predictors: dominant subtype + prior-season AUC (which carries each country's reporting scale).
mf <- lm(lauc ~ dominant + prior_lauc, data=train)
pr <- predict(mf, newdata=tst, interval="prediction", level=0.95)
tst$pred_lauc <- pr[,"fit"]; tst$pi_lo <- pr[,"lwr"]; tst$pi_hi <- pr[,"upr"]
cat("\npredictive model (log AUC ~ subtype + prior-AUC): coefficients\n"); print(round(coef(mf),3))
tst$in_pi <- tst$lauc >= tst$pi_lo & tst$lauc <= tst$pi_hi
rmse <- sqrt(mean((tst$lauc-tst$pred_lauc)^2)); r_all <- cor(tst$lauc, tst$pred_lauc)
# within-country deviation cor (the hard part: did we get the season-to-season move right?)
tst <- tst %>% group_by(country) %>% mutate(rc=lauc-mean(lauc), pc=pred_lauc-mean(pred_lauc)) %>% ungroup()
r_within <- if (sum(tst$rc!=0)>2) cor(tst$rc, tst$pc) else NA_real_
cat(sprintf("\nCROSS-VAL (log AUC): RMSE=%.2f | cor(real,pred)=%.2f | within-country dev cor=%.2f | %.0f%% inside 95%% PI\n",
            rmse, r_all, r_within, 100*mean(tst$in_pi)))
print(tst %>% transmute(country, season, real_lauc=round(lauc,2), pred_lauc=round(pred_lauc,2),
                        pi=sprintf("[%.2f, %.2f]",pi_lo,pi_hi), in_pi) %>% as.data.frame(), row.names=FALSE)
write.csv(tst %>% select(country,season,lauc,pred_lauc,pi_lo,pi_hi,in_pi), "output/precovid_crossval.csv", row.names=FALSE)

lims <- range(c(tst$lauc, tst$pred_lauc))
ggsave("output/precovid_crossval.png",
  ggplot(tst, aes(pred_lauc, lauc, color=season))+
    geom_abline(slope=1,intercept=0,color="grey60",linetype="dashed")+
    geom_errorbarh(aes(xmin=pi_lo,xmax=pi_hi),height=0,alpha=0.4)+ geom_point(size=2.4)+
    geom_text(aes(label=country),size=2.6,vjust=-0.7,show.legend=FALSE)+
    scale_color_manual(values=c("2023/2024"="#1b9e77","2024/2025"="#d95f02"))+
    coord_equal(xlim=lims,ylim=lims)+
    labs(title="Out-of-sample: pre-COVID model predicting post-COVID burden (log AUC)",
         subtitle=sprintf("Predict 2023/24 from 2022/23 prior, 2024/25 from 2023/24 prior. RMSE=%.2f, cor=%.2f, %.0f%% in 95%% PI. Bars = 95%% predictive interval.",
                          rmse, r_all, 100*mean(tst$in_pi)),
         x="predicted log(AUC)", y="observed log(AUC)")+theme_minimal(base_size=11), width=8, height=8, dpi=110)

# ---- |-baseline (country-only) vs full model: per-country-season prediction comparison ----
# Baseline knows ONLY the country -> predicts that country's mean log(AUC) (an intercept-only model;
# lm(~country) here; Bayesian partial pooling is near-identical at ~4 training seasons/country). The full
# model adds dominant subtype + prior-season AUC. The gap shows what those buy beyond the country scale.
mf0 <- lm(lauc ~ country, data=train)
pr0 <- predict(mf0, newdata=tst, interval="prediction", level=0.95)
tst$base_pred <- pr0[,"fit"]; tst$base_lo <- pr0[,"lwr"]; tst$base_hi <- pr0[,"upr"]
rmse0 <- sqrt(mean((tst$lauc-tst$base_pred)^2)); r0 <- cor(tst$lauc, tst$base_pred)
cat(sprintf("\nBASELINE (country only): RMSE=%.2f cor=%.2f  vs  FULL (+subtype+prior-AUC): RMSE=%.2f cor=%.2f\n", rmse0, r0, rmse, r_all))

ord <- tst %>% mutate(cs=paste(country, season)) %>% arrange(lauc) %>% pull(cs)
cmp <- bind_rows(
  tst %>% transmute(cs=paste(country,season), model="baseline (country only)",       pred=base_pred, lo=base_lo, hi=base_hi),
  tst %>% transmute(cs=paste(country,season), model="full (+ subtype + prior-AUC)",  pred=pred_lauc, lo=pi_lo,   hi=pi_hi)) %>%
  mutate(cs=factor(cs, levels=ord))
obs_df <- tst %>% transmute(cs=factor(paste(country,season), levels=ord), obs=lauc)
ggsave("output/precovid_crossval_compare.png",
  ggplot(cmp, aes(pred, cs))+
    geom_linerange(aes(xmin=lo, xmax=hi, color=model), position=position_dodge(0.55), linewidth=0.7, alpha=0.55)+
    geom_point(aes(color=model), position=position_dodge(0.55), size=2.3)+
    geom_point(data=obs_df, aes(x=obs, y=cs), shape=18, size=3.2, color="black", inherit.aes=FALSE)+
    scale_color_manual(values=c("baseline (country only)"="grey60","full (+ subtype + prior-AUC)"="#d95f02"))+
    labs(title="Post-COVID burden: country-only baseline vs full model",
         subtitle=sprintf("Diamond = observed; bars = 95%% predictive interval. RMSE %.2f (country only) -> %.2f (full: + subtype + prior-AUC).", rmse0, rmse),
         x="predicted log(AUC)   (diamond = observed)", y=NULL, color=NULL)+
    theme_minimal(base_size=11)+theme(legend.position="top", panel.grid.minor=element_blank()), width=9.5, height=5.5, dpi=110)
write.csv(tst %>% select(country,season,lauc,base_pred,base_lo,base_hi,pred_lauc,pi_lo,pi_hi), "output/precovid_crossval_compare.csv", row.names=FALSE)
cat("\nfigures -> output/precovid3_whisker.png ; output/precovid_crossval.png ; output/precovid_crossval_compare.png\n")
