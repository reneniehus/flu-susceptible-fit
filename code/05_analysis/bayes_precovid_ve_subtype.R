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
# Season-level VE against the DOMINANT subtype (all-ages primary care, end-of-season where published,
# else best-available early/interim; full rows + sources in data/external/vaccine_effectiveness.csv):
#   2014/15 A(H3N2) 14.4  I-MOVE eos            [ES.2016.21.7.30139]
#   2015/16 A(H1N1) 32.9  I-MOVE eos            [PMC6005601]
#   2016/17 A(H3N2) 25.7  I-MOVE early, target  [ES.2017.22.7.30464] (no all-ages eos published)
#   2017/18 B       45    interim B midpoint    [ES.2018.23.9.18-00086] (36-54; Yamagata mismatch)
#   2018/19 A(H1N1) 38    interim A midpoint    [ES.2019.24.1900121]    (32-43; H1N1 the higher end)
# Run from the repo root:  Rscript code/05_analysis/bayes_precovid_ve_subtype.R
suppressMessages({library(dplyr); library(lme4)}); set.seed(1)

pre <- c("2014/2015","2015/2016","2016/2017","2017/2018","2018/2019")
ve_by_season <- c("2014/2015"=14.4, "2015/2016"=32.9, "2016/2017"=25.7, "2017/2018"=45, "2018/2019"=38)

sub <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE) %>% transmute(season, dominant)
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

gibbs <- function(y, X, g, n_iter=6000, n_burn=2000, chains=3){
  n<-length(y); p<-ncol(X); G<-max(g); XtX<-crossprod(X); keep<-list()
  for (ch in 1:chains){
    beta<-rep(0,p); u<-rep(0,G); s2<-var(y); t2<-var(y)/2
    M<-matrix(NA,n_iter-n_burn,p)
    for (it in 1:n_iter){
      V<-chol2inv(chol(XtX/s2 + diag(1/100,p))); m<-V%*%(crossprod(X,y-u[g])/s2)
      beta<-as.numeric(m + t(chol(V))%*%rnorm(p)); e<-as.numeric(y - X%*%beta)
      for (gi in 1:G){ idx<-which(g==gi); vc<-1/(length(idx)/s2 + 1/t2); u[gi]<-rnorm(1, vc*sum(e[idx])/s2, sqrt(vc)) }
      r<-as.numeric(y - X%*%beta - u[g]); s2<-1/rgamma(1,0.01+n/2,0.01+sum(r^2)/2)
      t2<-1/rgamma(1,0.01+G/2,0.01+sum(u^2)/2)
      if (it>n_burn) M[it-n_burn,]<-beta
    }
    keep[[ch]]<-M
  }
  keep
}
rhat1 <- function(chs,col){ L<-nrow(chs[[1]]); cm<-sapply(chs,function(M)mean(M[,col]))
  B<-L*var(cm); W<-mean(sapply(chs,function(M)var(M[,col]))); sqrt(((L-1)/L*W+B/L)/W) }
q95 <- function(x) quantile(x, c(.5,.025,.975))

# ---- |-fit subtype + one continuous predictor, report contrasts + slope per outcome ----
fit_model <- function(dat, xvar, xlabel){
  dd <- dat[is.finite(dat[[xvar]]), ]; g <- as.integer(factor(dd$country))
  outs <- c(lauc="AUC (log)", onset_week="onset week", lpk="peak height (log)")
  res <- list()
  for (o in names(outs)){
    y <- as.numeric(scale(dd[[o]])); xz <- as.numeric(scale(dd[[xvar]]))
    X <- cbind(model.matrix(~ dominant, dd), xz)                    # [1, H3N2, B, predictor]
    ch <- gibbs(y, X, g); A <- do.call(rbind, ch)
    rows <- list(`H3N2 - H1N1`=A[,2], `B - H1N1`=A[,3], `B - H3N2`=A[,3]-A[,2], predictor=A[,4])
    labs <- c(`H3N2 - H1N1`="subtype: H3N2-H1N1", `B - H1N1`="subtype: B-H1N1",
              `B - H3N2`="subtype: B-H3N2", predictor=paste0("slope: ", xlabel))
    for (nm in names(rows)){ qq <- q95(rows[[nm]])
      res[[length(res)+1]] <- data.frame(outcome=outs[o], term=labs[[nm]], est=round(qq[1],2),
        lo=round(qq[2],2), hi=round(qq[3],2), excl0=ifelse(qq[2]>0|qq[3]<0,"*",""),
        rhat=round(rhat1(ch, if(nm=="predictor") 4 else if(nm=="B - H3N2") 3 else which(names(rows)==nm)+1),3)) }
    m <- suppressWarnings(lmer(y ~ dominant + xz + (1|country), dd, REML=TRUE))
    cat(sprintf("  [%-16s | %-14s] lme4 slope(%s)=%.2f | gibbs=%.2f\n", xlabel, outs[o], xvar, fixef(m)["xz"], mean(A[,4])))
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
ps <- m2[grepl("slope:", m2$term), c("outcome","est","lo","hi","excl0")]
print(ps, row.names=FALSE)

# ---- |-forest plot of both models ----
suppressMessages(library(ggplot2))
m1$model <- "Model 1: subtype + VE"; m2$model <- "Model 2: subtype + protection (VE x coverage)"
pd <- rbind(m1, m2); pd$term <- factor(pd$term, levels=rev(unique(pd$term)))
p <- ggplot(pd, aes(est, term, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) + facet_grid(model ~ outcome, scales="free_y", space="free_y") +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Pre-COVID within-country: subtype + VE (M1) and subtype + protection (M2)",
       subtitle="Protection = VE x 65+ coverage varies by country x season -> better identified than season-level VE. SD units, 95% CrI.",
       x="effect (SD units)", y=NULL) + theme_minimal(base_size=9)
ggsave("output/precovid_ve_subtype.png", p, width=12, height=6, dpi=110)
cat("\nfigures/tables -> output/precovid_ve_subtype.png ; output/precovid_ve_subtype_model{1,2}.csv\n")
