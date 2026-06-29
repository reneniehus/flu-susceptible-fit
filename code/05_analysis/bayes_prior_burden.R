# bayes_prior_burden.R
#
# Does the PREVIOUS season's burden predict the current season's descriptors, within country?
# (Susceptible-depletion / immunity-carryover: a big season should leave fewer susceptibles -> a
# smaller / later / less-steep next season.) Predictor = log AUC of the immediately preceding season
# in the same country; lagged by exactly one year, so the COVID gap auto-excludes invalid lags
# (2023/24 <- 2022/23 is dropped). Everything WITHIN country.
#
# Bayesian hierarchical model borrowing strength across all countries: y = a_c + b_c * x_w + eps,
# with country random intercept a_c ~ N(mu_a, tau_a^2) AND random slope b_c ~ N(mu_b, tau_b^2); the
# shared slope mu_b (the EU/EEA-wide within-country effect) learns from every country, country slopes
# shrink toward it. Fit by Gibbs (3 chains), cross-validated against lme4. x and y standardised (z;
# log for AUC/peak), x country-demeaned -> mu_b is the within-country slope in SD units.
#
# CAVEAT: burden magnitude is also reporting-persistent year to year (a country reporting high one
# season tends to again) -> a POSITIVE prior->current AUC/peak slope can be reporting, a NEGATIVE one
# is the immunity signal; timing/steepness outcomes are cleaner. Run from the repo root.
suppressMessages({library(dplyr); library(lme4); library(ggplot2)}); set.seed(1)

d <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE) %>% mutate(syr=as.integer(substr(season,1,4)))
prior <- d %>% transmute(country, syr_next=syr+1, prior_lauc=log(auc))
d <- d %>% left_join(prior, by=c("country"="country","syr"="syr_next")) %>% filter(is.finite(prior_lauc))
cat(sprintf("valid lagged country-seasons: %d (%d countries); seasons: %s\n",
            nrow(d), n_distinct(d$country), paste(sort(unique(d$season)), collapse=" ")))

gibbs_rs <- function(y, x, g, n_iter=8000, n_burn=3000, chains=3){
  n<-length(y); G<-max(g); draws<-vector("list",chains)
  for (ch in 1:chains){
    a<-rep(0,G); b<-rep(0,G); mua<-0; mub<-0; ta2<-var(y); tb2<-0.5; s2<-var(y)
    M<-matrix(NA,n_iter-n_burn,2)
    for (it in 1:n_iter){
      for (c in 1:G){ idx<-which(g==c); Z<-cbind(1,x[idx])
        Pp<-diag(c(1/ta2,1/tb2)); V<-chol2inv(chol(crossprod(Z)/s2+Pp))
        m<-V%*%(crossprod(Z,y[idx])/s2+Pp%*%c(mua,mub))
        ab<-as.numeric(m+t(chol(V))%*%rnorm(2)); a[c]<-ab[1]; b[c]<-ab[2] }
      va<-1/(G/ta2+1/100); mua<-rnorm(1,va*sum(a)/ta2,sqrt(va))
      vb<-1/(G/tb2+1/100); mub<-rnorm(1,vb*sum(b)/tb2,sqrt(vb))
      ta2<-1/rgamma(1,0.01+G/2,0.01+sum((a-mua)^2)/2)
      tb2<-1/rgamma(1,0.01+G/2,0.01+sum((b-mub)^2)/2)
      r<-y-(a[g]+b[g]*x); s2<-1/rgamma(1,0.01+n/2,0.01+sum(r^2)/2)
      if (it>n_burn) M[it-n_burn,]<-c(mub,sqrt(tb2))
    }
    draws[[ch]]<-M
  }
  draws
}
rhat1 <- function(chs,col=1){ L<-nrow(chs[[1]]); cm<-sapply(chs,function(M)mean(M[,col]))
  B<-L*var(cm); W<-mean(sapply(chs,function(M)var(M[,col]))); sqrt(((L-1)/L*W+B/L)/W) }

g <- as.integer(factor(d$country)); xz <- as.numeric(scale(d$prior_lauc))
d$xw <- xz - ave(xz, d$country)                                   # within-country prior burden (SD units)
outs <- c(auc="current AUC (log)", peak_height="current peak height (log)", peak_week="current peak week",
          onset_week="current onset week", steepness="current steepness")
res <- list()
for (o in names(outs)){
  y <- as.numeric(scale(if (o %in% c("auc","peak_height")) log(d[[o]]) else d[[o]]))
  ch <- gibbs_rs(y, d$xw, g); A <- do.call(rbind, ch)
  q <- quantile(A[,1], c(.5,.025,.975))
  res[[length(res)+1]] <- data.frame(outcome=outs[o], beta=round(q[1],2), lo=round(q[2],2), hi=round(q[3],2),
    excl0=ifelse(q[2]>0|q[3]<0,"*",""), het_sd=round(median(A[,2]),2), rhat=round(rhat1(ch),3))
  m <- suppressWarnings(lmer(y ~ xw + (xw|country), d, REML=TRUE))
  cat(sprintf("  [%-12s] lme4 slope %.2f | gibbs mu_b %.2f\n", o, fixef(m)["xw"], mean(A[,1])))
}
out <- do.call(rbind, res); rownames(out)<-NULL
cat("\n=== prior-season burden -> current descriptor: within-country shared slope (SD units, 95% CrI) ===\n")
print(out, row.names=FALSE); write.csv(out, "output/bayes_prior_burden.csv", row.names=FALSE)

p <- ggplot(out, aes(beta, outcome, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Does last season's burden predict this season? (within-country, partial pooling)",
       subtitle=sprintf("shared EU/EEA slope from %d lagged country-seasons; negative AUC/peak = immunity, positive = reporting persistence", nrow(d)),
       x="within-country shared slope per +1 SD prior log-burden (95% CrI)", y=NULL) + theme_minimal(base_size=11)
ggsave("output/bayes_prior_burden.png", p, width=11, height=4, dpi=110)
cat("figure -> output/bayes_prior_burden.png\n")
