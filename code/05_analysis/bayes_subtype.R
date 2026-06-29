# bayes_subtype.R
#
# Bayesian hierarchical model: does the dominant influenza subtype predict each descriptor, WITHIN
# country? Model:  y_i = X_i beta + u_country[i] + eps_i ,  u ~ N(0, tau^2), eps ~ N(0, sigma^2);
# X = intercept + dominant-subtype dummies (population-level / shared fixed effects), country random
# intercept = partial pooling. Weak priors (beta ~ N(0,100), sigma^2/tau^2 ~ InvGamma(0.01,0.01)).
# Fit by a Gibbs sampler (3 chains); reported as posterior subtype CONTRASTS with 95% credible
# intervals, cross-validated against lme4 fixed effects.
#
# HEAVY CAVEAT: subtype is ~season-determined here (2023/24=H1N1, 2024/25=B, 2025/26=H3N2 across
# almost all countries), so a "subtype effect" is almost inseparable from a SEASON effect (only 3
# post-COVID seasons). Read these as among-season differences labelled by subtype, NOT causal subtype
# effects. Outcomes standardised (z; log for AUC/peak) -> contrasts in SD units. Run from repo root.
suppressMessages({library(dplyr); library(lme4)})
set.seed(1)

d <- read.csv("output/descriptors_subtype.csv", stringsAsFactors=FALSE) %>%
  filter(!is.na(dominant)) %>%
  mutate(auc=log(auc), peak_height=log(peak_height),
         dominant=factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")))
cat(sprintf("n = %d country-seasons, %d countries, subtypes: %s\n",
            nrow(d), n_distinct(d$country), paste(table(d$dominant), names(table(d$dominant)), collapse="  ")))

gibbs <- function(y, X, g, n_iter=6000, n_burn=2000, chains=3){
  n<-length(y); p<-ncol(X); G<-max(g); XtX<-crossprod(X)
  keep<-list()
  for (ch in 1:chains){
    beta<-rep(0,p); u<-rep(0,G); s2<-var(y); t2<-var(y)/2
    M<-matrix(NA,n_iter-n_burn,p)
    for (it in 1:n_iter){
      V<-chol2inv(chol(XtX/s2 + diag(1/100,p))); m<-V%*%(crossprod(X,y-u[g])/s2)
      beta<-as.numeric(m + t(chol(V))%*%rnorm(p))
      e<-as.numeric(y - X%*%beta)
      for (c in 1:G){ idx<-which(g==c); vc<-1/(length(idx)/s2 + 1/t2); u[c]<-rnorm(1, vc*sum(e[idx])/s2, sqrt(vc)) }
      r<-as.numeric(y - X%*%beta - u[g]); s2<-1/rgamma(1,0.01+n/2,0.01+sum(r^2)/2)
      t2<-1/rgamma(1,0.01+G/2,0.01+sum(u^2)/2)
      if (it>n_burn) M[it-n_burn,]<-beta
    }
    keep[[ch]]<-M
  }
  keep
}
rhat <- function(chs){ L<-nrow(chs[[1]]); m<-length(chs); cm<-sapply(chs,colMeans)
  B<-L*apply(cm,1,var); W<-rowMeans(sapply(chs,function(M) apply(M,2,var)))
  sqrt(((L-1)/L*W + B/L)/W) }

g <- as.integer(factor(d$country))
outs <- c(auc="AUC (log)", peak_height="peak height (log)", peak_week="peak week",
          onset_week="onset week", steepness="steepness")
res <- list()
for (o in names(outs)){
  y <- as.numeric(scale(d[[o]])); X <- model.matrix(~ dominant, d)   # ref = A(H1N1)
  ch <- gibbs(y, X, g); A <- do.call(rbind, ch); rh <- rhat(ch)
  # contrasts: cols 2=H3N2-H1N1, 3=B-H1N1 ; derive B-H3N2
  draws <- cbind(`H3N2 - H1N1`=A[,2], `B - H1N1`=A[,3], `B - H3N2`=A[,3]-A[,2])
  for (cn in colnames(draws)){
    q <- quantile(draws[,cn], c(.5,.025,.975))
    res[[length(res)+1]] <- data.frame(outcome=outs[o], contrast=cn, est=round(q[1],2),
      lo=round(q[2],2), hi=round(q[3],2), excl0=ifelse(q[2]>0|q[3]<0,"*",""),
      rhat=round(max(rh),3))
  }
  # lme4 cross-check
  m <- suppressWarnings(lmer(y ~ dominant + (1|country), d, REML=TRUE))
  cat(sprintf("  [%s] lme4 fixef (H3N2-H1N1, B-H1N1): %.2f, %.2f | gibbs: %.2f, %.2f\n", o,
      fixef(m)[2], fixef(m)[3], mean(A[,2]), mean(A[,3])))
}
out <- do.call(rbind, res); rownames(out)<-NULL
cat("\n=== Bayesian within-country subtype contrasts (SD units, 95% CrI) ===\n")
print(out, row.names=FALSE)
write.csv(out, "output/bayes_subtype_contrasts.csv", row.names=FALSE)

# forest plot of the subtype contrasts
suppressMessages(library(ggplot2))
ff <- read.csv("output/bayes_subtype_contrasts.csv", stringsAsFactors=FALSE)
ff$contrast <- factor(ff$contrast, levels=rev(c("H3N2 - H1N1","B - H1N1","B - H3N2")))
p <- ggplot(ff, aes(est, contrast, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) + facet_wrap(~outcome, nrow=1) +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Bayesian within-country subtype contrasts (SD units, 95% credible intervals)",
       subtitle="CAVEAT: subtype is ~season-determined (3 post-COVID seasons) -> these are among-SEASON differences labelled by subtype",
       x="contrast (SD units)", y=NULL) + theme_minimal(base_size=10)
ggsave("output/bayes_subtype.png", p, width=13, height=3.6, dpi=110)
cat("figure -> output/bayes_subtype.png\n")
