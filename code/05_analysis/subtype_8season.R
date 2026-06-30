# subtype_8season.R
#
# Extends the dominant-subtype -> descriptor analysis to ALL EIGHT panel seasons by attaching a
# CONTINENTAL dominant subtype per season (data/external/dominant_subtype_by_season.csv: pre-COVID
# 2014/15-2018/19 from ECDC/Eurosurveillance literature, post-COVID 2023/24-2025/26 from ERVISS typing).
# The user-sanctioned assumption is that the dominant subtype is ~identical across EU/EEA countries
# within a season; we broadcast the continental label to every country in that season.
#
# WHY THIS HELPS: the original subtype analysis (bayes_subtype.R) used only 3 post-COVID seasons, each
# one subtype -> subtype was perfectly confounded with season. Across 8 seasons each subtype now RECURS
# in BOTH eras (A(H1N1): 2015/16, 2018/19, 2023/24 | A(H3N2): 2014/15, 2016/17, 2025/26 | B: 2017/18,
# 2024/25), so subtype is no longer collinear with either the calendar or the RespiCompass->ERVISS
# source change. We add an `era` covariate to absorb the source shift and read subtype contrasts net of it.
#
# Model (within country, partial pooling): y = X beta + u_country + eps, X = intercept + subtype dummies
# (+ era), u ~ N(0, tau^2). Gibbs sampler (3 chains), cross-checked against lme4. Outcomes standardised
# (z; log for AUC/peak) -> contrasts in SD units. CAVEAT still applies: only 8 seasons, continental
# subtype, so these remain among-season differences labelled by subtype, now much less season-confounded.
# Run from the repo root:  Rscript code/05_analysis/subtype_8season.R
suppressMessages({library(dplyr); library(lme4)}); set.seed(1)

desc <- read.csv("output/descriptors.csv", stringsAsFactors=FALSE)
sub  <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE) %>%
  transmute(season, dominant)
d <- desc %>% inner_join(sub, by="season") %>%
  mutate(era      = ifelse(source=="RespiCompass","pre","post"),
         auc      = log(auc), peak_height = log(peak_height),
         dominant = factor(dominant, levels=c("A(H1N1)","A(H3N2)","B")),
         era      = factor(era, levels=c("pre","post")))
cat(sprintf("n = %d country-seasons | %d countries | %d seasons\n", nrow(d), n_distinct(d$country), n_distinct(d$season)))
cat("subtype x era table (country-seasons):\n"); print(table(d$dominant, d$era))

gibbs <- function(y, X, g, n_iter=6000, n_burn=2000, chains=3){
  n<-length(y); p<-ncol(X); G<-max(g); XtX<-crossprod(X); keep<-list()
  for (ch in 1:chains){
    beta<-rep(0,p); u<-rep(0,G); s2<-var(y); t2<-var(y)/2
    M<-matrix(NA,n_iter-n_burn,p)
    for (it in 1:n_iter){
      V<-chol2inv(chol(XtX/s2 + diag(1/100,p))); m<-V%*%(crossprod(X,y-u[g])/s2)
      beta<-as.numeric(m + t(chol(V))%*%rnorm(p))
      e<-as.numeric(y - X%*%beta)
      for (gi in 1:G){ idx<-which(g==gi); vc<-1/(length(idx)/s2 + 1/t2); u[gi]<-rnorm(1, vc*sum(e[idx])/s2, sqrt(vc)) }
      r<-as.numeric(y - X%*%beta - u[g]); s2<-1/rgamma(1,0.01+n/2,0.01+sum(r^2)/2)
      t2<-1/rgamma(1,0.01+G/2,0.01+sum(u^2)/2)
      if (it>n_burn) M[it-n_burn,]<-beta
    }
    keep[[ch]]<-M
  }
  keep
}
rhat <- function(chs){ L<-nrow(chs[[1]]); cm<-sapply(chs,colMeans)
  B<-L*apply(cm,1,var); W<-rowMeans(sapply(chs,function(M) apply(M,2,var))); sqrt(((L-1)/L*W + B/L)/W) }

g <- as.integer(factor(d$country))
outs <- c(auc="AUC (log)", peak_height="peak height (log)", peak_week="peak week",
          onset_week="onset week", steepness="steepness")
res <- list()
for (o in names(outs)){
  y <- as.numeric(scale(d[[o]])); X <- model.matrix(~ dominant + era, d)   # ref = A(H1N1), era controls source shift
  ch <- gibbs(y, X, g); A <- do.call(rbind, ch); rh <- rhat(ch)
  ci <- grep("^dominant", colnames(X))           # columns 2 = H3N2-H1N1, 3 = B-H1N1
  draws <- cbind(`H3N2 - H1N1`=A[,ci[1]], `B - H1N1`=A[,ci[2]], `B - H3N2`=A[,ci[2]]-A[,ci[1]])
  for (cn in colnames(draws)){
    q <- quantile(draws[,cn], c(.5,.025,.975))
    res[[length(res)+1]] <- data.frame(outcome=outs[o], contrast=cn, est=round(q[1],2),
      lo=round(q[2],2), hi=round(q[3],2), excl0=ifelse(q[2]>0|q[3]<0,"*",""), rhat=round(max(rh),3))
  }
  m <- suppressWarnings(lmer(y ~ dominant + era + (1|country), d, REML=TRUE))
  cat(sprintf("  [%-12s] lme4 (H3N2-H1N1, B-H1N1): %.2f, %.2f | gibbs: %.2f, %.2f\n",
              o, fixef(m)[2], fixef(m)[3], mean(A[,ci[1]]), mean(A[,ci[2]])))
}
out <- do.call(rbind, res); rownames(out)<-NULL
cat("\n=== Within-country subtype contrasts across 8 seasons, net of era (SD units, 95% CrI) ===\n")
print(out, row.names=FALSE)
write.csv(out, "output/subtype_8season_contrasts.csv", row.names=FALSE)

# forest plot
suppressMessages(library(ggplot2))
out$contrast <- factor(out$contrast, levels=rev(c("H3N2 - H1N1","B - H1N1","B - H3N2")))
p <- ggplot(out, aes(est, contrast, color=excl0=="*")) + geom_vline(xintercept=0, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi)) + facet_wrap(~outcome, nrow=1) +
  scale_color_manual(values=c("FALSE"="grey55","TRUE"="#d95f02"), guide="none") +
  labs(title="Within-country dominant-subtype contrasts across 8 seasons (net of era)",
       subtitle="Continental subtype per season; each subtype recurs in both eras -> far less season-confounded than the 3-season analysis",
       x="contrast (SD units, 95% CrI)", y=NULL) + theme_minimal(base_size=10)
ggsave("output/subtype_8season.png", p, width=13, height=3.6, dpi=110)
cat("figure -> output/subtype_8season.png\n")
