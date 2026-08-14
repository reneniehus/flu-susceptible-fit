# analysis_helpers.R
#
# Shared machinery for the code/05_analysis driver models. Before this file the SAME Gibbs
# sampler, rhat and q95 helpers were copy-pasted across bayes_subtype.R, subtype_8season.R and
# bayes_precovid_ve_subtype.R, and the season-level VE values were hard-coded (twice) instead of
# derived from the provenance CSV. Base R + dplyr only; scripts source this AFTER set.seed(), and
# the samplers are byte-identical to the copies they replace, so seeded draws are unchanged.

# ---- |-varying-intercept Gibbs: y = X beta + u_g + eps, u ~ N(0, t2), 3 chains ----
# Weak priors: beta ~ N(0,100), s2/t2 ~ InvGamma(0.01, 0.01). Returns a list of per-chain
# post-burn-in beta draw matrices (rows = draws, cols = columns of X).
gibbs_ri <- function(y, X, g, n_iter=6000, n_burn=2000, chains=3){
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

# ---- |-crossed random intercepts: y = X beta + u_g1 + w_g2 + eps ----
# Same weak priors, plus w ~ N(0, t2b). Built for season-level predictors (subtype, VE): with a
# country intercept ONLY, the model treats the ~20-25 countries sharing one season-level value as
# independent replicates and the CrIs are anti-conservative -- the season intercept restores the
# honest effective replication (the number of seasons). Returns per-chain beta draws.
gibbs_ri2 <- function(y, X, g1, g2, n_iter=6000, n_burn=2000, chains=3){
  n<-length(y); p<-ncol(X); G1<-max(g1); G2<-max(g2); XtX<-crossprod(X); keep<-list()
  for (ch in 1:chains){
    beta<-rep(0,p); u<-rep(0,G1); w<-rep(0,G2); s2<-var(y); t2<-var(y)/2; t2b<-var(y)/2
    M<-matrix(NA,n_iter-n_burn,p)
    for (it in 1:n_iter){
      V<-chol2inv(chol(XtX/s2 + diag(1/100,p))); m<-V%*%(crossprod(X, y-u[g1]-w[g2])/s2)
      beta<-as.numeric(m + t(chol(V))%*%rnorm(p))
      e<-as.numeric(y - X%*%beta - w[g2])
      for (gi in 1:G1){ idx<-which(g1==gi); vc<-1/(length(idx)/s2 + 1/t2); u[gi]<-rnorm(1, vc*sum(e[idx])/s2, sqrt(vc)) }
      e<-as.numeric(y - X%*%beta - u[g1])
      for (gi in 1:G2){ idx<-which(g2==gi); vc<-1/(length(idx)/s2 + 1/t2b); w[gi]<-rnorm(1, vc*sum(e[idx])/s2, sqrt(vc)) }
      r<-as.numeric(y - X%*%beta - u[g1] - w[g2]); s2<-1/rgamma(1,0.01+n/2,0.01+sum(r^2)/2)
      t2 <-1/rgamma(1,0.01+G1/2,0.01+sum(u^2)/2)
      t2b<-1/rgamma(1,0.01+G2/2,0.01+sum(w^2)/2)
      if (it>n_burn) M[it-n_burn,]<-beta
    }
    keep[[ch]]<-M
  }
  keep
}

# ---- |-split-free Gelman-Rubin rhat over all columns of per-chain draw matrices ----
rhat_cols <- function(chs){ L<-nrow(chs[[1]]); cm<-sapply(chs,colMeans)
  B<-L*apply(cm,1,var); W<-rowMeans(sapply(chs,function(M) apply(M,2,var)))
  sqrt(((L-1)/L*W + B/L)/W) }

# ---- |-rhat of ONE derived quantity: fn maps a chain's draw matrix -> a draw vector ----
# Use for contrasts (e.g. B - H3N2 = col3 - col2): the rhat must be computed on the CONTRAST
# draws themselves, not borrowed from one of the ingredient columns.
rhat_derived <- function(chs, fn){ ds<-lapply(chs, fn); L<-length(ds[[1]])
  cm<-sapply(ds, mean); B<-L*var(cm); W<-mean(sapply(ds, var)); sqrt(((L-1)/L*W + B/L)/W) }

# ---- |-posterior median + 95% interval ----
q95 <- function(x) quantile(x, c(.5,.025,.975))

# ---- |-season-level VE against the DOMINANT subtype, derived from the provenance CSV ----
# One explicit, reproducible rule (no hard-coded vectors): among primary-care rows for the
# season's dominant subtype with age_group 'all' (or 'target_group' as fallback), prefer
#   end-of-season point > interim point > end-of-season study-range midpoint > interim midpoint;
# if the dominant subtype has no row at all, fall back to the same rule over 'influenza_A' /
# 'any' rows (subtype-nonspecific estimates). Returns a named vector season -> VE.
ve_vs_dominant <- function(seasons_dominant,
                           ve_csv = "data/external/vaccine_effectiveness.csv"){
  ve <- read.csv(ve_csv, stringsAsFactors=FALSE)
  subtype_rows <- function(rows, dom){
    key <- c("A(H1N1)"="A(H1N1)pdm09", "A(H3N2)"="A(H3N2)", "B"="B")[dom]
    exact <- rows[rows$subtype == key, ]
    if (dom == "B" && !nrow(exact)) exact <- rows[rows$subtype %in% c("B/Victoria","B/Yamagata"), ]
    if (!nrow(exact)) exact <- rows[rows$subtype %in% c("influenza_A","any"), ]   # subtype-nonspecific fallback
    exact
  }
  pick <- function(season, dom){
    rows <- ve[ve$season == season & ve$setting == "primary_care" & ve$age_group %in% c("all","target_group"), ]
    rows <- subtype_rows(rows, dom)
    if (!nrow(rows)) return(NA_real_)
    rows$val  <- ifelse(rows$value_type == "point", rows$ve_point, (rows$ci_low + rows$ci_high)/2)  # study_range -> midpoint
    rows$pref <- with(rows, 1L*(value_type=="point" & timing=="end_of_season") * 8 +
                            1L*(value_type=="point" & timing!="end_of_season") * 4 +
                            1L*(value_type!="point" & timing=="end_of_season") * 2 +
                            1L*(value_type!="point" & timing!="end_of_season") * 1 +
                            1L*(age_group=="all"))                                # tie-break: 'all' over 'target_group'
    rows$val[which.max(rows$pref)]
  }
  out <- vapply(seq_len(nrow(seasons_dominant)),
                function(i) pick(seasons_dominant$season[i], seasons_dominant$dominant[i]), numeric(1))
  stats::setNames(out, seasons_dominant$season)
}
