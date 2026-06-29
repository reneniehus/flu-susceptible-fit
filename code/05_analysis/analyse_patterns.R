# analyse_patterns.R -- core statistical patterns among the descriptive features (and vs vaccination).
# AUC/peak_height are reporting-scale dependent -> analysed WITHIN country (log + country-demeaned);
# steepness/onset_week/peak_week are scale-free -> also compared ACROSS countries. Run from repo root.
suppressMessages({library(dplyr); library(tidyr)})
d <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>%
  mutate(lauc=log(auc), lpk=log(peak_height), era=ifelse(source=="RespiCompass","pre","post"))
V <- c("lauc","lpk","onset_week","peak_week","steepness")

cat("=== n =", nrow(d), "country-seasons,", n_distinct(d$country),"countries ===\n")

# within-country: demean each variable by country, then correlate (season-to-season co-variation)
dw <- d %>% group_by(country) %>% mutate(across(all_of(V), ~ .x - mean(.x, na.rm=TRUE))) %>% ungroup()
cat("\n--- WITHIN-country correlations (country-demeaned) ---\n")
print(round(cor(dw[,V], use="pairwise"),2))

# across-country: country means of the scale-free features
cm <- d %>% group_by(country) %>% summarise(across(c(onset_week,peak_week,steepness), mean, na.rm=TRUE), .groups="drop")
cat("\n--- ACROSS-country correlations (country means, scale-free only) ---\n")
print(round(cor(cm[,c("onset_week","peak_week","steepness")]),2))

# headline questions
wc <- function(a,b) round(cor(dw[[a]], dw[[b]], use="pairwise"),2)
cat(sprintf("\nQ steep incline ~ high burden? within-country cor(steepness, log AUC) = %.2f ; cor(steepness, log peak) = %.2f\n", wc("steepness","lauc"), wc("steepness","lpk")))
cat(sprintf("Q early season meaning? within-country cor(onset_week, log AUC)=%.2f  cor(onset_week, steepness)=%.2f  cor(onset_week, peak_week)=%.2f  cor(onset_week, log peak)=%.2f\n",
            wc("onset_week","lauc"), wc("onset_week","steepness"), wc("onset_week","peak_week"), wc("onset_week","lpk")))

# era / source confound: do scale-free descriptors shift pre vs post?
cat("\n--- ERA check (source confound): scale-free descriptor means pre(RespiCompass) vs post(ERVISS) ---\n")
print(d %>% group_by(era) %>% summarise(n=n(), onset=round(mean(onset_week),1), peak_wk=round(mean(peak_week),1),
                                        steep=round(mean(steepness),2), .groups="drop") %>% as.data.frame())

# vaccination (pre-COVID block only, 97 obs): cross-country (country means) and within-country
dv <- d %>% filter(is.finite(vax_cov_65))
cat(sprintf("\n=== Vaccination linkage: %d country-seasons (%d countries) ===\n", nrow(dv), n_distinct(dv$country)))
vm <- dv %>% group_by(country) %>% summarise(vax=mean(vax_cov_65), onset=mean(onset_week), steep=mean(steepness),
                                             peak_wk=mean(peak_week), .groups="drop")
cat("ACROSS-country (country means): cor(vax, onset)=", round(cor(vm$vax,vm$onset),2),
    " cor(vax, steepness)=", round(cor(vm$vax,vm$steep),2), " cor(vax, peak_week)=", round(cor(vm$vax,vm$peak_wk),2), "\n")
dvw <- dv %>% group_by(country) %>% mutate(across(c(vax_cov_65,lauc,lpk,onset_week,peak_week,steepness), ~ .x-mean(.x,na.rm=TRUE))) %>% ungroup()
cat("WITHIN-country (demeaned): cor(vax, log AUC)=", round(cor(dvw$vax_cov_65,dvw$lauc,use="pairwise"),2),
    " cor(vax, steepness)=", round(cor(dvw$vax_cov_65,dvw$steepness,use="pairwise"),2),
    " cor(vax, onset)=", round(cor(dvw$vax_cov_65,dvw$onset_week,use="pairwise"),2), "\n")
