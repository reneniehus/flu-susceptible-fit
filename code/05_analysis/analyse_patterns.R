# analyse_patterns.R -- core statistical patterns among the descriptive features (and vs vaccination).
# AUC/peak_height are reporting-scale dependent -> analysed WITHIN country (log + country-demeaned);
# steepness/onset_week/peak_week are scale-free -> also compared ACROSS countries. Run from repo root.
suppressMessages({library(dplyr)})
descriptors <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE) %>%
  mutate(log_auc=log(auc), log_peak=log(peak_height))   # covid_era (by season) + source come from prepare_descriptors.R
feature_vars <- c("log_auc","log_peak","onset_week","peak_week","steepness")

cat("=== n =", nrow(descriptors), "country-seasons,", n_distinct(descriptors$country),"countries ===\n")

# within-country: demean each variable by country, then correlate (season-to-season co-variation)
within_country <- descriptors %>% group_by(country) %>% mutate(across(all_of(feature_vars), ~ .x - mean(.x, na.rm=TRUE))) %>% ungroup()
cat("\n--- WITHIN-country correlations (country-demeaned) ---\n")
print(round(cor(within_country[,feature_vars], use="pairwise"),2))

# across-country: country means of the scale-free features
country_means <- descriptors %>% group_by(country) %>% summarise(across(c(onset_week,peak_week,steepness), ~ mean(.x, na.rm=TRUE)), .groups="drop")
cat("\n--- ACROSS-country correlations (country means, scale-free only) ---\n")
print(round(cor(country_means[,c("onset_week","peak_week","steepness")]),2))

# headline questions
within_cor <- function(x,y) round(cor(within_country[[x]], within_country[[y]], use="pairwise"),2)
cat(sprintf("\nQ steep incline ~ high burden? within-country cor(steepness, log AUC) = %.2f ; cor(steepness, log peak) = %.2f\n", within_cor("steepness","log_auc"), within_cor("steepness","log_peak")))
cat(sprintf("Q early season meaning? within-country cor(onset_week, log AUC)=%.2f  cor(onset_week, steepness)=%.2f  cor(onset_week, peak_week)=%.2f  cor(onset_week, log peak)=%.2f\n",
            within_cor("onset_week","log_auc"), within_cor("onset_week","steepness"), within_cor("onset_week","peak_week"), within_cor("onset_week","log_peak")))

# era vs source: the two groupings diverge exactly on 2023/24 (post-COVID but RespiCompass-sourced
# for 20/22 countries), so print BOTH -- the season-era contrast is the substantive question, the
# source contrast is the measurement-shift check (they used to be conflated in one 'era' variable,
# which filed most of 2023/24 under 'pre' and understated the post-COVID onset shift).
cat("\n--- COVID-ERA check (by season; post = 2023/24 onward): scale-free descriptor means ---\n")
print(descriptors %>% group_by(covid_era) %>% summarise(n=n(), onset=round(mean(onset_week),1), peak_wk=round(mean(peak_week),1),
                                        steep=round(mean(steepness),2), .groups="drop") %>% as.data.frame())
cat("\n--- SOURCE check (measurement stream): the same means by RespiCompass vs ERVISS ---\n")
print(descriptors %>% group_by(source) %>% summarise(n=n(), onset=round(mean(onset_week),1), peak_wk=round(mean(peak_week),1),
                                        steep=round(mean(steepness),2), .groups="drop") %>% as.data.frame())

# vaccination (pre-COVID block only; live count printed below): cross-country (country means) and within-country
vax_seasons <- descriptors %>% filter(is.finite(vax_cov_65))
cat(sprintf("\n=== Vaccination linkage: %d country-seasons (%d countries) ===\n", nrow(vax_seasons), n_distinct(vax_seasons$country)))
vax_country_means <- vax_seasons %>% group_by(country) %>% summarise(vax=mean(vax_cov_65), onset=mean(onset_week), steep=mean(steepness),
                                             peak_wk=mean(peak_week), .groups="drop")
cat("ACROSS-country (country means): cor(vax, onset)=", round(cor(vax_country_means$vax,vax_country_means$onset),2),
    " cor(vax, steepness)=", round(cor(vax_country_means$vax,vax_country_means$steep),2), " cor(vax, peak_week)=", round(cor(vax_country_means$vax,vax_country_means$peak_wk),2), "\n")
vax_within <- vax_seasons %>% group_by(country) %>% mutate(across(c(vax_cov_65,log_auc,log_peak,onset_week,peak_week,steepness), ~ .x-mean(.x,na.rm=TRUE))) %>% ungroup()
cat("WITHIN-country (demeaned): cor(vax, log AUC)=", round(cor(vax_within$vax_cov_65,vax_within$log_auc,use="pairwise"),2),
    " cor(vax, steepness)=", round(cor(vax_within$vax_cov_65,vax_within$steepness,use="pairwise"),2),
    " cor(vax, onset)=", round(cor(vax_within$vax_cov_65,vax_within$onset_week,use="pairwise"),2), "\n")
