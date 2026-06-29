# prepare_descriptors.R
#
# Build the analysis table: the descriptive-method features for every country-season in the combined
# panel (data/slim_flu_iliplus.csv), merged with observed 65+ influenza vaccination coverage.
# Writes output/descriptors.csv (features only) and output/descriptors_vax.csv (features + coverage).
#
# NOTE on coverage availability: observed 65+ coverage runs 2012/13-2021/22, so it overlaps the
# analysis seasons only for the PRE-COVID block (2014/15-2018/19); the post-COVID seasons
# (2023/24-2025/26) carry NA coverage. Coverage is 65+ only, while ILI+ is all-age.
#
# Run from the repo root:  Rscript code/05_analysis/prepare_descriptors.R

source("code/01_main_supporting/setup.R")
source("code/02_settings/settings_version0.R"); params <- settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_sir_deterministic.R")
source("code/01_main_supporting/methods/method_sir_ekf.R")
source("code/01_main_supporting/methods/method_descriptive.R")
source("code/01_main_supporting/methods_registry.R")
models_in <- readRDS("output/models_in.rds")

# ---- descriptive features per country-season ----
slim <- read.csv("data/slim_flu_iliplus.csv", stringsAsFactors=FALSE)
src  <- distinct(slim, country_short, season, source)
rows <- list()
for (cc in sort(unique(slim$country_short))){
  sl  <- load_flu_iliplus_slim(cc)
  s   <- summarise_method_fit(run_method("descriptive", sl, params))
  s$n_weeks <- sapply(sl$ylist, function(y) sum(is.finite(y)))
  rows[[cc]] <- s[, c("country","season","auc","peak_height","peak_week","onset_week","steepness","cor","n_weeks")]
}
d <- do.call(rbind, rows) %>%
  left_join(src, by=c("country"="country_short","season"))
write.csv(d, "output/descriptors.csv", row.names=FALSE)

# ---- merge observed 65+ vaccination coverage ----
vax <- models_in$data_timeseries_long %>%
  filter(indicator=="vaccine_coverage", stream=="vaccination_history_65plus", scenario=="observed_history") %>%
  transmute(country=country_short, season, vax_cov_65=as.numeric(value))
dv <- d %>% left_join(vax, by=c("country","season"))
write.csv(dv, "output/descriptors_vax.csv", row.names=FALSE)

cat(sprintf("wrote output/descriptors.csv (%d rows) and output/descriptors_vax.csv\n", nrow(d)))
cat(sprintf("country-seasons with vaccination coverage: %d (seasons: %s)\n",
            sum(is.finite(dv$vax_cov_65)),
            paste(sort(unique(dv$season[is.finite(dv$vax_cov_65)])), collapse=" ")))
cat(sprintf("vaccination coverage (65+) range: %.2f - %.2f\n", min(dv$vax_cov_65,na.rm=TRUE), max(dv$vax_cov_65,na.rm=TRUE)))
