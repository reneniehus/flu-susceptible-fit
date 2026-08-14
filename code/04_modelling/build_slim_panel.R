# build_slim_panel.R
#
# Build the committed combined ILI+ panel (data/slim_flu_iliplus.csv) from models_in. Every data
# decision (COVID exclusion, >=15-positive-week rule, source stitch + overrides, unit scaling)
# lives in ONE place: stitch_iliplus_panel() in code/01_main_supporting/stitch_iliplus.R --
# this script just applies it with the analysis defaults and writes the CSV.
# tests/testthat/test-slim-panel.R requires the rebuild to reproduce the committed file exactly.
#
# Run from the repo root:  Rscript code/04_modelling/build_slim_panel.R

source("code/01_main_supporting/setup.R")
source("code/01_main_supporting/stitch_iliplus.R")
models_in <- readRDS("output/models_in.rds")

slim <- stitch_iliplus_panel(models_in, exclude_covid = TRUE, min_wk = 15)

write.csv(slim, "data/slim_flu_iliplus.csv", row.names=FALSE)
cat(sprintf("wrote data/slim_flu_iliplus.csv: %d rows | %d countries | %d seasons | %d country-seasons\n",
            nrow(slim), n_distinct(slim$country_short), n_distinct(slim$season),
            nrow(distinct(slim, country_short, season))))
slim %>% distinct(country_short, season, source) %>% count(season, source) %>%
  tidyr::pivot_wider(names_from=source, values_from=n, values_fill=0) %>% as.data.frame() %>% print(row.names=FALSE)
