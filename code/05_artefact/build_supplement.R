# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Build the RespiCast-value supplement: embed the figures, write a print HTML ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Injects the five PNGs (base64) into the supplement template so the page is fully
# self-contained, then a headless-Chrome step (see README / the shell call) prints it
# to output/respicast_value_supplement.pdf.
# Run:  Rscript code/05_figures/fig_respicast_value.R && Rscript code/05_artefact/build_supplement.R

source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()

TEMPLATE <- "code/05_artefact/respicast_supplement.html"
OUT_HTML <- file.path(params$output_dir, "respicast_value_supplement.html")
fig_dir  <- file.path(params$figure_dir, "respicast")
figs     <- c(FIG1 = "fig1_capacity.png", FIG2 = "fig2_preference.png", FIG3 = "fig3_engagement.png",
              FIG4 = "fig4_decisions.png", FIG5 = "fig5_integration.png")

data_uri <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  paste0("data:image/png;base64,", jsonlite::base64_enc(raw))
}

html <- paste(readLines(TEMPLATE, warn = FALSE), collapse = "\n")
for (tok in names(figs)) {
  html <- gsub(paste0("__", tok, "__"), data_uri(file.path(fig_dir, figs[[tok]])), html, fixed = TRUE)
}
writeLines(html, OUT_HTML)
cat(sprintf("wrote %s (%.0f KB)\n", OUT_HTML, file.info(OUT_HTML)$size / 1024))
