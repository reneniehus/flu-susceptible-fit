# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Build the dashboard page: inject the data JSON into the HTML template ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The template (code/05_artefact/dashboard_template.html) carries all the CSS + SVG
# rendering JS and a single placeholder token where the embedded data blob goes.
# We splice output/artefact_data.json in (string paste, not regex -- the JSON has
# backslashes/unicode that sub() would mangle) and write the self-contained page.
# Run:  Rscript code/05_artefact/build_pages.R

TEMPLATE <- "code/05_artefact/dashboard_template.html"
DATA     <- "output/artefact_data.json"
OUT      <- "artefact/dashboard.html"
TOKEN    <- "__VIZ_DATA__"

tpl  <- paste(readLines(TEMPLATE, warn = FALSE), collapse = "\n")
data <- paste(readLines(DATA,     warn = FALSE), collapse = "\n")

stopifnot(grepl(TOKEN, tpl, fixed = TRUE))
if (grepl("</script", data, fixed = TRUE))
  stop("data JSON contains a </script sequence -- would break the embedding.")

parts <- strsplit(tpl, TOKEN, fixed = TRUE)[[1]]
page  <- paste0(parts[1], data, parts[2])

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
writeLines(page, OUT)
cat(sprintf("built %s  (%.0f KB; data %.0f KB)\n",
            OUT, file.info(OUT)$size / 1024, file.info(DATA)$size / 1024))
