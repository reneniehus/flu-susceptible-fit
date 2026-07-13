# Splice the extracted data into the HTML template to produce the self-contained guidebook.
# Run from the project root, after make_guidebook_data.R:
#   Rscript guidebook/make_guidebook_data.R   # runs one default pandemic -> guidebook/guidebook_data.json
#   Rscript guidebook/build.R                 # -> guidebook/guidebook.html (open in a browser)
# We split on the __DATA__ placeholder and paste (rather than sub/gsub) so nothing in the JSON is ever
# interpreted as a regex backreference.
tpl   <- paste(readLines("guidebook/guidebook_template.html", warn = FALSE), collapse = "\n")
data  <- paste(readLines("guidebook/guidebook_data.json", warn = FALSE), collapse = "\n")
parts <- strsplit(tpl, "__DATA__", fixed = TRUE)[[1]]
if (length(parts) != 2) stop("build.R: expected exactly one __DATA__ placeholder in the template")
writeLines(paste0(parts[1], data, parts[2]), "guidebook/guidebook.html")
cat("wrote guidebook/guidebook.html (", file.size("guidebook/guidebook.html"), "bytes )\n")
