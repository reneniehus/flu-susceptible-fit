# validate.R
#
# Lightweight contracts that shield the pipeline from upstream data changes (ECDC renaming
# or dropping a column, removing a file). The philosophy is deliberately cheap: a couple of
# column/stream checks that fail LOUDLY and clearly at load time, rather than a heavy schema
# framework. The same helpers back the tests in tests/testthat/.
#
# Adding a column upstream is harmless (rate streams select named columns; detailed streams
# keep everything). Renaming/removing a required column, or a stream disappearing, is what we
# want to catch -- that is what these checks do.

# ---- |-fail clearly if a data frame is missing required columns ----
check_required_columns = function(df, required, what="data frame"){
  missing = setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required column(s): %s\n  (has: %s)",
                 what, paste(missing, collapse=", "), paste(names(df), collapse=", ")),
         call.=FALSE)
  }
  invisible(df)
}

# ---- |-required RAW columns every ERVISS file must have ----
# these are the columns standardise_erviss() needs for every stream. The richer pathogen
# columns (pathogen / pathogentype / pathogensubtype) are NOT universal -- severity,
# sequencing and variants files do not carry them -- so they are checked where they are
# actually consumed instead (extract_typing(), via check_required_columns()).
erviss_required_cols = function(schema=NULL){   # schema is accepted for call-site symmetry but unused: these 5 base cols are required for every schema
  c("countryname", "yearweek", "age", "indicator", "value")
}

# ---- |-post-load contract over the assembled data list ----
# checks that every expected epi stream is present and standardised to the canonical columns.
validate_data = function(data){
  expected_streams = c("erviss_ili_ari", "erviss_sari_rates",
                       "erviss_typing_sentinel", "erviss_typing_nonsentinel", "erviss_typing_sari",
                       "erviss_flu_type_subtype", "erviss_severity_nonsentinel",
                       "erviss_sequencing", "erviss_variants", "respicompass_iliplus")
  missing_streams = setdiff(expected_streams, names(data$epi))
  if (length(missing_streams) > 0) {
    stop("data$epi is missing expected stream(s): ", paste(missing_streams, collapse=", "), call.=FALSE)
  }
  # the slim "rates" streams must carry the canonical analysis columns
  for (nm in c("erviss_ili_ari", "respicompass_iliplus")) {
    check_required_columns(data$epi[[nm]], c("country_short", "date", "target", "agegroup", "value"),
                           paste0("data$epi$", nm))
  }
  invisible(TRUE)
}

# ---- |-contract over the gen_model_input() output (models_in) ----
validate_models_in = function(models_in){
  check_required_columns(models_in$data_timeseries_long,
                         c("country_short", "season", "date", "source", "stream", "pathogen",
                           "indicator", "agegroup", "value", "observed"),
                         "models_in$data_timeseries_long")
  check_required_columns(models_in$data_season_summary,
                         c("country_short", "season", "stream", "indicator",
                           "completeness", "completeness_active"),
                         "models_in$data_season_summary")
  # the ILI+ family must cover every pathogen we expect to construct it for
  iliplus_pathogens = models_in$data_timeseries_long %>%
    filter(indicator == "ili_plus") %>% pull(pathogen) %>% unique()
  missing_pathogens = setdiff(c("Influenza", "SARS-CoV-2", "RSV"), iliplus_pathogens)
  if (length(missing_pathogens) > 0) {
    warning("ILI+ is missing pathogen(s): ", paste(missing_pathogens, collapse=", "), call.=FALSE)
  }
  invisible(TRUE)
}
