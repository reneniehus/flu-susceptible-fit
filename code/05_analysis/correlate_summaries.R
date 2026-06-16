# correlate_summaries.R  --  SCAFFOLD (not implemented yet)
#
# Downstream analysis layer. It consumes the cross-method per-season summary table produced by
# run_all_methods() (code/01_main_supporting/methods_registry.R):
#
#   summaries: one row per (method, country, season) with columns
#              method, country, season, S0, R_eff, c, peak_week, onset_week, cor, process_noise, ...
#
# Planned flow (to be built later, once enough countries x seasons are fitted):
#   1. PLOT the summary statistics across countries and seasons (e.g. S0 by season, faceted by
#      country; peak/onset week heatmaps) to eyeball structure.
#   2. JOIN external factors per (country, season): vaccination coverage (data/vax_flu_*.csv via the
#      data layer) and dominant influenza subtype (erviss_flu_type_subtype / variants).
#   3. CORRELATE the fitted summaries (esp. S0 / R_eff) with those factors -- visually (scatter +
#      smooth) and as statistics (cor.test / rank correlation), keeping it simple and transparent.
#   4. REPORT detected correlations (figures + a short statistical summary) and interpret them
#      (e.g. "lower fitted susceptibility in seasons with higher vaccination coverage").
#
# Deliberately left as a scaffold: the method/summary contract above is the stable interface it
# will build on, so methods can keep changing without touching this layer.

# ---- |-entry point (stub) ----
analyse_summaries = function(summaries, params = NULL, data = NULL){
  stop("correlate_summaries.R: not implemented yet -- see the planned flow in the header.")
}
