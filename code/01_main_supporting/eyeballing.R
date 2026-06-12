# eyeballing.R
#
# Visual "eyeballing" of the canonical tables produced by gen_model_input(). It returns a
# manifest of figures, each bundled with a short title, subtitle and a few narrative bullets,
# so a report (code/03_report/eyeballing_report.Rmd) can render every figure large with its
# findings underneath. Two families of figures:
#   - data-quality figures  : how complete / reliable each indicator is, across countries & seasons
#   - temporal-dynamics figs : how the indicators move over time, by disease and country
#
# The bullets are STATIC narrative placeholders on purpose -- edit them to match what you see.
# Default scope is every country and every season available in the data.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Shared look & feel ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# consistent colour per pathogen across every figure
pathogen_colours = c("Influenza"="#1b9e77", "SARS-CoV-2"="#7570b3", "RSV"="#d95f02")

# shared minimal theme (note: we call ggplot2::ggplot directly to bypass the global
# ggplot() override from setup.R, so colour/fill scales are fully under our control here)
theme_eyeball = function(base_size=14){
  ggplot2::theme_minimal(base_size=base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text      = ggplot2::element_text(face="bold"),
      plot.title      = ggplot2::element_text(face="bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Individual figure builders ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-Quality: completeness heatmaps (country x season x indicator) ----
# shared selection of series + one tile plot, reused for both completeness measures
completeness_panel_df = function(season_summary){
  # curated set of series to show, with a readable panel label and a fixed ordering
  quality_series = tribble(
    ~stream,             ~indicator,            ~pathogen,     ~panel,
    "ili_ari",           "ILIconsultationrate", NA,            "ILI consultation rate",
    "ili_ari",           "ARIconsultationrate", NA,            "ARI consultation rate",
    "typing_sentinel",   "positivity",          "Influenza",   "Influenza positivity (sentinel)",
    "typing_sentinel",   "positivity",          "SARS-CoV-2",  "SARS-CoV-2 positivity (sentinel)",
    "typing_sentinel",   "positivity",          "RSV",         "RSV positivity (sentinel)",
    "ili_plus_sentinel", "ili_plus",            "Influenza",   "ILI+ Influenza (sentinel)"
  )
  key_of = function(stream, indicator, pathogen) paste(stream, indicator, ifelse(is.na(pathogen), "", pathogen), sep="|")
  selected = quality_series %>% mutate(.k=key_of(stream, indicator, pathogen))
  season_summary %>%
    filter(summary_level=="all_agegroups", temporal_resolution=="weekly") %>%
    mutate(.k=key_of(stream, indicator, pathogen)) %>%
    inner_join(selected %>% select(.k, panel), by=".k") %>%
    mutate(panel=factor(panel, levels=quality_series$panel))
}

plot_completeness = function(df, fill_var, title, fill_lab){
  ggplot2::ggplot(df, ggplot2::aes(season, country_short, fill=.data[[fill_var]])) +
    ggplot2::geom_tile(colour="white", linewidth=0.2) +
    ggplot2::facet_wrap(~panel, ncol=3) +
    ggplot2::scale_fill_viridis_c(limits=c(0,1), labels=scales::percent, na.value="grey92") +
    ggplot2::labs(title=title, x=NULL, y=NULL, fill=fill_lab) +
    theme_eyeball() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
}

fig_completeness = function(season_summary){
  plot_completeness(completeness_panel_df(season_summary), "completeness",
                    "Season-window completeness across countries and seasons",
                    "Reported weeks / weeks in season")
}

# alternative measure: completeness within the active reporting span (off-season ignored)
fig_completeness_active = function(season_summary){
  plot_completeness(completeness_panel_df(season_summary), "completeness_active",
                    "Active-span completeness across countries and seasons",
                    "Reported weeks / weeks in active span")
}

# ---- |-Quality 2: typing test volume (country x season) ----
fig_testing_volume = function(timeseries_long){
  # tests are a shared denominator across pathogens, so pick one pathogen to avoid triple-counting
  df = timeseries_long %>%
    filter(indicator=="tests", pathogen=="Influenza",
           stream %in% c("typing_sentinel", "typing_nonsentinel")) %>%
    summarise(total_tests=sum(value, na.rm=TRUE),
              .by=c(country_short, season, stream)) %>%
    mutate(total_tests=ifelse(total_tests <= 0, NA_real_, total_tests),   # 0 tests -> grey, not log10(0)
           stream=recode(stream, typing_sentinel="Sentinel", typing_nonsentinel="Non-sentinel"))

  ggplot2::ggplot(df, ggplot2::aes(season, country_short, fill=total_tests)) +
    ggplot2::geom_tile(colour="white", linewidth=0.2) +
    ggplot2::facet_wrap(~stream, ncol=2) +
    ggplot2::scale_fill_viridis_c(trans="log10", labels=scales::label_comma(), na.value="grey92") +
    ggplot2::labs(title="Typing test volume by country and season",
                  x=NULL, y=NULL, fill="Specimens tested (season total, log scale)") +
    theme_eyeball() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
}

# ---- |-helper: continuous-time, one-panel-per-country line figure ----
# used by the temporal-dynamics figures so that all seasons sit on one continuous axis
# (small multiples over country x season would be unreadable for ~30 countries x ~12 seasons).
facet_country_lines = function(df, colour_var, title, y_lab, colour_lab,
                               colours=NULL, linetype_var=NULL){
  aes_args = list(x=quote(date), y=quote(value), colour=as.name(colour_var))
  if (!is.null(linetype_var)) aes_args$linetype = as.name(linetype_var)
  p = ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_line(na.rm=TRUE, linewidth=0.4) +
    ggplot2::facet_wrap(~country_short, ncol=6, scales="free_y") +
    ggplot2::labs(title=title, x=NULL, y=y_lab, colour=colour_lab, linetype=NULL) +
    theme_eyeball(base_size=13)
  if (!is.null(colours)) p = p + ggplot2::scale_colour_manual(values=colours, na.value="grey60")
  p
}

# ---- |-Dynamics 1: ILI+ by pathogen ----
fig_iliplus_dynamics = function(timeseries_long){
  df = timeseries_long %>%
    filter(indicator=="ili_plus", stream=="ili_plus_sentinel", agegroup=="age_total")
  facet_country_lines(df, colour_var="pathogen",
                      title="ILI+ dynamics by pathogen (ILI rate x sentinel positivity)",
                      y_lab="ILI+", colour_lab="Pathogen", colours=pathogen_colours)
}

# ---- |-Dynamics 2: test positivity by pathogen ----
fig_positivity_dynamics = function(timeseries_long){
  df = timeseries_long %>%
    filter(indicator=="positivity", stream=="typing_sentinel", agegroup=="age_total")
  facet_country_lines(df, colour_var="pathogen",
                      title="Sentinel test positivity by pathogen",
                      y_lab="Positivity", colour_lab="Pathogen", colours=pathogen_colours)
}

# ---- |-Dynamics 3: ILI vs ARI consultation rates ----
fig_syndromic_dynamics = function(timeseries_long){
  df = timeseries_long %>%
    filter(stream=="ili_ari", agegroup=="age_total",
           indicator %in% c("ILIconsultationrate", "ARIconsultationrate"))
  facet_country_lines(df, colour_var="indicator_label",
                      title="Syndromic consultation rates: ILI vs ARI",
                      y_lab="Consultation rate", colour_lab="Indicator")
}

# ---- |-Dynamics 4: timing-aligned ILI+ (countries stacked on a shared time axis) ----
# one country per row, all sharing the same x (time) axis, so peaks line up vertically and
# the relative timing of waves across countries (and pathogens) is read off at a glance.
# Each country x pathogen curve is lightly smoothed (centred 3-week rolling mean) and then
# rescaled to its own peak (=1) so only TIMING is visible, not magnitude (a small RSV wave
# and a large flu wave both reach the top of the panel).
roll_mean_centred = function(x, k=3){
  # centred rolling mean over a window of k points; averages whatever is non-missing in the
  # window (so short gaps and the series edges are handled without dropping points)
  n = length(x); half = (k - 1) %/% 2; out = rep(NA_real_, n)
  for (i in seq_len(n)){
    w = x[max(1, i - half):min(n, i + half)]
    if (any(!is.na(w))) out[i] = mean(w, na.rm=TRUE)
  }
  out
}
fig_iliplus_aligned = function(timeseries_long){
  df = timeseries_long %>%
    filter(indicator=="ili_plus", stream=="ili_plus_sentinel", agegroup=="age_total") %>%
    arrange(country_short, pathogen, date) %>%
    mutate(.by = c(country_short, pathogen),
           value = roll_mean_centred(value, k=3)) %>%             # smooth before scaling
    mutate(.by = c(country_short, pathogen),
           value = { peak = suppressWarnings(max(value, na.rm=TRUE))   # per country x pathogen peak
                     if (!is.finite(peak) || peak <= 0) NA_real_ else value / peak })
  ggplot2::ggplot(df, ggplot2::aes(date, value, colour=pathogen)) +
    ggplot2::geom_line(na.rm=TRUE, linewidth=0.4) +
    ggplot2::facet_grid(country_short ~ ., switch="y") +
    ggplot2::scale_colour_manual(values=pathogen_colours, na.value="grey60") +
    ggplot2::labs(title="ILI+ timing aligned across countries (3-week smoothed, each curve scaled to its own peak)",
                  x=NULL, y="ILI+ (smoothed, scaled to peak = 1, per country & pathogen)", colour="Pathogen") +
    theme_eyeball(base_size=12) +
    ggplot2::theme(strip.text.y.left=ggplot2::element_text(angle=0),
                   axis.text.y=ggplot2::element_blank(), panel.grid.major.y=ggplot2::element_blank())
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Supplement figure builders (extra views worth keeping handy) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-Supplement: age-specific ILI consultation rate ----
fig_age_dynamics = function(timeseries_long){
  df = timeseries_long %>%
    filter(stream=="ili_ari", indicator=="ILIconsultationrate", agegroup!="age_total")
  facet_country_lines(df, colour_var="agegroup",
                      title="Age-specific ILI consultation rate",
                      y_lab="ILI consultation rate", colour_lab="Age group")
}

# ---- |-Supplement: influenza ILI+ from each source/stream (do they agree?) ----
fig_iliplus_source_compare = function(timeseries_long){
  df = timeseries_long %>%
    filter(indicator=="ili_plus", pathogen=="Influenza", agegroup=="age_total",
           stream %in% c("ili_plus_sentinel", "ili_plus_nonsentinel", "ili_plus_respicompass"))
  facet_country_lines(df, colour_var="stream",
                      title="Influenza ILI+ by source: sentinel vs non-sentinel vs RespiCompass",
                      y_lab="ILI+ (Influenza)", colour_lab="Source / stream")
}

# ---- |-Supplement: detection counts by pathogen (sentinel) ----
fig_detections_dynamics = function(timeseries_long){
  df = timeseries_long %>%
    filter(indicator=="detections", stream=="typing_sentinel", agegroup=="age_total")
  facet_country_lines(df, colour_var="pathogen",
                      title="Sentinel detection counts by pathogen",
                      y_lab="Detections", colour_lab="Pathogen", colours=pathogen_colours)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### eyeballing(): assemble the figure manifest ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Returns list(meta, figures); each figure = list(section, title, subtitle, bullets, plot).
# section is one of "quality" / "dynamics" / "supplement" so the report can group them.
# countries / seasons default to everything present in the data.
eyeballing = function(models_in, params=NULL, data=NULL, countries=NULL, seasons=NULL){
  long    = models_in$data_timeseries_long
  summary = models_in$data_season_summary
  if (is.null(countries)) countries = sort(unique(long$country_short))
  if (is.null(seasons))   seasons   = sort(unique(long$season))

  long    = long    %>% filter(country_short %in% countries, season %in% seasons)
  summary = summary %>% filter(country_short %in% countries, season %in% seasons)

  figures = list(

    # ---- quality figures ----
    completeness = list(
      section  = "quality",
      title    = "Season-window completeness across countries and seasons",
      subtitle = paste("Definition: weeks with a non-missing value / number of ISO (Wednesday) weeks in the full Aug-Jul season.",
                       "A reported 0 counts as observed; positivity needs both detections and tests.",
                       "Caveat: the off-season is in the denominator, so winter-only streams cap well below 100% even when fully reported - compare with the active-span version in the supplement."),
      bullets  = c(
        "ILI/ARI consultation rates are the most complete streams; most countries report the large majority of in-season weeks recently. *(placeholder — edit)*",
        "Typing positivity is patchier and uneven across countries; several show large gaps in one or more seasons. *(placeholder — edit)*",
        "Pre-2020 seasons carry no SARS-CoV-2 or RSV typing — those blanks are expected, not missing data. *(placeholder — edit)*",
        "Values plateau near ~45-50% for cleanly-reported winter streams: that is the off-season denominator, not poor reporting. *(placeholder — edit)*"
      ),
      plot     = fig_completeness(summary)
    ),

    testing_volume = list(
      section  = "quality",
      title    = "Typing test volume by country and season",
      subtitle = "Definition: total specimens tested per season (sum of weekly sentinel/non-sentinel tests), log colour scale. Zero/absent = grey. A proxy for how trustworthy the positivity (and hence ILI+) estimates are.",
      bullets  = c(
        "Sentinel testing volumes are typically smaller and noisier than non-sentinel. *(placeholder — edit)*",
        "Very low season totals make positivity (and hence ILI+) unstable for those country/seasons. *(placeholder — edit)*",
        "Large cross-country differences in volume reflect surveillance system size, not just disease burden. *(placeholder — edit)*"
      ),
      plot     = fig_testing_volume(long)
    ),

    # ---- temporal-dynamics figures ----
    iliplus_dynamics = list(
      section  = "dynamics",
      title    = "ILI+ dynamics by pathogen",
      subtitle = "ILI consultation rate x sentinel positivity, one panel per country, continuous time (free y-axis per country).",
      bullets  = c(
        "Influenza ILI+ shows the classic sharp winter wave in most countries. *(placeholder — edit)*",
        "SARS-CoV-2 ILI+ appears from 2020/21 with less regular, multi-peak timing. *(placeholder — edit)*",
        "RSV ILI+ tends to lead the influenza peak by several weeks where both are observed. *(placeholder — edit)*",
        "Flat or absent curves indicate the country/season lacks ILI or positivity data (see the quality figures). *(placeholder — edit)*"
      ),
      plot     = fig_iliplus_dynamics(long)
    ),

    iliplus_aligned = list(
      section  = "dynamics",
      title    = "ILI+ timing aligned across countries",
      subtitle = "Same ILI+ series, but countries stacked in a single column on one shared time axis, lightly smoothed with a centred 3-week rolling mean, and each country x pathogen curve rescaled to its own peak (=1). This removes week-to-week noise and magnitude so only wave timing remains: a small RSV wave and a large influenza wave both reach the top of the panel.",
      bullets  = c(
        "Read vertically: a near-vertical alignment of peaks means countries waved synchronously that season. *(placeholder — edit)*",
        "Look for west-to-east or north-to-south lags in the influenza peak. *(placeholder — edit)*",
        "RSV (orange) consistently sitting left of influenza (green) confirms its earlier seasonal onset. *(placeholder — edit)*"
      ),
      plot     = fig_iliplus_aligned(long)
    ),

    positivity_dynamics = list(
      section  = "dynamics",
      title    = "Test positivity by pathogen",
      subtitle = "Sentinel positivity (detections / tests) per pathogen, one panel per country.",
      bullets  = c(
        "Positivity isolates pathogen circulation timing independent of consultation behaviour. *(placeholder — edit)*",
        "Influenza and RSV positivity are strongly seasonal; SARS-CoV-2 is more persistent year-round. *(placeholder — edit)*",
        "Spiky positivity in low-volume weeks (see test volume) should be read with caution. *(placeholder — edit)*"
      ),
      plot     = fig_positivity_dynamics(long)
    ),

    syndromic_dynamics = list(
      section  = "dynamics",
      title    = "Syndromic consultation rates: ILI vs ARI",
      subtitle = "Raw ERVISS consultation rates, one panel per country.",
      bullets  = c(
        "ARI sits above ILI by construction (broader case definition). *(placeholder — edit)*",
        "The ILI/ARI gap and its seasonality varies by country, hinting at differing case definitions. *(placeholder — edit)*",
        "These raw rates are the syndromic backbone every ILI+ series is built on. *(placeholder — edit)*"
      ),
      plot     = fig_syndromic_dynamics(long)
    ),

    # ---- supplement figures (recommended extras) ----
    completeness_active = list(
      section  = "supplement",
      title    = "Active-span completeness across countries and seasons",
      subtitle = paste("Definition: weeks with a non-missing value / number of weeks between the first and last reported week (inclusive).",
                       "Unlike season-window completeness this ignores the off-season and instead measures gaps WITHIN each stream's active reporting period: ~100% means no mid-season holes."),
      bullets  = c(
        "Most consultation-rate streams approach 100% here — confirming the low season-window values are off-season, not gaps. *(placeholder — edit)*",
        "Streams that stay low even on this measure have genuine mid-season reporting holes. *(placeholder — edit)*"
      ),
      plot     = fig_completeness_active(summary)
    ),

    iliplus_source_compare = list(
      section  = "supplement",
      title    = "Influenza ILI+ by source (sentinel vs non-sentinel vs RespiCompass)",
      subtitle = "Three independent constructions of influenza ILI+ overlaid per country; close agreement is reassuring, divergence flags a source/definition issue.",
      bullets  = c(
        "The sentinel ILI+ should track the RespiCompass ILI+ closely (same construction). *(placeholder — edit)*",
        "Non-sentinel ILI+ can differ in level where its testing population differs. *(placeholder — edit)*"
      ),
      plot     = fig_iliplus_source_compare(long)
    ),

    age_dynamics = list(
      section  = "supplement",
      title    = "Age-specific ILI consultation rate",
      subtitle = "ERVISS ILI consultation rate by age group, one panel per country (age_total excluded).",
      bullets  = c(
        "Young children (0-4, 5-14) usually carry the highest ILI consultation rates. *(placeholder — edit)*",
        "Age coverage is itself uneven across countries — some report only the total. *(placeholder — edit)*"
      ),
      plot     = fig_age_dynamics(long)
    ),

    detections_dynamics = list(
      section  = "supplement",
      title    = "Sentinel detection counts by pathogen",
      subtitle = "Weekly sentinel detections per pathogen, one panel per country (counts, not rates).",
      bullets  = c(
        "Detection counts combine circulation and testing effort; pair with positivity and test volume. *(placeholder — edit)*",
        "Useful for spotting which pathogen dominated a given season in absolute terms. *(placeholder — edit)*"
      ),
      plot     = fig_detections_dynamics(long)
    )
  )

  list(
    meta    = list(countries=countries, seasons=seasons,
                   n_countries=length(countries), n_seasons=length(seasons)),
    figures = figures
  )
}
