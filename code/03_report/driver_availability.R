# driver_availability.R
#
# Companion availability panel to the ILI+ data-availability plot (code/03_report/data_availability.R):
# shows which country-seasons carry the three externally-pulled candidate drivers -- dominant subtype,
# vaccine effectiveness (VE), and 65+ vaccination coverage. Designed to sit DIRECTLY BELOW the ILI+
# panel: same 12-season chronological x-axis, same COVID greying, same RespiCompass->ERVISS handover line.
#
# Key point made visually: subtype & VE are CONTINENTAL / season-level (one "EU/EEA" row, identical
# across countries), whereas 65+ coverage is COUNTRY-specific (a full country x season grid with gaps).
# Fill encodes availability + provenance (already in repo / newly pulled / interim-only / missing / COVID).
# Run from the repo root:  Rscript code/03_report/driver_availability.R
suppressMessages({library(ggplot2)})

seasons_all <- c("2014/2015","2015/2016","2016/2017","2017/2018","2018/2019",
                 "2019/2020","2020/2021","2021/2022","2022/2023",
                 "2023/2024","2024/2025","2025/2026")          # same axis as the ILI+ panel
covid <- c("2019/2020","2020/2021","2021/2022","2022/2023")
pre   <- c("2014/2015","2015/2016","2016/2017","2017/2018","2018/2019")
post  <- c("2023/2024","2024/2025","2025/2026")
countries <- sort(unique(read.csv("data/slim_flu_iliplus.csv", stringsAsFactors=FALSE)$country_short))

# ---- |-subtype (continental: literature pre-COVID, ERVISS post-COVID) ----
sub <- read.csv("data/external/dominant_subtype_by_season.csv", stringsAsFactors=FALSE)
subtype_df <- data.frame(driver="Dominant\nsubtype", row="EU/EEA", season=sub$season,
  status=ifelse(sub$basis=="continental_literature","newly pulled","already in repo"), stringsAsFactors=FALSE)

# ---- |-vaccine effectiveness (continental; interim-only seasons flagged) ----
ve <- read.csv("data/external/vaccine_effectiveness.csv", stringsAsFactors=FALSE)
ve_s <- aggregate(timing ~ season, ve, function(x) if (any(x=="end_of_season")) "newly pulled" else "newly pulled (interim)")
ve_df <- data.frame(driver="Vaccine\neffectiveness", row="EU/EEA", season=ve_s$season, status=ve_s$timing, stringsAsFactors=FALSE)

# ---- |-65+ coverage (country x season: in-repo pre-COVID + externally pulled post-COVID) ----
dv <- read.csv("output/descriptors_vax.csv", stringsAsFactors=FALSE)
cov_pre <- subset(dv, is.finite(vax_cov_65) & season %in% pre, select=c(country, season)); cov_pre$status <- "already in repo"
cp <- read.csv("data/external/vaccination_coverage_65plus_postcovid.csv", stringsAsFactors=FALSE)
cov_post <- subset(cp, panel_country %in% c(TRUE,"TRUE") & season %in% post, select=c(country_short, season))
names(cov_post)[1] <- "country"; cov_post$status <- "newly pulled"; cov_post <- unique(cov_post)
have <- rbind(cov_pre, cov_post)
grid <- merge(expand.grid(country=countries, season=c(pre,post), stringsAsFactors=FALSE), have,
              by=c("country","season"), all.x=TRUE)
grid$status[is.na(grid$status)] <- "not available"
coverage_df <- data.frame(driver="65+ vaccination\ncoverage", row=grid$country, season=grid$season,
                          status=grid$status, stringsAsFactors=FALSE)

# ---- |-COVID-excluded cells (grey block, matches the ILI+ panel) ----
covid_df <- rbind(
  data.frame(driver="Dominant\nsubtype",       row="EU/EEA", season=covid, status="excluded: COVID"),
  data.frame(driver="Vaccine\neffectiveness",   row="EU/EEA", season=covid, status="excluded: COVID"),
  data.frame(driver="65+ vaccination\ncoverage", row=rep(countries, each=length(covid)),
             season=rep(covid, times=length(countries)), status="excluded: COVID"))

d <- rbind(subtype_df, ve_df, coverage_df, covid_df)
d$driver <- factor(d$driver, levels=c("Dominant\nsubtype","Vaccine\neffectiveness","65+ vaccination\ncoverage"))
d$season <- factor(d$season, levels=seasons_all)
d$row    <- factor(d$row, levels=c(rev(sort(countries)), "EU/EEA"))
d$status <- factor(d$status, levels=c("already in repo","newly pulled","newly pulled (interim)","not available","excluded: COVID"))
cols <- c("already in repo"="#1b9e77","newly pulled"="#d95f02","newly pulled (interim)"="#e6ab02",
          "not available"="grey92","excluded: COVID"="grey65")
xdiv <- match("2023/2024", seasons_all) + 0.5

p <- ggplot(d, aes(season, row, fill=status)) +
  geom_tile(color="white", linewidth=0.5, width=0.95, height=0.9) +
  geom_vline(xintercept=xdiv, linetype="dashed", color="grey20") +
  facet_grid(driver ~ ., scales="free_y", space="free_y") +
  scale_fill_manual(name=NULL, values=cols, drop=FALSE) +
  labs(title="External driver data availability (companion to the ILI+ panel)",
       subtitle="Subtype & VE are continental (season-level -> one 'EU/EEA' row, same for all countries); 65+ coverage is country-specific.\nDashed line = RespiCompass->ERVISS handover; grey = COVID seasons excluded.",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(),
        strip.text.y=element_text(angle=0, face="bold"), legend.position="top")

dir.create("output", showWarnings=FALSE)
ggsave("output/driver_availability.png", p, width=12, height=8.5, dpi=110)
cat("figure -> output/driver_availability.png\n")

# brief coverage summary to stdout
cat(sprintf("coverage cells: in-repo(pre)=%d, newly-pulled(post)=%d, missing(non-COVID)=%d\n",
            sum(coverage_df$status=="already in repo"), sum(coverage_df$status=="newly pulled"),
            sum(coverage_df$status=="not available")))
