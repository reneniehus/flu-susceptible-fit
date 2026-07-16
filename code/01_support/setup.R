# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Session setup: libraries, console helpers, shared plotting palette ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sourced first by code/00_main.R. Keeps the noisy library() calls, the dplyr
# re-masking and the project-wide colour palette in one place so every script
# downstream starts from the same, quiet, predictable session.

# ---- |-quiet library() ----
# wrap base::library so a dozen attach messages don't drown the run log
library <- function(...) suppressPackageStartupMessages(base::library(..., quietly = TRUE))

library(readxl)       # read the .xlsx survey export
library(dplyr)        # data wrangling grammar
library(tidyr)        # pivot / unnest
library(readr)        # fast typed CSV read/write
library(stringr)      # string helpers (str_detect, str_squish, ...)
library(purrr)        # map_* over the many forecast files
library(lubridate)    # ISO-week dates on the forecast timeline
library(data.table)   # fread(): the ~2,100 forecast CSVs load in seconds
library(jsonlite)     # write the artefact-ready JSON
library(ggplot2)      # static figure companions to the HTML artefacts
library(patchwork)    # compose multi-panel figures ( `/` and `|` operators )
library(scales)       # axis formatting for the figures
library(here)         # repo-root-relative paths, so scripts run from anywhere

# ---- |-re-mask the verbs the tidyverse loses to data.table / stats ----
# data.table also exports first()/last()/between(); pin the dplyr ones we rely on.
select  <- dplyr::select
filter  <- dplyr::filter
between <- dplyr::between
first   <- dplyr::first
last    <- dplyr::last

# ---- |-console helpers ----
# tiny status printer so each stage announces itself the same way (no colour dep)
say  <- function(...) cat(paste0("  ", ..., "\n"))
step <- function(...) cat(paste0("\n== ", ..., " ==\n"))

# max() over a possibly-empty vector, without the "-Inf / no non-missing" warning
safe_max <- function(x) if (length(x)) max(x) else NA_integer_

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Shared palette ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ONE colour system for the whole project, so the R figures and the HTML artefacts
# read as siblings. These are the validated categorical/sequential hues from the
# data-viz design system (light-surface steps); the artefact re-declares the same
# hex as CSS custom properties. Assign categorical hues in FIXED order, never cycled.

palette_cat <- c(                       # categorical identity slots, in safety order
  blue    = "#2a78d6",
  green   = "#008300",
  magenta = "#e87ba4",
  yellow  = "#eda100",
  aqua    = "#1baf7a",
  orange  = "#eb6834",
  violet  = "#4a3aa7",
  red     = "#e34948"
)

# indicator -> colour, fixed so COVID / ILI / ARI keep their identity across every figure
palette_indicator <- c(
  "COVID-19 hospitalisations" = "#2a78d6",   # blue  (RespiCast-Covid19 hub)
  "ILI incidence"             = "#eb6834",   # orange
  "ARI incidence"             = "#4a3aa7"    # violet
)

# the ordered Likert ramp (diverging: disagree/unlikely <-> agree/likely, gray neutral)
palette_likert <- c(
  "Very unlikely"              = "#184f95",   # ends use the blue<->red diverging poles
  "Unlikely"                   = "#6da7ec",
  "Unsure"                     = "#c3c2b7",   # neutral gray midpoint
  "Likely"                     = "#eb8c6a",
  "Very likely"                = "#c0392b"
)
palette_agree <- c(
  "Strongly disagree"          = "#184f95",
  "Disagree"                   = "#6da7ec",
  "Neither agree nor disagree" = "#c3c2b7",
  "Agree"                      = "#eb8c6a",
  "Strongly agree"             = "#c0392b"
)

# a light ggplot theme reused by every figure (recessive grid, no chart junk)
theme_project <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, colour = "#e1e0d9"),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "#52514e"),
      axis.title       = element_text(colour = "#52514e"),
      legend.position  = "top"
    )
}
