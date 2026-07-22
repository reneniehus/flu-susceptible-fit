#!/usr/bin/env bash
# ============================================================================
# Reproduce the interactive dashboard (and the static figures) end to end.
#
#   ./reproduce.sh          full rebuild: data -> artefact JSON -> dashboard + figures
#   ./reproduce.sh page     ONLY rebuild artefact/dashboard.html from the committed
#                           output/artefact_data.json (no hub clones needed)
#
# The committed output/artefact_data.json means the dashboard reproduces WITHOUT the
# ~4 GB of hub submissions; the full rebuild regenerates that JSON from source, and
# needs the five hubs cloned into ../hubs (or $RESPICAST_HUBS_DIR) -- see README.md.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"

if [[ "${1:-}" == "page" ]]; then
  Rscript code/05_artefact/build_pages.R
  echo "dashboard -> artefact/dashboard.html"
  exit 0
fi

Rscript code/00_main.R                    # survey + 5 hubs -> output/*.csv + artefact_data.json
Rscript code/05_artefact/build_pages.R    # artefact_data.json + template -> artefact/dashboard.html
Rscript code/05_figures/fig_survey.R      # -> output/figures/survey.png
Rscript code/05_figures/fig_coverage.R    # -> output/figures/coverage.png
echo
echo "dashboard -> artefact/dashboard.html   (open in a browser, or publish as an Artifact)"
