#!/bin/bash
# SessionStart hook for Claude Code on the web: make tests + report rendering work.
#
# This repo uses renv, but CRAN is unreachable under the web network policy, so
# renv::restore() cannot fetch packages. The R packages instead ship in the apt base image,
# so we install the few system libraries the pipeline needs and populate the renv project
# library from the system packages with renv::hydrate() (fully offline). Idempotent.
set -euo pipefail

# Only run in the remote (web) environment; locally the user manages their own renv.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

export DEBIAN_FRONTEND=noninteractive

# System libraries: libmagick++-dev (summarytools), pandoc (renders the eyeballing report),
# and the apt builds of the packages the pipeline + tests need.
apt-get update -qq
apt-get install -y -qq libmagick++-dev pandoc r-cran-renv r-cran-testthat r-cran-emayili

# Populate the renv project library (renv/library/, gitignored) from the installed system
# packages, so renv is active and functional offline. The data pipeline and the full test
# suite work from this. NOTE: a few plotting-only packages (ggdist, arrow) are not in the
# apt repos and cannot be fetched offline (CRAN is blocked here); rendering the eyeballing
# report's plots needs them installed manually if/when CRAN is reachable.
Rscript -e 'options(renv.consent = TRUE); renv::hydrate()'
