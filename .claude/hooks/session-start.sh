#!/bin/bash
# ==============================================================================
# SessionStart hook — make this repo runnable in a fresh Claude Code container.
# ==============================================================================
# The base image ships R 4.3, pandoc, LibreOffice and Chromium, but NO R packages,
# so without this hook `Rscript code/00_main.R` and every build script fail on the
# first library() call.
#
# Two things worth knowing about this environment (both verified, not assumed):
#   1. Direct CRAN (cloud.r-project.org) is NOT reachable through the agent proxy,
#      so plain install.packages() fails. Ubuntu's r-cran-* apt binaries are both
#      the fastest and the reliable route; Posit's P3M binary repo IS reachable and
#      is used only as a fallback for anything apt is missing.
#   2. The forecasting-hub clones (~4 GB) are NOT fetched here -- far too slow for
#      session start. The committed output/artefact_data.json means the dashboard
#      and every document still rebuild without them (`./reproduce.sh page`).
#      Only a full re-derivation from raw hub data needs them; see README.md.
#
# Idempotent: on resume/clear/compact it detects the packages and exits in ~1s.
set -euo pipefail

# Local machines already have their own R setup -- only touch the cloud container.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# The 13 packages code/01_support/setup.R loads (magick/pdftools are NOT needed --
# they were only ever used ad hoc, and they drag in heavy system libraries).
R_PKGS=(readxl dplyr tidyr readr stringr purrr lubridate data.table jsonlite ggplot2 patchwork scales here)

# Dev-only: lintr powers `Rscript -e 'lintr::lint("file.R")'`. Installed best-effort --
# a failure here must never block a session, so it is kept out of the hard check below.
DEV_PKGS=(lintr)

# ---- fast path: everything already present (resume / clear / compact) ----
# Checks the dev package too, so a container that has the runtime but not lintr
# still gets topped up rather than skipped.
if Rscript -e 'q(status = as.integer(!all(sapply(commandArgs(TRUE), requireNamespace, quietly = TRUE))))' \
     "${R_PKGS[@]}" "${DEV_PKGS[@]}" >/dev/null 2>&1; then
  echo "[session-start] R packages already installed - nothing to do."
  exit 0
fi

SUDO=""
[ "$(id -u)" -ne 0 ] && SUDO="sudo"
export DEBIAN_FRONTEND=noninteractive

echo "[session-start] Installing R packages (apt binaries; no compilation)..."
$SUDO apt-get update -qq || true   # stale lists shouldn't abort the session

# apt names are lower-case: r-cran-data.table, r-cran-ggplot2, ...
APT_PKGS=()
for p in "${R_PKGS[@]}" "${DEV_PKGS[@]}"; do APT_PKGS+=("r-cran-${p,,}"); done
$SUDO apt-get install -y -qq --no-install-recommends "${APT_PKGS[@]}" || true

# ---- fallback: anything apt couldn't supply, take from the P3M binary repo ----
# (binary, so still no compilation; CRAN proper is unreachable behind the proxy)
MISSING=$(Rscript -e 'cat(paste(Filter(function(p) !requireNamespace(p, quietly = TRUE), commandArgs(TRUE)), collapse = " "))' \
            "${R_PKGS[@]}" 2>/dev/null || true)
if [ -n "${MISSING// /}" ]; then
  echo "[session-start] apt missed:$MISSING - fetching binaries from Posit P3M..."
  Rscript -e 'pkgs <- commandArgs(TRUE); install.packages(pkgs, repos = "https://packagemanager.posit.co/cran/__linux__/noble/latest", quiet = TRUE)' \
    $MISSING || true
fi

# ---- verify: fail loudly rather than let the session start half-broken ----
if ! Rscript -e 'bad <- Filter(function(p) !requireNamespace(p, quietly = TRUE), commandArgs(TRUE));
                 if (length(bad)) { cat("MISSING:", paste(bad, collapse = ", "), "\n"); quit(status = 1) }' \
     "${R_PKGS[@]}"; then
  echo "[session-start] ERROR: some R packages are still missing (see above)." >&2
  exit 1
fi

echo "[session-start] Ready: R packages installed. Rebuild the dashboard with ./reproduce.sh page"
