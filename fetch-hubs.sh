#!/usr/bin/env bash
# ==============================================================================
# Clone the five forecasting hubs needed for a FULL re-derivation from raw data.
# ==============================================================================
# Only needed to regenerate output/artefact_data.json from source. The dashboard
# and all documents rebuild from committed data without this (./reproduce.sh page),
# which is why the SessionStart hook deliberately does NOT run it: ~4 GB of clones
# is far too slow for session start.
#
#   ./fetch-hubs.sh              -> clones into ../hubs
#   RESPICAST_HUBS_DIR=/data ./fetch-hubs.sh   -> clones somewhere else
#
# Shallow clones (--depth 1) are enough: the current tree already contains every
# historical submission. Re-running is safe -- existing clones are skipped.
set -euo pipefail

HUBS_DIR="${RESPICAST_HUBS_DIR:-$(cd "$(dirname "$0")/.." && pwd)/hubs}"
ORG=https://github.com/european-modelling-hubs
REPOS=(
  covid19-forecast-hub-europe_archive   # COVID cases / hospitalisations / deaths, 2021-2024 (legacy format, ~2.4 GB)
  RespiCast-Covid19                     # COVID-19 hospitalisations, 2024-
  RespiCast-SyndromicIndicators         # ILI + ARI incidence, 2024-
  flu-forecast-hub_archive              # ILI, 2023/24 season
  ari-forecast-hub_archive              # ARI, 2023/24 season
)

mkdir -p "$HUBS_DIR"
echo "cloning into $HUBS_DIR (this pulls several GB and takes a while)"
for r in "${REPOS[@]}"; do
  if [ -d "$HUBS_DIR/$r/.git" ]; then
    echo "  [skip] $r already cloned"
  else
    echo "  [clone] $r"
    git clone --depth 1 "$ORG/$r.git" "$HUBS_DIR/$r"
  fi
done

echo
echo "done. Now run the full rebuild:  ./reproduce.sh"
[ -n "${RESPICAST_HUBS_DIR:-}" ] && echo "(remember to keep RESPICAST_HUBS_DIR set: $HUBS_DIR)"
