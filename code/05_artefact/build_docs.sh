#!/usr/bin/env bash
# ==============================================================================
# Build the two documents: the RespiCast decision note (PDF + Word) and the
# survey supplement (PDF), from their committed sources.
# ==============================================================================
#   decision note : code/05_artefact/respicast_decision_note.md   (single source)
#                     -> pandoc  -> output/respicast_decision_note.docx
#                     -> pandoc + Chromium print-to-pdf -> ...note.pdf
#   supplement    : code/05_artefact/respicast_supplement.html + the five figures
#                     -> build_supplement.R -> self-contained HTML -> PDF
#
# LibreOffice cannot convert docx->pdf in this sandbox, so both formats are
# rendered from the SAME markdown/HTML instead -- that keeps them identical.
#
# Run from the repo root (after the figures exist):
#   Rscript code/05_figures/fig_respicast_value.R
#   Rscript code/05_figures/fig_respicast_decision.R
#   ./code/05_artefact/build_docs.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

# ---- locate a Chromium that can print to PDF (path is version-stamped, so glob) ----
CHROME="${CHROME_BIN:-}"
if [ -z "$CHROME" ]; then
  for c in /opt/pw-browsers/chromium_headless_shell-*/chrome-linux/headless_shell \
           /opt/pw-browsers/chromium-*/chrome-linux/chrome \
           "$(command -v chromium || true)" "$(command -v google-chrome || true)"; do
    [ -n "$c" ] && [ -x "$c" ] && { CHROME="$c"; break; }
  done
fi
[ -n "$CHROME" ] || { echo "ERROR: no Chromium found; set CHROME_BIN=/path/to/chrome" >&2; exit 1; }
echo "chromium: $CHROME"

topdf() {  # topdf <input.html> <output.pdf>
  "$CHROME" --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
            --print-to-pdf="$2" "file://$PWD/$1" 2>/dev/null
}

# ---- 1. decision note: one markdown source -> .docx and .pdf ----
NOTE=code/05_artefact/respicast_decision_note.md
pandoc "$NOTE" -o output/respicast_decision_note.docx
pandoc "$NOTE" --embed-resources --standalone \
       -H code/05_artefact/decision_note_style.html \
       -o output/respicast_decision_note.html
topdf output/respicast_decision_note.html output/respicast_decision_note.pdf

# ---- 2. survey supplement: embed the figures, then print ----
Rscript code/05_artefact/build_supplement.R
topdf output/respicast_value_supplement.html output/respicast_value_supplement.pdf

echo
echo "documents ->"
ls -la output/respicast_decision_note.{pdf,docx} output/respicast_value_supplement.pdf
