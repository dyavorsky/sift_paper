#!/bin/sh
# Render every comparison notebook. Run from the repo root:
#   sh comparisons/render_all.sh
# Experiments are NOT re-run -- do that first if you want fresh numbers:
#   for f in comparisons/exp_*.R; do Rscript "$f"; done
for f in comparisons/index.qmd comparisons/compare-*.qmd; do
  echo "--- $f"
  quarto render "$f" --to html || exit 1
done
