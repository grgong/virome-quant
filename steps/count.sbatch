#!/usr/bin/env bash
# Counts per vOTU (CoverM), detection, result tables and heatmap (count_viral.py, plot_heatmap.py).
# Outputs: results/{coverm_contig.tsv, viral_long.tsv, votus.tsv, viral_rpm.tsv, viral_rpm.pdf/png}
source "${VP}/steps/common.sh"
python3 -I "${VP}/bin/count_viral.py" "${OUT}"
python3 -I "${VP}/bin/plot_heatmap.py" "${OUT}"
touch "${OUT}/done/count"
log "done"
