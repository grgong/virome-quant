#!/usr/bin/env bash
# Counts per vOTU (CoverM), detection and result tables (bin/count_viral.py).
# Outputs: results/{coverm_contig.tsv, viral_long.tsv, votus.tsv, viral_rpm.tsv}
source "${VP}/steps/common.sh"
python3 -I "${VP}/bin/count_viral.py" "${OUT}"
touch "${OUT}/done/count"
log "done"
