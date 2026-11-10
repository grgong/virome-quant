"""Check the samplesheet (and the vOTU panel) and write the manifest the steps read.

Usage: python3 -I check_inputs.py <samples.tsv> <manifest.tsv> [<panel_dir>]
samples.tsv: tab-separated, header with at least sample, fastq_1, fastq_2, flowcell_lane; any other
column (stage, tissue, ...) is metadata carried to the result tables. Several lanes of one sample:
comma-separated files in fastq_1 / fastq_2 (same order) and flowcell:lane items in flowcell_lane.
A relative FASTQ path is relative to the samplesheet's directory.
manifest.tsv: sample, fastq_1, fastq_2 (absolute paths, no header), in samplesheet order.
panel_dir: panel.fa, virus_units.tsv and quality_summary.tsv, which must name the same vOTUs.
"""
import csv
import os
import re
import sys

src, dst = sys.argv[1:3]
panel = sys.argv[3] if len(sys.argv) > 3 else None
base = os.path.dirname(os.path.abspath(src))
if "\r" in open(src, newline="").read():
    sys.exit(f"ERROR: {src} has Windows line endings; convert it: sed -i 's/\\r$//' {src}")
rows = list(csv.DictReader(open(src, newline=""), delimiter="\t"))
errors = []
if not rows:
    sys.exit(f"ERROR: {src} has no sample rows")
missing = {"sample", "fastq_1", "fastq_2", "flowcell_lane"} - set(rows[0])
if missing:
    sys.exit(f"ERROR: {src} lacks the column(s) {', '.join(sorted(missing))} (tab-separated?)")
seen, manifest = set(), []
for i, r in enumerate(rows, start=2):
    if any(v is None for v in r.values()) or None in r:
        errors.append(f"line {i}: number of columns differs from the header")
        continue
    s = r["sample"]
    if not re.fullmatch(r"[A-Za-z0-9._-]+", s):
        errors.append(f"line {i}: sample name '{s}' (letters, digits, . _ - only)")
    if s in seen:
        errors.append(f"line {i}: sample '{s}' twice")
    seen.add(s)
    f1, f2 = ([os.path.join(base, f) for f in r[c].split(",")] for c in ("fastq_1", "fastq_2"))
    if len(f1) != len(f2):
        errors.append(f"line {i}: {len(f1)} fastq_1 files but {len(f2)} fastq_2 files")
    for f in f1 + f2:
        if re.search(r"\s", f):
            errors.append(f"line {i}: space in the path {f}")
        elif not os.path.isfile(f) or os.path.getsize(f) == 0:
            errors.append(f"line {i}: missing or empty file {f}")
    if not all(re.fullmatch(r"[^:,\s]+:[^:,\s]+", x) for x in r["flowcell_lane"].split(",")):
        errors.append(f"line {i}: flowcell_lane '{r['flowcell_lane']}' (flowcell:lane, comma-separated; from the read headers)")
    manifest.append(f"{s}\t{','.join(f1)}\t{','.join(f2)}\n")

if panel:
    try:
        fa = {l[1:].split()[0] for l in open(f"{panel}/panel.fa") if l.startswith(">")}
        units = {r["ref"] for r in csv.DictReader(open(f"{panel}/virus_units.tsv"), delimiter="\t")
                 if r["is_representative"] == "1"}
        checkv = {r["contig_id"] for r in csv.DictReader(open(f"{panel}/quality_summary.tsv"), delimiter="\t")}
    except (OSError, KeyError) as e:
        sys.exit(f"ERROR: panel {panel}: {e} (needs panel.fa, virus_units.tsv, quality_summary.tsv)")
    if not fa or not (fa == units == checkv):
        errors.append(f"panel {panel}: panel.fa ({len(fa)}), virus_units.tsv ({len(units)} representatives) "
                      f"and quality_summary.tsv ({len(checkv)}) do not name the same vOTUs")

if errors:
    sys.exit("ERROR in the inputs:\n  " + "\n  ".join(errors))
with open(dst, "w") as fo:
    fo.writelines(manifest)
print(f"{len(rows)} samples, metadata columns: "
      + (", ".join(c for c in rows[0] if c not in ("sample", "fastq_1", "fastq_2", "flowcell_lane")) or "none")
      + (f"; panel: {len(fa)} vOTUs" if panel else ""))
