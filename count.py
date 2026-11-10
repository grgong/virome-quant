"""Count reads per vOTU in each library with CoverM.

    python3 count.py samples.tsv panel_dir outdir

Reads must be at least 90 % identical to the vOTU. A vOTU counts as present in a library when
those reads cover at least 75 % of it (Roux et al. 2017), otherwise its count is 0. RPM is per
million reads that pass fastp. Writes outdir/results/viral_long.tsv (every library and vOTU) and
viral_rpm.tsv (the high-quality vOTUs found in at least one library).
"""
import csv
import json
import os
import subprocess
import sys

MIN_ID, MIN_COVERED = 90, 0.75
samples_tsv, panel, out = sys.argv[1:]


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


samples = read_tsv(samples_tsv)
meta = [c for c in samples[0] if c not in ("sample", "fastq_1", "fastq_2")]
rep = {r["ref"]: r for r in read_tsv(f"{panel}/virus_units.tsv") if r["is_representative"] == "1"}
quality = {r["contig_id"]: r["miuvig_quality"] for r in read_tsv(f"{panel}/quality_summary.tsv")}
assert set(quality) == set(rep), "virus_units.tsv and quality_summary.tsv name different vOTUs"

bams = [f"{out}/map/{s['sample']}.bam" for s in samples]
coverm = subprocess.run(["coverm", "contig", "-b", *bams, "--min-read-percent-identity", str(MIN_ID),
                         "--contig-end-exclusion", "0", "-m", "count", "covered_fraction",
                         "--output-format", "sparse", "-t", "4"],
                        capture_output=True, text=True, check=True).stdout
reads, covered = {}, {}
for r in csv.DictReader(coverm.splitlines(), delimiter="\t"):
    key = (r["Sample"], r["Contig"])
    reads[key], covered[key] = int(r["Read Count"]), float(r["Covered Fraction"])

os.makedirs(f"{out}/results", exist_ok=True)
rpm = {}
with open(f"{out}/results/viral_long.tsv", "w") as f:
    f.write("\t".join(["sample", *meta, "votu", "display_name", "checkv_quality", "reads_id90",
                       "covered_fraction", "detected", "reads", "reads_total", "rpm"]) + "\n")
    for s in samples:
        sid = s["sample"]
        total = json.load(open(f"{out}/fastp/{sid}.json"))["summary"]["after_filtering"]["total_reads"]
        for u in sorted(rep):
            detected = covered[sid, u] >= MIN_COVERED
            n = reads[sid, u] if detected else 0
            rpm[sid, u] = n / total * 1e6
            f.write("\t".join([sid, *(s[c] for c in meta), u, rep[u]["display_name"], quality[u],
                               str(reads[sid, u]), f"{covered[sid, u]:.4f}", str(int(detected)), str(n),
                               str(total), f"{rpm[sid, u]:.3f}"]) + "\n")

hq = [u for u in sorted(rep) if quality[u] == "High-quality" and any(rpm[s["sample"], u] for s in samples)]
with open(f"{out}/results/viral_rpm.tsv", "w") as f:
    f.write("\t".join(["sample", *meta, *(f"{rep[u]['display_name']} [{u}]" for u in hq)]) + "\n")
    for s in samples:
        f.write("\t".join([s["sample"], *(s[c] for c in meta), *(f"{rpm[s['sample'], u]:.3f}" for u in hq)]) + "\n")
print(f"{len(rep)} vOTUs, {len(hq)} high-quality ones found in at least one library")
