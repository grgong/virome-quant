"""Reads per vOTU in every sample (CoverM), detection, and the result tables.

Usage: python3 -I count_viral.py <run_dir>

Inputs (under run_dir): samples.tsv, panel/{virus_units.tsv, quality_summary.tsv} (from the discovery side),
map/<sample>.bam and .nreads
Reads: primary alignments with >= 80 % of the trimmed read aligned (map step) and >= 90 % identity
(CoverM --min-read-percent-identity 90). A vOTU is detected in a sample when >= 75 % of its
reference is covered by these reads (CoverM covered_fraction, --contig-end-exclusion 0);
otherwise its reads count 0 (Roux et al. 2017, PeerJ 5:e3817). RPM = counted reads per million
reads passing fastp, both mates counted.
Outputs (results/):
  coverm_contig.tsv   CoverM output (sample x vOTU representative)
  viral_long.tsv      sample x vOTU, every vOTU: the sample's metadata columns, vOTU name and CheckV
                      tier, reads >= 90 % identity, covered fraction, detected, counted reads, RPM,
                      the highest RPM of the vOTU in another sample, that sample, and this
                      sample's RPM as % of it
  votus.tsv           one row per vOTU: name, length, contigs, CheckV tier and completeness,
                      samples detected, highest RPM
  viral_rpm.tsv       samples (rows, with metadata) x MIUViG high-quality vOTUs detected in at
                      least one sample (columns): RPM
Descriptive only, no test.
"""
import collections
import csv
import os
import subprocess
import sys

MIN_ID, MIN_COVERED = 90, 0.75
RUN = sys.argv[1]
OUTD = f"{RUN}/results"
os.makedirs(OUTD, exist_ok=True)

samples = list(csv.DictReader(open(f"{RUN}/samples.tsv"), delimiter="\t"))
meta_cols = [c for c in samples[0] if c not in ("sample", "fastq_1", "fastq_2")]
units = list(csv.DictReader(open(f"{RUN}/panel/virus_units.tsv"), delimiter="\t"))
rep = {r["ref"]: r for r in units if r["is_representative"] == "1"}
votus = sorted(rep)
n_contigs = collections.Counter(r["unit"] for r in units)
checkv = {r["contig_id"]: r for r in csv.DictReader(open(f"{RUN}/panel/quality_summary.tsv"), delimiter="\t")}
assert set(checkv) == set(votus), "virus_units.tsv and quality_summary.tsv name different vOTUs"

bams = [f"{RUN}/map/{s['sample']}.bam" for s in samples]
cov = subprocess.run(["coverm", "contig", "-b", *bams, "--min-read-percent-identity", str(MIN_ID),
                      "--contig-end-exclusion", "0", "-m", "count", "covered_fraction",
                      "--output-format", "sparse", "-t", os.environ.get("SLURM_CPUS_PER_TASK", "4")],
                     capture_output=True, text=True, check=True).stdout
with open(f"{OUTD}/coverm_contig.tsv", "w") as fo:
    fo.write(cov)
reads_id, covered = {}, {}
for r in csv.DictReader(cov.splitlines(), delimiter="\t"):
    key = (r["Sample"], r["Contig"])
    reads_id[key], covered[key] = int(r["Read Count"]), float(r["Covered Fraction"])
for s in samples:
    s["total"] = int(list(csv.DictReader(open(f"{RUN}/map/{s['sample']}.nreads"), delimiter="\t"))[0]["reads_total"])
    assert all((s["sample"], u) in reads_id for u in votus), f"CoverM output lacks {s['sample']}"

counted = {k: n if covered[k] >= MIN_COVERED else 0 for k, n in reads_id.items()}
by_sample = {s["sample"]: s for s in samples}
rpm = {(sid, u): n / by_sample[sid]["total"] * 1e6 for (sid, u), n in counted.items()}
by_unit = collections.defaultdict(list)  # vOTU -> [(rpm, sample)], highest first
for (sid, u), v in rpm.items():
    if v:
        by_unit[u].append((v, sid))
for vs in by_unit.values():
    vs.sort(reverse=True)


def highest_other(sid, u):
    return next(((v, o) for v, o in by_unit[u] if o != sid), (0.0, ""))


def name(u):
    return rep[u]["display_name"]


with open(f"{OUTD}/viral_long.tsv", "w") as fo:
    fo.write("\t".join(["sample", *meta_cols, "votu", "display_name", "checkv_quality", "reads_id90",
                        "covered_fraction", "detected", "reads", "reads_total", "rpm", "max_rpm_other_sample",
                        "max_other_sample", "pct_of_max_other"]) + "\n")
    for s in samples:
        sid = s["sample"]
        for u in votus:
            k, r = (sid, u), rpm[(sid, u)]
            m, mo = highest_other(sid, u)
            cols = [f"{m:.3f}", mo, f"{100 * r / m:.4f}" if m else "NA"]
            fo.write("\t".join([sid, *(s[c] for c in meta_cols), u, name(u), checkv[u]["miuvig_quality"],
                                str(reads_id[k]), f"{covered[k]:.4f}", str(int(covered[k] >= MIN_COVERED)),
                                str(counted[k]), str(s["total"]), f"{r:.3f}", *cols]) + "\n")

with open(f"{OUTD}/votus.tsv", "w") as fo:
    fo.write("votu\tdisplay_name\tlength\tn_contigs\tmiuvig_quality\tcheckv_completeness\tcompleteness_method\t"
             "samples_detected\tmax_rpm\tlca_lineage\n")
    for u in votus:
        c = checkv[u]
        fo.write("\t".join([u, name(u), rep[u]["length"], str(n_contigs[u]), c["miuvig_quality"], c["completeness"],
                            c["completeness_method"], str(len(by_unit[u])),
                            f"{by_unit[u][0][0]:.1f}" if by_unit[u] else "0", rep[u]["lca_lineage"]]) + "\n")

hq = [u for u in votus if checkv[u]["miuvig_quality"] == "High-quality" and by_unit[u]]
with open(f"{OUTD}/viral_rpm.tsv", "w") as fo:
    fo.write("\t".join(["sample", *meta_cols, *(f"{name(u)} [{u}]" for u in hq)]) + "\n")
    for s in samples:
        fo.write("\t".join([s["sample"], *(s[c] for c in meta_cols), *(f"{rpm[(s['sample'], u)]:.3f}" for u in hq)]) + "\n")
print(f"{len(votus)} vOTUs, {sum(checkv[u]['miuvig_quality'] == 'High-quality' for u in votus)} MIUViG high-quality, "
      f"{len(hq)} of them detected in >= 1 sample -> {OUTD}/viral_rpm.tsv")
