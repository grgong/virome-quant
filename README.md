# Virome quantification (RNA-seq of one host species)

Counts the viruses of a set of RNA-seq libraries in every library. The work is split in two:

- **Discovery** (done separately, on a cluster that holds nr, from 2 M read pairs per library): host
  removal, de novo assembly of the pooled non-host reads (MEGAHIT and rnaviralSPAdes), DIAMOND
  blastx against nr, taxonomy (CAT/BAT 50 % majority rule), clustering of the viral contigs into
  vOTUs (95 % ANI over 85 % of the shorter sequence, MIUViG; Roux et al. 2019) and genome quality
  (CheckV; Nayfach et al. 2021). You receive the result as a **panel**: one sequence per vOTU.
- **This directory** (you, on your cluster): draw the 2 M read pairs per library and send them
  (`run.sh subsample`); then map all reads of every library to the panel and count them
  (`run.sh quantify`). No database is needed here.

## Setup (once)

```bash
. /local/env/envconda.sh                                   # GenOuest; elsewhere: your conda
conda env create -p /scratch/$USER/envs/virome_quant -f environment.yml
cp config.sh my_config.sh                                  # then set CONDA_ENV (and the partition) in it
```

## 1. Samplesheet

Tab-separated, one row per library: `sample`, `fastq_1`, `fastq_2`, `flowcell_lane`, plus any
columns of your own (`stage`, `tissue`, ...), which are copied to the result tables. Sample
names: letters, digits, `.`, `_`, `-`. A library sequenced on several lanes: comma-separated files
in `fastq_1` and `fastq_2` (same order) and `FLOWCELL:LANE` items in `flowcell_lane`. The flowcell
and lane are fields 3 and 4 of a read name: `zcat x_R1.fastq.gz | head -1 | cut -d: -f3,4`.
Example: `samples.example.tsv`. The same samplesheet is used in steps 2 and 4.

## 2. Subsample and send

```bash
bash run.sh subsample my_config.sh samples.tsv /scratch/$USER/virome_2M
```

When the jobs have finished (`squeue -u $USER`), send the whole output directory: the reads
(`reads/`, about 280 MB per library), `samples_2M.tsv` and the md5 sums (`reads/*.md5`). The
subsample is a random draw (seqtk, seed 11) of 2,000,000 read pairs per library, the same draw
the discovery would make from the full files.

## 3. Receive the panel

A directory with three files: `panel.fa` (the vOTU sequences), `virus_units.tsv` (the viral
contigs of every vOTU, its name and lineage) and `quality_summary.tsv` (CheckV).

## 4. Quantify

```bash
bash run.sh quantify my_config.sh samples.tsv <panel_dir> /scratch/$USER/virome_quant
```

Every library: all reads -> fastp (adapters, poly-G, low complexity, >= 50 bp) -> minimap2 to the
panel. Counting (CoverM): reads at >= 90 % identity with >= 80 % of the read aligned; a vOTU is
detected in a library when these reads cover >= 75 % of it (Roux et al. 2017), otherwise it counts
0. RPM = reads per million reads passing fastp. Results in `<outdir>/results/` (below).

Time: about 10 min per 100 M reads for mapping, at most `MAX_PARALLEL` libraries at a time.
Space: each running mapping job writes up to ~20 GB to `TMPDIR_JOBS` (removed at its end); the
BAMs in `<outdir>/map/` hold only reads on the panel, ~30 bytes per read (a library with 165 M
viral reads: 4.9 GB). The BAMs are needed only to count again; delete them once the results are
written if space is short.

**If a job fails**: read its log in `<outdir>/logs/`, fix the cause, then
`bash run.sh resume <outdir>` (finished libraries are not redone).

## First run: libraries of known status

Before the full set, run step 4 on a few libraries whose viruses are known (for instance one
without, one with a moderate and one with a high titre of a virus the discovery found) with the
panel made from them, and check that the counts are those the discovery side gives. The test
set of the original project (`test/`: samplesheet, panel and the expected counts) is supplied
with the package by the discovery side and is not part of the public repository.

## Results (`<outdir>/results/`)

| File | Content |
|---|---|
| `viral_rpm.tsv`, `viral_rpm.pdf` | libraries x high-quality vOTUs detected in at least one library: RPM (0 = not detected) |
| `votus.tsv` | every vOTU of the panel: name (closest relative, "-like" below 90 % aa identity), length, CheckV tier and completeness, libraries detected, highest RPM, lineage |
| `viral_long.tsv` | library x vOTU, every vOTU: reads at >= 90 % identity, covered fraction, detected, RPM, and the carry-over columns below |

## Reading the results

- **Only high-quality vOTUs are in `viral_rpm.tsv`** (CheckV MIUViG high-quality, > 90 %
  complete). Viruses assembled only in part are in `votus.tsv` and `viral_long.tsv`: a virus
  whose longest contig is 3.4 kb of a ~10 kb genome (CheckV 70 %) is not in `viral_rpm.tsv`
  although it is present. CheckV's database under-represents RNA viruses, so a genuine but
  divergent virus can be "not determined".
- **Detection floor.** Covering 75 % of a 10 kb genome takes about 90 reads, ~1 RPM in a
  library of 100 M reads. Below that a vOTU counts 0 even with reads on it; `reads_id90` and
  `covered_fraction` in `viral_long.tsv` show them (10 reads covering 11 % of a 10 kb genome count 0).
- **Carry-over between libraries.** Reads of a virus can pass from one library to another on the
  same flowcell at up to ~0.1 % (Pettersson 2019, French 2022, Waller 2024; 0.3 % within a lane,
  Asplund 2019). `pct_of_max_same_flowcell` gives a library's RPM as a percentage of the highest
  RPM of that vOTU in another library on its flowcell; below 0.1 % is the level of carry-over. A
  library with 93 % of its reads from one virus (~930,000 RPM) puts up to ~930 RPM into the
  libraries on its flowcell, the RPM of a moderately infected library: on a shared flowcell such
  a library cannot be told from carry-over by RPM alone.
- **Composition.** RPM is per million of all reads. In a library where most reads are one virus,
  every other virus has a lower RPM for the same number of copies per host cell.
- New libraries later: map them to the same panel (step 4). A virus that is not in the panel is
  not counted; to look for new viruses, send their 2 M subsamples (step 2).
- The numbers are descriptive; nothing here makes a statistical test.

Origin: a viral screen of an aphid host-shift RNA-seq study (v5, 2026-10/11); same tools,
versions and thresholds.
