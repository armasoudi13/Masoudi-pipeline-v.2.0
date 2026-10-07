```
            .-~~~~~~~-.
         .'   _        '.
        /    (o)    _.-~'
       |          .-'              .    .    .    .    .    .    .
       |          '-.
        \           '-._
         '.             .'
           '-.._____..-'
 __  __     _     ____    ___   _   _  ____   ___
|  \/  |   / \   / ___|  / _ \ | | | ||  _ \ |_ _|
| |\/| |  / _ \  \___ \ | | | || | | || | | | | |
| |  | | / ___ \  ___) || |_| || |_| || |_| | | |
|_|  |_|/_/   \_\|____/  \___/  \___/ |____/ |___|
          ~ P I P E L I N E ~   v1.0.0
```

# masoudi-pipeline — Complete Command-Line Guide

This guide explains **every command and option** of `masoudi_pipeline.py`, with examples you can copy and adapt.

Contents

1. Before you start
2. Basic syntax
3. Input options (`-i`, `-r`, `-o`, `-t`, merging, resume)
4. Step 1 options: QC and trimming
5. Step 2 options: mapping and BAM filtering
6. Step 3 options: duplicates and SNP calling
7. Step 4 options: consensus
8. Ready-to-use scenarios
9. Batch processing, servers, and automation
10. Reading the results
11. Running single steps manually
12. Troubleshooting with commands
13. Option quick-reference table

---

## 1. Before you start

```bash
# activate the environment (every new terminal)
conda activate masoudi-pipeline

# go to your working folder
cd /path/to/my_project

# check that the pipeline and tools are found
python masoudi_pipeline.py --version
python masoudi_pipeline.py -h
```

Prepare:

| Item | Requirement |
|---|---|
| FASTQ | `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz` (basecalled ONT reads) |
| Reference | one FASTA file; can contain several contigs/segments |
| Disk space | roughly 5–10× the size of your FASTQ (SAM + BAM + trimmed reads) |

---

## 2. Basic syntax

```bash
python masoudi_pipeline.py -i <FASTQ files or folders> -r <reference.fasta> -o <output folder> [options]
```

Minimum command (all defaults):

```bash
python masoudi_pipeline.py -i reads.fastq.gz -r reference.fasta -o results
```

Tip for long commands: use `\` at the end of a line to continue on the next line.

```bash
python masoudi_pipeline.py \
  -i reads.fastq.gz \
  -r reference.fasta \
  -o results \
  -t 16
```

---

## 3. Input and general options

### `-i / --input`  (required)

One or more FASTQ files, wildcards, or folders. **Each file or folder becomes one sample**, and the sample name is the file name without `.fastq.gz`.

```bash
# one file  -> sample "reads"
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results

# several files -> samples S1, S2, S3
python masoudi_pipeline.py -i S1.fastq.gz S2.fastq.gz S3.fastq.gz -r ref.fasta -o results

# wildcard (all .fastq.gz in a folder)
python masoudi_pipeline.py -i fastq/*.fastq.gz -r ref.fasta -o results

# barcode folders: all FASTQ inside each folder are concatenated -> samples barcode01, barcode02
python masoudi_pipeline.py -i fastq_pass/barcode01 fastq_pass/barcode02 -r ref.fasta -o results

# all barcode folders at once
python masoudi_pipeline.py -i fastq_pass/barcode* -r ref.fasta -o results
```

Two files or folders with the same name are rejected (to avoid overwriting). Rename one, or use `--merge-fastq`.

### `-r / --reference`  (required)

Your custom reference FASTA. The pipeline copies it to `00_reference/reference.fasta` and creates the indexes there, so your original folder stays clean.

```bash
python masoudi_pipeline.py -i reads.fastq.gz -r /data/refs/Sabin2_VP1.fasta -o results
```

### `-o / --outdir`

Output folder (created if missing). Default: `masoudi_results`.

### `-t / --threads`

Threads for all tools. Default 8. Use the number of CPU cores you have (check with `nproc`).

```bash
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results -t 24
```

### `--merge-fastq` and `--sample-name`

Put **all inputs into one sample** (for example a run split into many FASTQ chunks).

```bash
python masoudi_pipeline.py \
  -i run/*.fastq.gz \
  --merge-fastq --sample-name patient_07 \
  -r ref.fasta -o results_patient07
```

Without `--sample-name` the merged sample is called `merged`.

### `--resume`

Skip every step whose output file already exists. Use it after a crash or when you only changed later-step options.

```bash
# first run (crashed or stopped)
python masoudi_pipeline.py -i fastq/*.fastq.gz -r ref.fasta -o results -t 16

# run again with the SAME -o folder
python masoudi_pipeline.py -i fastq/*.fastq.gz -r ref.fasta -o results -t 16 --resume
```

Important: `--resume` does not know that you changed a parameter. If you change a filter (for example `--min-mapq`) and want that step recomputed, delete the matching output files first, or use a new `-o` folder.

```bash
# example: you changed --min-af, so redo VCF filtering and consensus only
rm results/07_variants/*.filtered.vcf results/08_consensus/*
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --min-af 0.05 --resume
```

---

## 4. Step 1 options: FastQC → Trim Galore → FastQC → MultiQC

What happens: FastQC on raw reads, then Trim Galore, then FastQC on trimmed reads (validation), then MultiQC combines the reports.

| Option | Default | Trim Galore flag | Meaning |
|---|---|---|---|
| `--trim-quality N` | 10 | `--quality` | trim low-quality bases from the 3' end |
| `--min-length N` | 200 | `--length` | discard reads shorter than N bp |
| `--max-length N` | off | `--max_length` | discard reads longer than N bp |
| `--stringency N` | 3 | `--stringency` | minimum adapter overlap to trim |
| `--clip5 N` | 0 | `--clip_R1` | remove N bp from the start of every read |
| `--clip3 N` | 0 | `--three_prime_clip_R1` | remove N bp from the end of every read |
| `--adapter SEQ` | auto | `--adapter` | use this adapter sequence |

Examples:

```bash
# stricter quality, longer minimum length
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --trim-quality 15 --min-length 500

# keep only reads of a defined size range (e.g. a 1.5 kb amplicon)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --min-length 1200 --max-length 1800

# remove 25 bp of primer sequence from both ends
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --clip5 25 --clip3 25

# custom adapter
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --adapter AGATCGGAAGAGC
```

How to choose:

| Situation | Suggested values |
|---|---|
| Standard R9.4.1 reads (Q10–Q15) | `--trim-quality 10 --min-length 200` (default) |
| High-accuracy reads (Q20+) | `--trim-quality 15 --min-length 300` |
| Amplicons of known size | set `--min-length` and `--max-length` around the amplicon size |
| Adapters already removed by Dorado/Guppy | defaults are fine |

Look at: `04_multiqc_trimming/masoudi_multiqc_qc_trimming.html` (compare raw vs trimmed).

---

## 5. Step 2 options: mapping and BAM filtering

What happens: mapper → SAM → `samtools flagstat` → filter + sort → BAM + BAI → `samtools flagstat` (validation).

| Option | Default | Meaning |
|---|---|---|
| `--mapper` | `minimap2` | `minimap2` (`-ax map-ont --secondary=no`) or `bwa-mem2` (`mem -x ont2d`) |
| `--minimap2-preset` | `map-ont` | minimap2 `-x` preset: `map-ont`, `lr:hq` (Q20+ reads), `splice` (not for DNA viruses) |
| `--mapper-extra "..."` | none | extra arguments passed to the mapper unchanged |
| `--min-mapq N` | 10 | keep alignments with mapping quality ≥ N (`samtools view -q`) |
| `--exclude-flags N` | 2308 | drop alignments with these SAM flags (`samtools view -F`) |
| `--no-sam` | off | delete the SAM file(s) at the end to save space |

About `--exclude-flags`: the default 2308 = 4 (unmapped) + 256 (secondary) + 2048 (supplementary). To keep supplementary alignments use 260 (4 + 256).

Examples:

```bash
# BWA-MEM2 mapping
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --mapper bwa-mem2

# accurate reads preset
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --minimap2-preset lr:hq

# more sensitive mapping for a divergent reference (smaller k-mer)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --mapper-extra="-k 11 -w 5"

# keep every mapped read regardless of MAPQ
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --min-mapq 0

# keep supplementary alignments
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --exclude-flags 260

# save disk space
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --no-sam
```

Check: `05_mapping/<sample>.sorted.flagstat.txt` (mapped reads) and the `pct_aligned` column of `11_summary/qc_summary.tsv`.

---

## 6. Step 3 options: duplicates and SNP calling

### Duplicates: `--dedup`

| Value | What it does |
|---|---|
| `remove` (default) | `samtools markdup -r`: duplicates are deleted from the final BAM |
| `mark` | duplicates are flagged but kept (most tools ignore flagged reads) |
| `none` | no duplicate handling |

```bash
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --dedup mark
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --dedup none
```

Warning: for amplicon or PCR-tiled data, reads with identical start/end are expected, and `remove` will cut coverage heavily. Use `mark` or `none` there.

### Callers: `--callers`

```bash
# both (default)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --callers lofreq bcftools

# only LoFreq
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --callers lofreq

# only bcftools (also use it for the consensus)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --callers bcftools --consensus-from bcftools
```

### Variant filters

| Option | Default | Where it is applied |
|---|---|---|
| `--min-bq N` | 7 | minimum base quality (LoFreq `-q/-Q`, bcftools mpileup `-Q`) |
| `--min-depth N` | 10 | VCF filter `INFO/DP>=N`; LoFreq `-C`; positions below N become **N** in the consensus |
| `--min-qual N` | 20 | VCF filter `QUAL>=N` |
| `--min-af X` | 0.10 | VCF filter allele frequency ≥ X |
| `--indels` | off | also call indels |

The same three conditions (QUAL, depth, allele frequency) are applied to both callers when the `*.filtered.vcf` files are made.

Examples:

```bash
# stricter calls for high-depth data
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --min-depth 50 --min-qual 40 --min-af 0.20

# sensitive minority-variant search (high depth only!)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --min-depth 200 --min-af 0.02 --min-bq 10

# include indels
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --indels
```

Guidance on `--min-af` for ONT:

| Goal | Value |
|---|---|
| Reliable consensus-level changes | 0.30–0.50 |
| Default screening | 0.10 |
| Minority variants (need >500x depth and confirmation) | 0.02–0.05 |

ONT systematic errors can reach several percent, so very low thresholds produce false positives.

---

## 7. Step 4 options: consensus

| Option | Default | Meaning |
|---|---|---|
| `--consensus-from` | `lofreq` | which caller's VCF is applied to the reference |
| `--consensus-af X` | 0.5 | only variants with AF ≥ X enter the consensus |

Rules:

- Variants must first pass the filters of Step 3 and then also AF ≥ `--consensus-af`.
- Positions with depth below `--min-depth` are written as **N**.
- If `--consensus-from` names a caller that was not run, the first available caller is used.

```bash
# majority-rule consensus (default)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results

# very strict: variant must be present in ≥70% of reads, depth ≥ 30 or N
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --consensus-af 0.7 --min-depth 30

# consensus from bcftools calls
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results \
  --callers lofreq bcftools --consensus-from bcftools
```

Output header: single-contig reference → `>sample`; multi-contig → `>sample|contig`.

---

## 8. Ready-to-use scenarios

### A. Standard single sample

```bash
python masoudi_pipeline.py \
  -i sample01.fastq.gz \
  -r reference.fasta \
  -o results_sample01 \
  -t 16
```

### B. Poliovirus/enterovirus run with barcodes (several samples, same reference)

```bash
python masoudi_pipeline.py \
  -i fastq_pass/barcode* \
  -r Sabin_reference.fasta \
  -o results_run01 -t 16 \
  --trim-quality 12 --min-length 300 \
  --min-depth 20 --min-af 0.15
```

### C. Amplicon sequencing (known amplicon size, no duplicate removal)

```bash
python masoudi_pipeline.py \
  -i amplicon.fastq.gz -r ref.fasta -o results_amplicon -t 16 \
  --min-length 1200 --max-length 1800 \
  --dedup none \
  --min-depth 50 --consensus-af 0.7
```

### D. High-accuracy (R10.4.1) reads

```bash
python masoudi_pipeline.py \
  -i reads.fastq.gz -r ref.fasta -o results_hq -t 16 \
  --minimap2-preset lr:hq --trim-quality 15 --min-bq 10
```

### E. Low-depth sample (relax filters)

```bash
python masoudi_pipeline.py \
  -i lowdepth.fastq.gz -r ref.fasta -o results_low -t 8 \
  --min-length 150 --min-mapq 0 --min-depth 5 --min-af 0.25
```

### F. One patient, many FASTQ chunks

```bash
python masoudi_pipeline.py \
  -i run/*.fastq.gz --merge-fastq --sample-name patient_07 \
  -r ref.fasta -o results_patient07
```

### G. Fast run, minimal disk use

```bash
python masoudi_pipeline.py \
  -i reads.fastq.gz -r ref.fasta -o results -t 16 \
  --callers bcftools --consensus-from bcftools --no-sam
```

### H. Compare two mappers on the same data

```bash
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o res_minimap2 --mapper minimap2
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o res_bwamem2  --mapper bwa-mem2

# compare the summaries
paste <(cut -f1,14,17,18 res_minimap2/11_summary/qc_summary.tsv) \
      <(cut -f14,17,18   res_bwamem2/11_summary/qc_summary.tsv) | column -t
```

---

## 9. Batch processing, servers, and automation

### Run in the background (the job continues if you close the terminal)

```bash
nohup python masoudi_pipeline.py -i fastq/*.fastq.gz -r ref.fasta \
  -o results -t 32 > run.log 2>&1 &

tail -f run.log        # watch progress (Ctrl+C stops watching, not the job)
```

### Each sample has its own reference

```bash
# samples.tsv  (sample_file <TAB> reference)
#   s1.fastq.gz   refA.fasta
#   s2.fastq.gz   refB.fasta
while IFS=$'\t' read -r fq ref; do
  name=$(basename "$fq" .fastq.gz)
  python masoudi_pipeline.py -i "$fq" -r "$ref" -o "results/$name" -t 8
done < samples.tsv
```

### Slurm cluster job

```bash
#!/bin/bash
#SBATCH --job-name=masoudi
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

source ~/miniforge3/etc/profile.d/conda.sh
conda activate masoudi-pipeline

python masoudi_pipeline.py -i fastq/*.fastq.gz -r ref.fasta -o results -t $SLURM_CPUS_PER_TASK
```

Submit with `sbatch run_masoudi.sh`.

### Keep a record of every run

```bash
cat results/commands.sh          # all commands that were executed
cat results/masoudi_pipeline.log # progress log including your command line
ls results/logs/                 # one log per tool and sample
```

---

## 10. Reading the results

```
results/
├── 01_qc_raw/ ... 03_qc_trimmed/      FastQC reports
├── 04_multiqc_trimming/               MultiQC for step 1
├── 05_mapping/                        SAM, sorted BAM + BAI, flagstat
├── 06_dedup/                          final BAM + BAI, markdup stats, flagstat
├── 07_variants/                       VCF files
├── 08_consensus/                      consensus FASTA
├── 09_coverage/                       depth table, coverage table, N-mask BED
├── 10_plots/                          PNG plots
├── 11_summary/                        qc_summary.tsv, overview plot, full MultiQC
├── logs/  commands.sh  masoudi_pipeline.log
```

Useful commands:

```bash
# QC summary as a readable table
column -t -s $'\t' results/11_summary/qc_summary.tsv | less -S

# number of variants in the final VCF
grep -vc '^#' results/07_variants/sample01.final.vcf

# first lines of the final VCF
grep -v '^##' results/07_variants/sample01.final.vcf | head

# consensus length and number of N
python3 - <<'EOF'
s=''.join(l.strip() for l in open('results/08_consensus/sample01.consensus.fasta') if not l.startswith('>'))
print(len(s), 'bp,', s.upper().count('N'), 'N')
EOF

# mapped reads per contig
cat results/05_mapping/sample01.idxstats.txt

# open reports (Linux)
xdg-open results/11_summary/masoudi_multiqc_full_report.html
xdg-open results/10_plots/sample01.coverage_depth.png
```

Main files per sample:

| File | Content |
|---|---|
| `07_variants/<s>.final.vcf` | variants used for the consensus |
| `07_variants/<s>.lofreq.filtered.vcf` | LoFreq calls after filtering |
| `07_variants/<s>.bcftools.filtered.vcf` | bcftools calls after filtering |
| `07_variants/<s>.*.raw.vcf` | unfiltered calls |
| `08_consensus/<s>.consensus.fasta` | consensus sequence |
| `06_dedup/<s>.final.bam` + `.bai` | final alignment |
| `09_coverage/<s>.depth.tsv` | depth at every position |
| `09_coverage/<s>.lowdepth_mask.bed` | regions masked with N |
| `10_plots/<s>.qc.png` | read length, quality, read fate |
| `10_plots/<s>.coverage_depth.png` | depth track, histogram, cumulative coverage |
| `10_plots/<s>.vcf.png` | variant positions, AF spectrum, types, substitutions |

`qc_summary.tsv` columns: `reads_raw`, `bases_raw`, `mean_len_raw`, `n50_raw`, `meanQ_raw`, the same four for trimmed reads, `pct_reads_retained`, `reads_aligned_filtered`, `pct_aligned`, `reads_final_bam`, `duplicates_flagged`, `mean_depth`, `median_depth`, `breadth_1x`, `breadth_min_depth`, `lofreq_variants`, `bcftools_variants`, `consensus_len`, `consensus_N`, `consensus_N_pct`.

Rough quality checks:

| Metric | Good sign |
|---|---|
| `pct_reads_retained` | > 70% |
| `pct_aligned` | > 80% for the correct reference |
| `breadth_min_depth` | close to 100% |
| `consensus_N_pct` | < 5% |

---

## 11. Running single steps manually

Every command the pipeline ran is saved in `results/commands.sh`. Examples of the main steps (replace names with yours):

```bash
# FastQC
fastqc -t 8 -o qc_out reads.fastq.gz

# Trim Galore
trim_galore --phred33 --gzip --quality 10 --length 200 --basename sample -o trimmed reads.fastq.gz

# minimap2 mapping to SAM
minimap2 -ax map-ont -t 8 --secondary=no ref.fasta sample_trimmed.fq.gz > sample.sam

# filter, sort, index, flagstat
samtools view -@ 8 -u -F 2308 -q 10 sample.sam | samtools sort -@ 8 -o sample.sorted.bam -
samtools index sample.sorted.bam
samtools flagstat sample.sorted.bam

# remove duplicates
samtools sort -n sample.sorted.bam | samtools fixmate -m - - | samtools sort - | samtools markdup -r -f markdup.txt - sample.final.bam
samtools index sample.final.bam

# depth with samtools mpileup
samtools mpileup -aa -d 0 -Q 0 -q 0 -B -f ref.fasta sample.final.bam | cut -f1,2,4 > depth.tsv

# LoFreq SNP calling
lofreq call-parallel --pp-threads 8 -f ref.fasta -o sample.lofreq.vcf -q 7 -Q 7 -m 10 -C 10 -B sample.final.bam

# bcftools SNP calling
bcftools mpileup -f ref.fasta -B -Q 7 -q 10 -d 1000000 -a FORMAT/AD,FORMAT/DP -Ou sample.final.bam \
  | bcftools call -mv --ploidy 1 -V indels -Ov -o sample.bcftools.vcf

# filter a VCF
bcftools view -i 'QUAL>=20 && INFO/DP>=10 && INFO/AF>=0.10' sample.lofreq.vcf -Ov -o sample.lofreq.filtered.vcf

# consensus
bcftools view -i 'INFO/AF>=0.5' sample.lofreq.filtered.vcf | bcftools norm -f ref.fasta -c w -Oz -o sites.vcf.gz
bcftools index sites.vcf.gz
bcftools consensus -f ref.fasta -m lowdepth_mask.bed -o consensus.fa sites.vcf.gz

# MultiQC
multiqc -o multiqc_out qc_out trimmed
```

---

## 12. Troubleshooting with commands

| Problem | Check / fix |
|---|---|
| `Missing tools in PATH` | `conda activate masoudi-pipeline` then `which samtools lofreq` |
| A sample failed | `ls results/logs/` then `tail -30 results/logs/<sample>.<step>.log` |
| `0 reads left after trimming` | rerun with `--trim-quality 7 --min-length 100` |
| Very low `pct_aligned` | wrong reference? try `--min-mapq 0`; BLAST a few reads to identify the virus |
| Consensus is mostly N | check `breadth_min_depth`; rerun with lower `--min-depth` (e.g. 5) or sequence deeper |
| LoFreq cannot be installed | run with `--callers bcftools --consensus-from bcftools` |
| Run stopped in the middle | rerun the same command with `--resume` |
| Out of disk space | add `--no-sam`; delete `00_input` links and old SAM files |
| Out of memory with BWA-MEM2 | use the default `--mapper minimap2` |
| Changed a parameter but results are identical | `--resume` skipped the step; delete that step's output files or use a new `-o` |

Check free disk space and size of the results folder:

```bash
df -h .
du -sh results/*
```

---

## 13. Option quick-reference

| Option | Default | Group |
|---|---|---|
| `-i, --input` | required | input |
| `-r, --reference` | required | input |
| `-o, --outdir` | `masoudi_results` | output |
| `-t, --threads` | 8 | general |
| `--merge-fastq` | off | input |
| `--sample-name` | `merged` | input |
| `--resume` | off | general |
| `--trim-quality` | 10 | step 1 |
| `--min-length` | 200 | step 1 |
| `--max-length` | off | step 1 |
| `--stringency` | 3 | step 1 |
| `--clip5` / `--clip3` | 0 / 0 | step 1 |
| `--adapter` | auto | step 1 |
| `--mapper` | `minimap2` | step 2 |
| `--minimap2-preset` | `map-ont` | step 2 |
| `--mapper-extra` | none | step 2 |
| `--min-mapq` | 10 | step 2 |
| `--exclude-flags` | 2308 | step 2 |
| `--no-sam` | off | step 2 |
| `--dedup` | `remove` | step 3 |
| `--callers` | `lofreq bcftools` | step 3 |
| `--min-bq` | 7 | step 3 |
| `--min-depth` | 10 | step 3 / 4 |
| `--min-qual` | 20 | step 3 |
| `--min-af` | 0.10 | step 3 |
| `--indels` | off | step 3 |
| `--consensus-from` | `lofreq` | step 4 |
| `--consensus-af` | 0.5 | step 4 |
| `--version` | – | info |
| `-h, --help` | – | info |
