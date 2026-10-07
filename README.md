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

# masoudi-pipeline

A single-file Python pipeline for **Oxford Nanopore (ONT)** reference-based analysis.

**Input:** one or many FASTQ(.gz) files (or barcode folders) + a custom reference FASTA
**Output:** filtered VCF, consensus FASTA, QC summary (TSV + MultiQC), QC plots, coverage/depth plots, VCF plots, BAM/BAI/SAM

```
FASTQ ─► [1] FastQC ─► Trim Galore ─► FastQC ─► MultiQC
              │
              ▼
        [2] minimap2 | BWA-MEM2 ─► SAM ─► filter ─► sort ─► BAM+BAI ─► samtools flagstat
              │
              ▼
        [3] samtools markdup ─► flagstat ─► LoFreq  +  bcftools mpileup/call ─► filtered VCFs
              │
              ▼
        [4] bcftools consensus (low-depth positions = N) ─► consensus.fasta
              │
              ▼
        QC summary TSV · MultiQC · QC/coverage/depth/VCF plots
```

---

## 1. Installation

### Linux / macOS (recommended: conda / mamba)

```bash   NO
# 1) get conda (skip if installed) – e.g. Miniforge
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

# 2) create the environment
conda env create -f environment.yml          # or: mamba env create -f environment.yml
conda activate masoudi-pipeline

# 3) check that everything is found
python masoudi_pipeline.py --version
for t in fastqc trim_galore cutadapt multiqc minimap2 bwa-mem2 samtools bcftools bgzip lofreq; do
  command -v $t >/dev/null && echo "OK   $t" || echo "MISSING $t"; done
```

Manual alternative (same result):

```bash    YES
# stop conda from ever using the blocked channel
conda config --remove channels defaults
conda config --set channel_priority strict

# create the environment using ONLY conda-forge and bioconda
conda create -n masoudi-pipeline --override-channels -c conda-forge -c bioconda \
  python=3.11 numpy matplotlib fastqc trim-galore cutadapt multiqc minimap2 bwa-mem2 \
  "samtools>=1.17" "bcftools>=1.17" htslib lofreq

conda activate masoudi-pipeline
```

### Windows
Bioconda tools do not run natively on Windows. Use **WSL2** (Ubuntu): `wsl --install -d Ubuntu`, then follow the Linux steps inside Ubuntu.

### Notes
* **Apple-silicon Mac:** LoFreq and BWA-MEM2 may be unavailable for `osx-arm64`. Create the env with `CONDA_SUBDIR=osx-64 conda env create -f environment.yml`, or run with `--callers bcftools --mapper minimap2`.
* Python packages needed by the script itself: only `numpy` and `matplotlib`.

---

## 2. Quick start

```bash
# one FASTQ
python masoudi_pipeline.py -i /home/alireza/run1/sample/barcode30_merged.fastq.gz -r /home/alireza/reference/reference.fasta -o /home/alireza/ngs-results -t 8

# several FASTQ files  -> each file is its own sample
python masoudi_pipeline.py -i S1.fastq.gz S2.fastq.gz S3.fastq.gz -r reference.fasta -o results -t 8

# glob
python masoudi_pipeline.py -i fastq/*.fastq.gz -r reference.fasta -o results

# ONT barcode folders (all FASTQ inside a folder are concatenated = one sample)
python masoudi_pipeline.py -i fastq_pass/barcode01 fastq_pass/barcode02 -r reference.fasta -o results

# many FASTQ files of ONE sample -> merge them into one
python masoudi_pipeline.py -i run/*.fastq.gz --merge-fastq --sample-name patient1 -r reference.fasta -o results

# BWA-MEM2 instead of minimap2, only bcftools, keep duplicates marked (not removed)
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --mapper bwa-mem2 --callers bcftools --dedup mark

# re-run after a crash / parameter change of later steps: skip finished steps
python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results --resume
```

`python masoudi_pipeline.py -h` prints all options.

---

## 3. All options and the filters they control

### Input / output
| Option | Default | Meaning |
|---|---|---|
| `-i, --input` | required | FASTQ file(s), globs or barcode folder(s) |
| `-r, --reference` | required | custom reference FASTA (≥1 contig/segment) |
| `-o, --outdir` | `masoudi_results` | output folder |
| `-t, --threads` | 8 | threads for every tool |
| `--merge-fastq` / `--sample-name` | off / `merged` | merge all inputs into one sample |
| `--resume` | off | skip steps whose output already exists |

### Step 1 – Trim Galore (`trim_galore --phred33 --gzip ...`)
| Option | Default | Trim Galore flag | Effect |
|---|---|---|---|
| `--trim-quality` | 10 | `--quality` | 3' quality trimming. ONT reads are ~Q10–Q20, so the Illumina default (20) is too strict |
| `--min-length` | 200 | `--length` | discard reads shorter than N bp |
| `--max-length` | off | `--max_length` | discard reads longer than N bp |
| `--stringency` | 3 | `--stringency` | minimum overlap with an adapter to trim |
| `--clip5` / `--clip3` | 0 / 0 | `--clip_R1` / `--three_prime_clip_R1` | fixed bp removed from read start / end |
| `--adapter` | auto | `--adapter` | custom adapter sequence |

FastQC is run on the raw reads and again on the trimmed reads (validation); both are merged by MultiQC (`04_multiqc_trimming/`).

### Step 2 – Mapping and BAM handling
| Option | Default | Effect |
|---|---|---|
| `--mapper` | `minimap2` | `minimap2 -ax map-ont --secondary=no` or `bwa-mem2 mem -x ont2d` |
| `--minimap2-preset` | `map-ont` | e.g. `lr:hq` for very accurate (Q20+) ONT reads |
| `--mapper-extra` | – | extra arguments passed verbatim to the mapper |
| `--min-mapq` | 10 | `samtools view -q` : drop alignments with MAPQ < N |
| `--exclude-flags` | 2308 | `samtools view -F` : removes unmapped (4) + secondary (256) + supplementary (2048) |
| `--no-sam` | off | delete the large SAM at the end (kept by default) |

Order: mapper → `.sam` → `samtools flagstat` (raw alignments) → `samtools view -F -q | samtools sort` → `.sorted.bam` → `samtools index` (`.bai`) → `samtools flagstat` (validation) + `idxstats`.

### Step 3 – Duplicates and SNP calling
| Option | Default | Effect |
|---|---|---|
| `--dedup` | `remove` | `remove` = `samtools markdup -r`; `mark` = flag only; `none` = skip |
| `--callers` | `lofreq bcftools` | any of them |
| `--min-bq` | 7 | min base quality (`lofreq -q/-Q`, `bcftools mpileup -Q`) |
| `--min-depth` | 10 | VCF filter `INFO/DP>=`; LoFreq `-C`; positions below it are masked **N** in the consensus |
| `--min-qual` | 20 | VCF filter `QUAL>=` |
| `--min-af` | 0.10 | VCF filter allele frequency |
| `--indels` | off | also call indels (LoFreq runs `indelqual --dindel` first; bcftools stops skipping indels) |

Commands used:

```bash
# depth / coverage
samtools mpileup -aa -d 0 -Q 0 -q 0 -B -f ref.fasta final.bam | cut -f1,2,4 > sample.depth.tsv
samtools coverage final.bam > sample.coverage.tsv

# LoFreq (+ indelqual when --indels)
lofreq call-parallel --pp-threads T -f ref.fasta -o raw.vcf -q 7 -Q 7 -m 10 -C 10 -B [--call-indels] final.bam

# bcftools (mpileup + call)
bcftools mpileup -f ref.fasta -B -Q 7 -q 10 -d 1000000 -L 1000000 -a FORMAT/AD,FORMAT/DP -Ou final.bam \
  | bcftools call -mv --ploidy 1 [-V indels] -Ov -o raw.vcf

# the same filter logic for both callers (bcftools view -i)
#   LoFreq  : QUAL>=20 && INFO/DP>=10 && INFO/AF>=0.10
#   bcftools: QUAL>=20 && INFO/DP>=10 && (DP4[2]+DP4[3])/sum(DP4)>=0.10
```

> Modern `samtools mpileup` only writes text pileup (no VCF). The pipeline therefore uses `samtools mpileup` for per-base depth and `bcftools mpileup | bcftools call` (its direct successor) for SNP calling.

### Step 4 – Consensus
| Option | Default | Effect |
|---|---|---|
| `--consensus-from` | `lofreq` | which caller's VCF is applied |
| `--consensus-af` | 0.5 | only variants with AF ≥ this are written into the consensus (majority allele) |

```bash
bcftools view -i 'INFO/AF>=0.5' sample.lofreq.filtered.vcf | bcftools norm -f ref.fasta -c w -Oz -o sites.vcf.gz
bcftools index sites.vcf.gz
bcftools consensus -f ref.fasta -m sample.lowdepth_mask.bed -o consensus.fa sites.vcf.gz
```

---

## 4. Output structure

```
results/
├── 00_reference/reference.fasta (+ .fai, BWA-MEM2 index)
├── 00_input/<sample>.fastq.gz
├── 01_qc_raw/          FastQC of raw reads
├── 02_trimmed/         <sample>_trimmed.fq.gz + Trim Galore reports
├── 03_qc_trimmed/      FastQC of trimmed reads (validation)
├── 04_multiqc_trimming/masoudi_multiqc_qc_trimming.html
├── 05_mapping/         <sample>.sam, .sorted.bam/.bai, *.flagstat.txt, idxstats
├── 06_dedup/           <sample>.final.bam/.bai, markdup stats, flagstat, samtools stats
├── 07_variants/        <sample>.final.vcf   <sample>.lofreq|bcftools.raw.vcf / .filtered.vcf / stats
├── 08_consensus/       <sample>.consensus.fasta
├── 09_coverage/        <sample>.depth.tsv, .coverage.tsv, .lowdepth_mask.bed
├── 10_plots/           <sample>.qc.png  <sample>.coverage_depth.png  <sample>.vcf.png
├── 11_summary/         qc_summary.tsv, run_overview.png, masoudi_multiqc_full_report.html
├── logs/               one log per tool call
└── commands.sh         every command that was executed
```

`qc_summary.tsv` columns: reads / bases / mean length / N50 / mean Q (raw and trimmed), % reads retained, % aligned, reads in final BAM, duplicates, mean & median depth, breadth ≥1x and ≥`--min-depth`, variants per caller, consensus length, N count and % N.

Plots: **qc.png** (read length, read quality, read fate), **coverage_depth.png** (depth along genome, depth histogram, cumulative breadth), **vcf.png** (variants along genome by AF, AF spectrum, variant types, substitution spectrum).

---

## 5. Practical advice for ONT data

* **Duplicate removal:** `markdup -r` is on because you asked for it, but for **amplicon / tiled-PCR** data many reads legitimately share start and end positions and removal will wreck your coverage. For amplicons use `--dedup mark` or `--dedup none`.
* **Basecaller trimming:** if Dorado/Guppy already trimmed adapters/barcodes, Trim Galore mostly does quality + length trimming.
* **Indels:** ONT errors concentrate in homopolymers, so indels are off by default. Enable `--indels` only if you will review them manually.
* **Minority variants:** the default `--min-af 0.10` filters out most ONT noise. Lower it only with high depth (>500x) and compare LoFreq vs bcftools.
* **Highly accurate reads (R10.4.1, Q20+):** try `--minimap2-preset lr:hq --trim-quality 15`.
* **Multi-segment references:** all contigs are handled; the consensus header becomes `>sample|contig`.
* **Disk space:** SAM files are large; add `--no-sam` if you do not need them.

## 6. Troubleshooting

| Problem | Fix |
|---|---|
| `Missing tools in PATH` | `conda activate masoudi-pipeline` |
| `0 reads left after trimming` | lower `--trim-quality` / `--min-length` |
| Very low % aligned | wrong reference, or MAPQ too strict: try `--min-mapq 0` |
| Consensus mostly N | depth below `--min-depth`: lower it or sequence deeper |
| LoFreq not installable | `--callers bcftools` |
| A step failed | read `logs/<sample>.<step>.log`, fix, rerun with `--resume` |
