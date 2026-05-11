# WGS-CNA

A Nextflow DSL2 pipeline for calling Copy Number Alterations (CNA) from whole-genome sequencing (WGS) and low-pass WGS (lp-WGS) data.

The pipeline preprocesses raw FASTQ files through trimming, alignment, and duplicate removal (via the [fastq-preprocess](https://github.com/nf-austin/fastq-preprocess) submodule), then calls CNAs using [HMMcopy readCounter](https://bioconductor.org/packages/HMMcopy/) and [ichorCNA](https://github.com/broadinstitute/ichorCNA).

---

## Pipeline steps

1. **Adapter trimming & QC** (`fastp`) — trims adapters, filters low-quality reads, and produces per-sample JSON/HTML reports.
2. **Reference genome download** (`wget`) — downloads and decompresses hg38 or hg19 from UCSC. Skipped if `--bwa_index` is provided.
3. **BWA-MEM2 indexing** — builds the BWA-MEM2 index from the downloaded FASTA. Skipped if `--bwa_index` is provided.
4. **Alignment** (`bwa-mem2 mem`) — aligns trimmed reads to the reference with read group tags. Each lane is aligned independently.
5. **Coordinate sort + fixmate** (`samtools`) — prepares the BAM for duplicate marking; steps are piped to avoid writing intermediate files.
6. **Duplicate removal** (`samtools markdup`) — removes PCR/optical duplicates.
7. **Lane merging** (`samtools merge`) — merges per-lane BAMs into a single BAM per sample. Single-lane samples skip this step.
8. **Indexing** (`samtools index`) — indexes the final coordinate-sorted BAM.
9. **Alignment QC** (`samtools stats/flagstat/idxstats`) — produces per-sample alignment statistics.
10. **Read counting** (`hmmcopy readCounter`) — bins reads into fixed-size genomic windows to produce a WIG file per sample.
11. **CNA calling** (`ichorCNA`) — estimates copy number states, tumor fraction, and ploidy from the read-count WIG files using a hidden Markov model.
12. **MultiQC report** — aggregates fastp and samtools QC outputs into a single HTML report.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04
- Docker **or** Singularity (all tools run in containers — nothing else needs to be installed)

---

## Running from GitHub (no clone needed)

Nextflow can pull and run this pipeline straight from GitHub — only Nextflow and Docker/Singularity need to be installed locally. The `fastq-preprocess` submodule is fetched automatically thanks to `manifest.recurseSubmodules = true` in `nextflow.config`.

```bash
nextflow run nf-austin/wgs-cna \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz"
```

Pin to a specific release tag for reproducibility:

```bash
nextflow run nf-austin/wgs-cna \
    -r v1.0.0 \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz"
```

Nextflow caches the pipeline under `~/.nextflow/assets/`. To pull a fresh copy after an update:

```bash
nextflow pull nf-austin/wgs-cna
```

> **Note**: `recurseSubmodules` requires Nextflow ≥ 22.10. If you are stuck on an older version, clone the repo manually as described below.

---

## Installation (local clone)

If you prefer a local checkout — e.g. to modify the code — clone the repository and initialise the `fastq-preprocess` submodule:

```bash
git clone https://github.com/<your-org>/wgs-cna.git
cd wgs-cna
git submodule update --init --recursive
```

---

## Usage

### Minimal run (auto-download hg38 reference, standard WGS)

```bash
nextflow run main.nf \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz"
```

### Supply a pre-built BWA-MEM2 index (skips download + indexing)

```bash
nextflow run main.nf \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --bwa_index "/path/to/bwa_index/"
```

### hg19 reference

```bash
nextflow run main.nf \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --genome hg19
```

### Low-pass WGS mode

Sets a 1 Mb window, applies the ichorCNA ultra-low-pass panel of normals, and broadens the normal fraction and ploidy search space:

```bash
nextflow run main.nf \
    -profile docker \
    --run_lp_wgs true \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --genome hg38
```

### Singularity

```bash
nextflow run main.nf \
    -profile singularity \
    --fastqs "data/*_{1,2}.fastq.gz"
```

### Custom output directory

```bash
nextflow run main.nf \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --outdir /path/to/results
```

### Provide pre-downloaded ichorCNA resource files

If you already have local copies of the WIG / panel-of-normals / centromere files:

```bash
nextflow run main.nf \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --bwa_index "/path/to/bwa_index/" \
    --gc_wig   "/path/to/gc_hg38_500kb.wig" \
    --map_wig  "/path/to/map_hg38_500kb.wig" \
    --pon_rds  "/path/to/HD_ULP_PoN_hg38_500kb.rds" \
    --centromere_txt "/path/to/GRCh38_centromere.txt"
```

---

## Parameters

### Core inputs

| Parameter | Default | Description |
|---|---|---|
| `--fastqs` | *required* | Glob pattern for paired FASTQ files (e.g. `"data/*_{1,2}.fastq.gz"`) |
| `--outdir` | `results` | Output directory |
| `--genome` | `hg38` | Reference genome (`hg38`, `hg19`). Ignored if `--bwa_index` is set. |
| `--bwa_index` | `null` | Path to a pre-built BWA-MEM2 index directory. Skips download and indexing. |

### Pipeline mode

| Parameter | Default | Description |
|---|---|---|
| `--run_lp_wgs` | `false` | Enable low-pass WGS mode (adjusts window size, ploidy priors, and normal fraction search) |

### CNA / ichorCNA tuning

| Parameter | Default (standard WGS) | Default (lp-WGS) | Description |
|---|---|---|---|
| `--cna_window_size` | `500000` | `1000000` | Read-count bin size in bp |
| `--cna_min_mapq` | `20` | `20` | Minimum mapping quality for read counting |
| `--ichorcna_ploidy` | `c(2)` | `c(2,3)` | Ploidy states to evaluate |
| `--ichorcna_normal` | `c(0.5)` | `c(0.5,0.6,0.7,0.8,0.9,0.95)` | Normal cell fraction states to evaluate |
| `--ichorcna_max_cn` | `5` | `5` | Maximum copy number state |
| `--ichorcna_txn_e` | `0.99` | `0.9999` | HMM transition probability (higher = fewer breakpoints) |
| `--ichorcna_txn_strength` | `10000` | `10000` | Transition prior strength |

### ichorCNA resource files (auto-downloaded if not provided)

| Parameter | Description |
|---|---|
| `--gc_wig` | GC content WIG file for the chosen genome/window combination |
| `--map_wig` | Mappability WIG file |
| `--pon_rds` | Panel of normals RDS file (lp-WGS only; ignored in standard WGS unless provided manually) |
| `--centromere_txt` | Centromere positions text file |

Resource files for `hg19` and `hg38` at both 500 kb and 1 Mb window sizes are auto-downloaded from the [ichorCNA GitHub repository](https://github.com/broadinstitute/ichorCNA) when not supplied. CNA calling currently requires `hg19` or `hg38`; the alignment step additionally supports `mm10` and `mm39` (but CNA will fail without manually supplied WIG files for those genomes).

---

## Output structure

```
results/
├── qc/
│   ├── fastp/              # per-sample .json and .html reports
│   ├── samtools/           # per-sample stats, flagstat, idxstats
│   └── multiqc/            # multiqc_report.html + multiqc_data/
├── aligned/                # final .bam and .bai per sample
├── reference/              # downloaded FASTA and BWA-MEM2 index (if built)
└── cna/
    ├── readcounts/         # per-sample .wig read-count files
    └── ichorcna/           # per-sample ichorCNA results:
        │                   #   {id}.cna.seg   — IGV-compatible CNA segments
        │                   #   {id}.seg.txt   — tabular segment summary
        │                   #   {id}.params.txt — estimated tumor fraction & ploidy
        │                   #   {id}.RData      — full R data object
        └──                 #   {id}_genomeWide.pdf — genome-wide CNA plot
```

---

## Multi-lane FASTQ naming

The pipeline expects Illumina-style lane-tagged filenames:

```
SAMPLE_S1_L001_R1_001.fastq.gz
SAMPLE_S1_L001_R2_001.fastq.gz
SAMPLE_S1_L002_R1_001.fastq.gz
SAMPLE_S1_L002_R2_001.fastq.gz
```

Each lane is trimmed and aligned independently, then per-lane BAMs are merged into a single `SAMPLE_S1.bam` before duplicate marking and QC. The `_L001`-style suffix is stripped to determine the sample name.