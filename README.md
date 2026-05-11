# WGS-CNA Workflow

A Nextflow pipeline that supports calling Copy Number Alterations (CNA) and Chromosomal Instability (CIN) from Whole Genome Sequencing (WGS) and low-pass WGS (lp-WGS) data.

## Usage

### Standard WGS

```
nextflow run main.nf \
    -profile docker \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --genome hg38
```

### Low-pass WGS

```
nextflow run main.nf \
    -profile docker \
    --run_lp_wgs true \
    --fastqs "data/*_{1,2}.fastq.gz" \
    --genome hg38

```