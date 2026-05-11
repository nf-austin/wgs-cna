#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_GENOME } from './modules/fastq-preprocess/subworkflows/prepare_genome'
include { QC             } from './modules/fastq-preprocess/subworkflows/qc'
include { ALIGN          } from './modules/fastq-preprocess/subworkflows/align'
include { MULTIQC        } from './modules/fastq-preprocess/modules/multiqc/main'
include { CNA            } from './subworkflows/cna'

workflow {
    if (!params.fastqs) { error "Please provide input FASTQ files using --fastqs" }

    // lp-WGS overrides the standard-WGS defaults for these tuning parameters.
    // Applied here (not in nextflow.config) because Nextflow 26+ rejects top-level
    // if-statements in config files.
    def cna_window_size = params.run_lp_wgs ? 1000000                          : params.cna_window_size
    def ichorcna_ploidy = params.run_lp_wgs ? "c(2,3)"                         : params.ichorcna_ploidy
    def ichorcna_normal = params.run_lp_wgs ? "c(0.5,0.6,0.7,0.8,0.9,0.95)"    : params.ichorcna_normal
    def ichorcna_txn_e  = params.run_lp_wgs ? 0.9999                           : params.ichorcna_txn_e

    ch_fastqs = Channel.fromFilePairs(params.fastqs, checkIfExists: true)
        .map { key, reads ->
            def sample = key.replaceAll(/_L\d{3}$/, '')
            [ [id: key, sample: sample], reads ]
        }

    if (params.bwa_index) {
        ch_bwa_index = Channel.fromPath(params.bwa_index, checkIfExists: true).first()
    } else {
        PREPARE_GENOME(Channel.value(params.genome))
        ch_bwa_index = PREPARE_GENOME.out.bwa_index
    }

    // Resolve resource paths (fallback to auto-download URLs)
    def window_str = cna_window_size.toString()
    def res        = params.ichorcna_resources?.get(params.genome)
    def url_gc     = res?.get(window_str)?.gc_wig
    def url_map    = res?.get(window_str)?.map_wig
    def url_pon    = params.run_lp_wgs ? res?.get(window_str)?.pon_rds : 'NO_PON'
    def url_centro = res?.centromere

    def path_gc     = params.gc_wig         ?: url_gc
    def path_map    = params.map_wig        ?: url_map
    def path_pon    = params.pon_rds        ?: url_pon    ?: 'NO_PON'
    def path_centro = params.centromere_txt ?: url_centro ?: 'NO_CENTRO'

    if (!path_gc || !path_map) {
        error "Could not resolve GC/MAP WIG files. Provide them manually or use supported genome/window combinations (hg19/hg38 and 500kb/1Mb)."
    }

    ch_gc_wig         = Channel.fromPath(path_gc).first()
    ch_map_wig        = Channel.fromPath(path_map).first()
    ch_pon_rds        = path_pon == 'NO_PON' ? Channel.value(file('NO_PON')) : Channel.fromPath(path_pon).first()
    ch_centromere_txt = path_centro == 'NO_CENTRO' ? Channel.value(file('NO_CENTRO')) : Channel.fromPath(path_centro).first()

    QC(ch_fastqs)
    ALIGN(QC.out.reads, ch_bwa_index)

    CNA(
        ALIGN.out.bam_bai,
        cna_window_size,
        params.cna_min_mapq,
        ch_gc_wig,
        ch_map_wig,
        ch_pon_rds,
        ch_centromere_txt,
        params.genome,
        ichorcna_ploidy,
        ichorcna_normal,
        params.ichorcna_max_cn,
        ichorcna_txn_e,
        params.ichorcna_txn_strength
    )

    ch_multiqc = QC.out.json
        .mix(ALIGN.out.qc_stats)
        .mix(ALIGN.out.qc_flagstat)
        .mix(ALIGN.out.qc_idxstats)
        .collect()

    MULTIQC(ch_multiqc)
}