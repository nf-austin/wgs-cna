include { HMMCOPY_READCOUNTER    } from '../modules/hmmcopy/readcounter/main'
include { ICHORCNA               } from '../modules/ichorcna/main'
include { ICHORCNA_POSTPROCESES  } from '../modules/ichorcna-postprocess/main'

workflow CNA {
    take:
    bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    window_size    // val
    min_mapq       // val
    gc_wig         // path
    map_wig        // path
    pon_rds        // path
    centromere_txt // path
    genome_build   // val
    ploidy         // val
    normal         // val
    max_cn         // val
    txn_e          // val
    txn_strength   // val

    main:
    HMMCOPY_READCOUNTER(bam_bai, window_size, min_mapq)

    ICHORCNA(
        HMMCOPY_READCOUNTER.out.wig,
        gc_wig,
        map_wig,
        pon_rds,
        centromere_txt,
        genome_build,
        ploidy,
        normal,
        max_cn,
        txn_e,
        txn_strength
    )

    ICHORCNA_POSTPROCESES(
        ICHORCNA.out.params_txt.map { meta, f -> f }.collect(),
        ICHORCNA.out.seg.map        { meta, f -> f }.collect(),
        ICHORCNA.out.plots.map      { meta, f -> f }.collect()
    )

    emit:
    cna        = ICHORCNA.out.cna
    seg        = ICHORCNA.out.seg
    params_txt = ICHORCNA.out.params_txt
    plots      = ICHORCNA.out.plots
    summary    = ICHORCNA_POSTPROCESES.out.summary
    report     = ICHORCNA_POSTPROCESES.out.report
}