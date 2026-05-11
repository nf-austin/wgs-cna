process HMMCOPY_READCOUNTER {
    tag "$meta.id"
    publishDir "${params.outdir}/cna/readcounts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    val window_size
    val min_mapq

    output:
    tuple val(meta), path("*.wig"), emit: wig

    script:
    """
    readCounter \\
        --window ${window_size} \\
        --quality ${min_mapq} \\
        ${bam} > ${meta.id}.bin${window_size}.wig
    """
}