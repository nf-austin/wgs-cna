process ICHORCNA {
    tag "$meta.id"
    publishDir "${params.outdir}/cna/ichorcna", mode: 'copy'

    input:
    tuple val(meta), path(wig)
    path gc_wig
    path map_wig
    path pon_rds
    path centromere_txt
    val genome_build
    val ploidy
    val normal
    val max_cn
    val txn_e
    val txn_strength

    output:
    tuple val(meta), path("*.cna.seg"), emit: cna
    tuple val(meta), path("*.seg.txt"), emit: seg
    tuple val(meta), path("*.params.txt"), emit: params_txt
    tuple val(meta), path("*.pdf"), emit: plots
    path "*.RData", emit: rdata

    script:
    def pon_arg = pon_rds.name != 'NO_PON' ? "--normalPanel ${pon_rds}" : ""
    def centro_arg = centromere_txt.name != 'NO_CENTRO' ? "--centromere ${centromere_txt}" : ""

    """
    runIchorCNA.R \\
        --id ${meta.id} \\
        --WIG ${wig} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        --genomeBuild ${genome_build} \\
        ${pon_arg} \\
        ${centro_arg} \\
        --ploidy "${ploidy}" \\
        --normal "${normal}" \\
        --maxCN ${max_cn} \\
        --txnE ${txn_e} \\
        --txnStrength ${txn_strength} \\
        --includeHOMD False \\
        --estimateNormal True \\
        --estimatePloidy True \\
        --estimateScPrevalence True \\
        --outDir ./
    """
}