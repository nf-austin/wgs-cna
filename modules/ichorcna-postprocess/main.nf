process ICHORCNA_POSTPROCESES {
    tag "all samples"
    publishDir "${params.outdir}/cna/ichorcna_postprocess", mode: 'copy'

    input:
    path params_txts  // collected *.params.txt from all samples
    path seg_txts     // collected *.seg.txt from all samples
    path pdfs         // collected PDF figures from all samples

    output:
    path "ichorCNA_summary.tsv",  emit: summary
    path "ichorCNA_report.html",  emit: report

    script:
    """
    export PATH=/opt/conda/bin:\$PATH
    collect_ichorcna_data.py \\
        --params ${params_txts} \\
        --segs   ${seg_txts} \\
        --pdfs   ${pdfs} \\
        --output-tsv  ichorCNA_summary.tsv \\
        --output-html ichorCNA_report.html
    """
}