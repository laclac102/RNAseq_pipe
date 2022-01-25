process multiqc {

    publishDir "${baseDir}/ercc_samples/output/multiqc", mode: 'copy'

    input:
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path(report_trim_qc)

    output:
    file("*")

    script:
    """
    multiqc . -m fastqc -m Trim_Galore -m star -m rseqc -m preseq -m picard -m qualimap -m featureCounts \\
    -m custom_content -m plot_sample_distance -m plot_gene_heatmap -m DESeq2 -m gProfiler
    """
}