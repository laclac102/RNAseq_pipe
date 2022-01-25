process deseq2 {
    tag "DESeq2"
    publishDir "${baseDir}/ercc_samples/output/deseq2", mode: 'copy'

    input:
    path(merge)
    path(params.design)
    path(params.compare)

    output:
    path "*.{xlsx,jpg}"
    path "*_DESeq_results.tsv", emit: deseq2_result
    path "*{heatmap,plot,matrix}.tsv", emit: report
    path "v_DESeq2.txt"

    script:
    """
    DESeq2.r ${merge} ${params.design} $params.deseq2_fdr $params.deseq2_lfc ${params.compare}
    Rscript -e "write(x=as.character(packageVersion('DESeq2')), file='v_DESeq2.txt')"
    """
}