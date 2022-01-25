process dupradar {
    tag "Dupradar from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/dupradar" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    path params.gtf

    output:
    path "*.pdf"
    path "*_dupMatrix.txt"
    path "*_intercept_slope.txt"
    path "*_mqc.txt", emit: report
    path "v_dupRadar.txt"

    script:
    """
    dupRadar.r ${bam} ${params.gtf} 0 paired 1
    Rscript -e "write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
    """
}