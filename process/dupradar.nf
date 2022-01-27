process dupradar {
    label 'mid_memory'
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

    script:
    """
    dupRadar.r ${bam} ${params.gtf} ${params.strandedness} paired 1
    Rscript -e "write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
    """
}

// Usage: dupRadar.r <input.bam> <annotation.gtf> <strandDirection:0=unstranded/1=forward/2=reverse> <paired/single> <nbThreads> <R-package-location (optional)>