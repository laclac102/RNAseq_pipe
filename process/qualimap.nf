process qualimap {
    tag "Qualimap from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/qualimap" , mode:'copy'
    
    input:
    tuple val(pair_id), path(bam)
    path(params.gtf)

    output:
    path("${pair_id}_Aligned.sortedByCoord.out_rnaseq_qc/${pair_id}.pdf") , emit: report

    script:
    """
    qualimap rnaseq -bam ${bam} -gtf ${params.gtf} -outfile ${pair_id}.pdf -p strand-specific-forward
    """
}