
process mark_dup {
    tag "Qualimap from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/markdup" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    

    output:
    tuple val(pair_id), path("${pair_id}.markDups.bam"), path("${pair_id}.markDups.bam.{bai,csi}")
    path "${pair_id}*", emit: report

    script:
    """
    picard MarkDuplicates \\
            INPUT=${bam} \\
            OUTPUT=${pair_id}.markDups.bam \\
            METRICS_FILE=${pair_id}.markDups_metrics.txt \\
            REMOVE_DUPLICATES=false \\
            ASSUME_SORTED=true \\
            PROGRAM_RECORD_ID='null' \\
            VALIDATION_STRINGENCY=LENIENT
        samtools index -c ${pair_id}.markDups.bam
        samtools index ${pair_id}.markDups.bam 
    """
}
process preseq {
    tag "Qualimap from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/preseq" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    
    output:
    path "${pair_id}*", emit: report

    script:
    """
    preseq lc_extrap -v -D -B ${bam} -o ${pair_id}.ccurve.txt
    """
}