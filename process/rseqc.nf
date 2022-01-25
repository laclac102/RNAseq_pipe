process rseqc {
    tag "RNA seq QC from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/rseqc" , mode:'copy'
    
    input:
    tuple val(pair_id), path(bam)
    path(bai)
    path(params.bed12)

    output:
    path "${pair_id}*" , emit: report

    script:
    """
    read_distribution.py -i ${bam} -r ${params.bed12} > "${pair_id}.rseqc.read_distribution.txt"
    infer_experiment.py -i ${bam} -r ${params.bed12}
    inner_distance.py -i ${bam} -o ${pair_id}.rseqc -r ${params.bed12}
    read_duplication.py -i ${bam} -o ${pair_id}.rseqc.read_duplication
    junction_annotation.py -i ${bam} -r ${params.bed12} -o ${pair_id}.rseqc.junction_annotation
    junction_saturation.py -i ${bam} -o ${pair_id}.rseqc -r ${params.bed12}
    bam_stat.py -i ${bam} > ${pair_id}.rseqc.bam_stat.txt
    """
}