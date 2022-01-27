  
process qc {
    tag "QC from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/QC", mode:'copy'
    
    input:
    tuple val(pair_id), path(read1), path(read2)

    output:
    path ('*_fastqc.{zip,html}'), emit: report

    script:
    """
    fastqc ${read1} --threads $task.cpus -q
    fastqc ${read2} --threads $task.cpus -q
    """
}
process trimming {
    tag "Trimming from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/trim", mode: 'copy'

    input:
    tuple val(pair_id), path(read1), path(read2)

    output:
    tuple val(pair_id), path("*fq"), emit: trim_out
    path '*trimming_report.txt', emit: report
    path '*_fastqc.{zip,html}', emit: report_trim_qc

    script:
    """
    trim_galore --paired --fastqc --length ${params.length} \\
    -j $task.cpus ${read1} ${read2} \\
    --basename ${pair_id}
    """
}

process mapping {

    tag "Mapping from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/map", mode: 'copy'
    
    input:
    tuple val(pair_id), path(trim_out)
    path params.genomeDir
    
    output:
    tuple val(pair_id), path("${pair_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    path("${pair_id}_Aligned.sortedByCoord.out.bam.bai"), emit:bai
    path("*.out"), emit: report


    script:
    """
    STAR --readFilesIn ${trim_out[0]} ${trim_out[1]}  \\
    --genomeDir ${params.genomeDir} \\
    --outFileNamePrefix ${pair_id}_ \\
    --runThreadN ${task.cpus} \\
    --outSAMtype BAM SortedByCoordinate \\
    --peOverlapNbasesMin ${params.star_min_overlap} \\
    --alignIntronMax ${params.ailgn_intron_max} \\
    --peOverlapMMp ${params.star_max_overlap_mismatch}

    samtools index ${pair_id}_Aligned.sortedByCoord.out.bam
    samtools index -c ${pair_id}_Aligned.sortedByCoord.out.bam
    """
}
    // --readFilesCommand zcat \\