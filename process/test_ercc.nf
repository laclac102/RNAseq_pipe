/*
*Scripts form Xingaulag
*Date: 10.01.2022 
*FastQC process 
*/

nextflow.enable.dsl=2

baseDir ="/home/xingau/workflow"
        params.reads="$baseDir/ercc_samples/samples/*.fastq.gz"
        params.length=20
        params.quality=20
        params.genome="$baseDir/ercc_samples/ref/chr22_ERCC92_transcripts.fa"
        params.out_Dir="$baseDir/ercc_samples/output"
        params.qc="${params.out_Dir}/QC"
        params.genomeDir="$baseDir/ercc_samples/index"
        params.bed12="baseDir/ercc_samples/output/ref/bed12/chr22_rRNA.bed12"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         length       : ${params.length}
         quality      : ${params.quality}
         Ref genome   : ${params.genome}
         Index genome : ${params.genomeDir}
         Output storage: ${params.out_Dir}
         """
         .stripIndent()

 Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set{read_pairs_ch}
    
process qc {
    tag "QC from $pair_id"
    storeDir="${params.out_Dir}/QC"
    
    input:
    tuple val(pair_id), path(reads)

    output:
    file("*.html")

    script:
    """
    fastqc ${reads} -o ${params.out_Dir}/QC --threads $task.cpus -q
    """
}
process trimming {
    tag "Trimming from $pair_id"
    storeDir="${params.out_Dir}/trim"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*_val_{1,2}.fq.gz"), emit: trim_out

    script:
    """
    trim_galore --paired --fastqc --retain_unpaired --length ${params.length} \\
    -j $task.cpus ${reads}
    
    """
}
process mapping {

    tag "Mapping from $pair_id"
    storeDir="${params.out_Dir}/map"
    
    input:
    tuple val(pair_id), file(trim_out)
    path genomeDir
    
    output:
    tuple val(pair_id), path("${pair_id}Aligned.sortedByCoord.out.bam"), path("*.{bai,csi}"), emit:bam

    script:
    """
    STAR --readFilesIn ${params.out_Dir}/trim/${pair_id}_1_val_1.fq.gz ${params.out_Dir}/trim/${pair_id}_2_val_2.fq.gz \\
    --genomeDir ${params.genomeDir} \\
    --outFileNamePrefix ${pair_id} \\
    --runThreadN ${task.cpus} \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand zcat

    samtools index ${pair_id}Aligned.sortedByCoord.out.bam
    samtools index -c ${pair_id}Aligned.sortedByCoord.out.bam
    """
}


workflow {
    qc(read_pairs_ch)
    trimming(read_pairs_ch.collect())
    mapping(trimming.out, params.genomeDir)
}