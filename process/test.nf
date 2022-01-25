/*
*Scripts form Xingaulag
*Date: 10.01.2022 
*FastQC process 
*/

nextflow.enable.dsl=2

baseDir ="/home/xingau/workflow"
params.reads="$baseDir/input/*_{1,2}.fq"
params.length=20
params.quality=20
params.genome="$baseDir/ref/transcriptome.fa"
params.out_Dir="$baseDir/output"
params.genomeDir="$baseDir/index"

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
    input:
    tuple val(pair_id), path(reads)

    output:
    file("*")

    script:
    """
    mkdir QC_pair_$pair_id
    fastqc ${reads} -o QC_pair_$pair_id --threads $task.cpus -q
    """
}
process trimming {
    tag "Trimming from $pair_id"
    storeDir="${params.out_Dir}"\

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*_val_{1,2}.fq"), emit: trim_out

    script:
    """
    trim_galore --paired --fastqc --retain_unpaired --length ${params.length} -j $task.cpus ${reads}
    """
}

workflow {

    qc(read_pairs_ch)
    trimming(read_pairs_ch)
}