nextflow.enable.dsl=2

baseDir ="/home/xingau/workflow"
params.reads="$baseDir/input/*_{1,2}.fq"
params.length=20
params.quality=20
params.genome="$baseDir/ref/transcriptome.fa"
params.out_Dir="$baseDir/output"
params.genomeDir="$baseDir/index"

process index {
    tag "Index ref Genome"
    storeDir="${params.genomeDir}"

    input:
    path genome

    output:
    file("*")

    script:
    """
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeFastaFiles ${params.genome} --genomeDir ${params.genomeDir}
    """
}
workflow {

    index(Channel.from(params.genome))
}


