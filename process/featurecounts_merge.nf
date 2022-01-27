process featurecounts {
    label 'low_memory'
    tag "FeatureCounts from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/featurecounts" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    path params.gtf

    output:
    path("${pair_id}_gene.featureCounts.txt"), emit: counts
    path "${pair_id}_gene.featureCounts.txt.summary", emit:report

    script:
    """
    featureCounts -a ${params.gtf} -F 'GTF' -g ${params.fc_group_features} -t 'exon' -o ${pair_id}_gene.featureCounts.txt --extraAttributes ${params.fc_extra_attributes} -p -s 1 ${bam} -T 1
    """
}

process merge_counts {
    publishDir "${baseDir}/ercc_samples/output/merge_counts", mode: 'copy'
    tag "Merge Featuresfor all samples"
    input:
    path(counts)

    output:
    path 'merged_gene_counts.txt', emit:merge
    file 'gene_lengths.txt'

    script:

    gene_ids = "<(tail -n +2 ${counts[0]} | cut -f1,7 )"
    merge_counts = counts.collect{filename ->
    "<(tail -n +2 ${filename} | sed 's:${params.bam_suffix}::' | cut -f8)"}.join(" ")
    """
    paste $gene_ids $merge_counts > merged_gene_counts.txt
    tail -n +2 ${counts[0]} | cut -f1,6 > gene_lengths.txt
    """
}