nextflow.enable.dsl=2

baseDir='/mnt/d/RNAseq'
    params.length=20
    params.quality=20
    params.genomeDir="${baseDir}/ercc_samples/ref/index"
    params.bed12="${baseDir}/ercc_samples/ref/bed12/chr22_with_ERCC92.bed12"
    params.csvDir ="${baseDir}/metadata/ercc_fullmeta.csv"
    params.gtf="${baseDir}/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"
    params.bam_suffix = "_Aligned.sortedByCoord.out.bam"
    params.design="${baseDir}/metadata/design.csv"
    params.compare="${baseDir}/metadata/comparison.csv"
    params.deseq2_fdr = 0.05
    params.deseq2_lfc = 0.585
    params.gprofiler_fdr = 0.05
    params.gprofiler_organism = 'hsapiens'

meta = Channel.from(file(params.csvDir))
    .splitCsv(header:true)
    .map{ row-> tuple("$row.pair_id"), file("$row.read1"), file("$row.read2") }
    .set{sample_ch}

include { check_design } from './process/check_design'
include { qc; trimming; mapping } from './process/map'
include { rseqc } from ('./process/rseqc')
include { qualimap } from ('./process/qualimap')
include { mark_dup; preseq } from ('./process/markdup_preseq')
include { dupradar } from ('./process/dupradar')
include { featurecounts; merge_counts } from ('./process/featurecounts_merge')
include { deseq2 } from ('./process/DESeq2.nf')
include { gprofiler } from ('./process/gprofiler')
include { multiqc } from ('./process/multiqc')

workflow {
    check_design(params.design, params.compare)
    qc(sample_ch)
    trimming(sample_ch)
    mapping(trimming.out.trim_out, params.genomeDir)
    rseqc(mapping.out.bam, mapping.out.bai, params.bed12)
    qualimap(mapping.out.bam, params.gtf)
    mark_dup(mapping.out.bam, mapping.out.bai)
    preseq(mapping.out.bam, mapping.out.bai)
    dupradar(mapping.out.bam, mapping.out.bai, params.gtf)
    featurecounts(mapping.out.bam, mapping.out.bai, params.gtf)
    merge_counts(featurecounts.out.counts.collect())
    deseq2(merge_counts.out.merge, params.design, params.compare)
    gprofiler(deseq2.out.deseq2_result)
    
    multiqc(qc.out.report.collect(), \
    trimming.out.report.collect(), \
    trimming.out.report_trim_qc.collect(), \
    mapping.out.report.collect(), \
    rseqc.out.report.collect(), \
    mark_dup.out.report.collect(), \
    preseq.out.report.collect(), \
    dupradar.out.report.collect(),\
    qualimap.out.report.collect(), \
    featurecounts.out.report.collect(), \
    deseq2.out.report.collect(), \
    gprofiler.out.report.collect())
} 