process gprofiler {
    publishDir "${baseDir}/ercc_samples/output/gprofiler", mode: 'copy'

    input:
    path(deseq2_result)

    output:
    path "*_gProfiler_results.tsv", emit: report
    path "*_gProfiler_results.xlsx"
    path "v_gProfiler.txt"

    script:
    """
    gProfiler.py ${deseq2_result} -o $params.gprofiler_organism -q $params.deseq2_fdr -p $params.gprofiler_fdr
    pip freeze | grep gprofiler > v_gProfiler.txt
    """
}