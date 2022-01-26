// params.ignore_R1 = false

// Sanity check design file
process check_design {
    tag "$design"
    publishDir "${baseDir}/ercc_samples/output/check_design", mode:'copy'

    input:
    path design
    path compare
    // path kitqcplotdatasets

    output:
    path "checked_${design}", emit: checked_design

    script:
    comparison_file = params.compare ? "-c $compare" : ''
    // kitqcplotset_file = params.kit_QC ? "-s $kitqcplotdatasets" : ''
    // ignorer1 = params.ignore_R1 ? "--ignore_r1" : "--no_ignore_r1"
    """
    check_design.py $comparison_file $design 
    """
}
// check_design.py $comparison_file $kitqcplotset_file $design $ignorer1