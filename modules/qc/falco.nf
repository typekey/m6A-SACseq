process FALCO_BEFORE {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.QC_BEFORE}/${group}_${replicate}"

    input:
    tuple val(group), val(replicate), path(fastq_1), path(fastq_2)
    
    output:
    path '*'

    script:
    """
    falco ${fastq_1} ${fastq_1}
    """
}

process FALCO_AFTER {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.QC_AFTER}/${sampleId}"

    input:
    path(cut_fastq1)
    path(cut_fastq2)
    path(reprot)
    
    output:
    path '*'

    script:
    { sampleId = cut_fastq1.toString().split("_")[0..1].join("_") }
    """
    falco ${cut_fastq2}
    """
}

// process FALCO_AFTER {
//     debug params.debug
//     label 'process_low'

//     publishDir "${params.outdir.QC_BEFORE}"

//     input:
//     tuple val(group), val(replicate), path(fastq_1), path(fastq_2)
    
//     output:
//     '*'

//     script:
//     """
//     falco ${fastq_1}
//     """
// }