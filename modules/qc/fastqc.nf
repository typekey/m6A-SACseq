process FASTQC_BEFORE {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.QC_BEFORE}/${group}_${replicate}"

    input:
    tuple val(group), val(replicate), path(fastq_1), path(fastq_2)
    
    output:
    path '*.zip'

    script:
    """
    fastqc -f fastq -q ${fastq_1} ${fastq_2}
    """
}

process FASTQC_AFTER {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.QC_AFTER}/${sampleId}"

    input:
    path(cut_fastq1)
    path(cut_fastq2)
    // path(reprot)
    
    output:
    path '*.zip'

    script:
    { sampleId = cut_fastq1.toString().split("_")[0..1].join("_") }
    """
    fastqc -f fastq -q ${cut_fastq1} ${cut_fastq2}
    """
}

process FASTQC_UNMAPPED {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.QC_UNMAPPED}/${sampleId}"

    input:
    tuple path(unmapped_fastq1), path(unmapped_fastq2)
    // path(reprot)
    
    output:
    path '*.zip'

    script:
    { sampleId = unmapped_fastq1.toString().split("_")[0..1].join("_") }
    """
    fastqc -f fastq -q ${unmapped_fastq1} ${unmapped_fastq2}
    """
}