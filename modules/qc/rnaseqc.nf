process RNASEQC {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.QC_RNASEQC}"

    input:
    tuple path(input_bam1), path(input_bam2)
    
    
    output:
    path '*.metrics.tsv', emit: tsv
    // path '*.gct', emit: gct

    script:
    { 
        sampleId1 = input_bam1.toString().replace(".bam","")
        sampleId2 = input_bam2.toString().replace(".bam","")
    }
    """
    ${params.tools.rnaseqc} ${params.references.GRCh38.collapse_gtf} ${input_bam1} . -s ${sampleId1} --coverage
    ${params.tools.rnaseqc} ${params.references.GRCh38.collapse_gtf} ${input_bam2} . -s ${sampleId2} --coverage
    """
}

// ${params.tools.rnaseqc} ${params.references.GRCh38.collapse_gtf} ${input_bam1} rnaseqc -s ${sampleId1} --coverage -v
// ${params.tools.rnaseqc} ${params.references.GRCh38.collapse_gtf} ${input_bam2} rnaseqc -s ${sampleId2} --coverage -v