process COLLAPSE_GTF {
    debug params.debug
    label 'process_low'
    publishDir "${params.genome_base}"

    input:
    path input_gtf

    output:
    path "*.collapse.gtf"

    script:
    { base_name = input_gtf.name.replace(".gtf", ".collapse.gtf")}
    """
    python3 ${params.bin_base}/collapse_annotation.py ${input_gtf} ${base_name}
    """
}