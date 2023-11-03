process BGZIP_CPMPRESS {
     publishDir "${params.outdir.MAPPING}/genome/Unmapped", pattern: "*.fq.gz"

    input:
    tuple path(input_file1), path(input_file2)

    output:
    path "*"

    script:
    """
    bgzip -@ ${params.threads.base} -l 9 ${input_file1} > ${input_file1}.fq.gz
    bgzip -@ ${params.threads.base} -l 9 ${input_file2} > ${input_file2}.fq.gz
    """
}