process BOWTIE2_INDEX {
    debug params.debug
    publishDir "${params.genome_base}/bowtie2_index", pattern: "*.bt2"
    publishDir "${params.genome_base}/bowtie2_index/log", pattern: "*.log"
    label "process_high"

    input:
    path genome

    output:
    path "*"

    script:
    {index_name = genome.name.toString().replace(".fa", "")}
    """
    ${params.tools.bowtie2_build} $genome $index_name > "${index_name}.log" 2>&1
    """
}