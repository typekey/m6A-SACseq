process STAR_INDEX {
    debug params.debug
    publishDir "${params.genome_base}/star_index", pattern: "GenomeDir/*"
    // publishDir "${params.genome_base}/star_index/log", pattern: "*.log"
    label "process_high"

    input:
    path genome
    path gtf

    output:
    path "*"

    script:
    {index_name = genome.name.toString().replace(".fa", "")}
    """
    STAR \
    --runMode genomeGenerate\
    --runThreadN ${params.threads.star}\
    --genomeFastaFiles  ${genome}\
    --sjdbGTFfile ${gtf} 1>>${index_name}.log 2>>${index_name}.log
    """
}