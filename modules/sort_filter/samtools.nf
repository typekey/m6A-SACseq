process SAMTOOLS_SORT_FILTER {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}"

    input:
    path input_sam

    output:
    path '*.bam'

    script:
    { 
        sampleId = input_sam.name.toString().replace("sam", "bam") 
        reftype =  input_sam.name.toString().split("_")[2].replace(".sam", "")
    }
    """
    ${params.tools.samtools} sort -@ ${params.threads.samtools}  --input-fmt-option 'filter=[NM]<=10' -m 2G -O BAM -o ${sampleId} ${input_sam}
    """
}

process SAMTOOLS_MERGE {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/combined_mapping"

    input:
    tuple val(group), path(input_bam)

    output:
    path '*'

    script:
    """
    ${params.tools.samtools} merge -@ ${params.threads.base} -o ${group}.bam ${input_bam}
    """
}

process SAMTOOLS_MERGE_GENOME {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/separate_genome_combined"

    input:
    path genome_bam

    output:
    path '*'

    script:
    { group = genome_bam.name.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.samtools} view -@ ${params.threads.base} --input-fmt-option 'filter=[NH]==1' -h ${genome_bam} -O BAM -o ${group}_unique.bam
    ${params.tools.samtools} view -@ ${params.threads.base} --input-fmt-option 'filter=[NH]>1' -h ${genome_bam} -O BAM -o ${group}_multi.bam
    """
}

process SAMTOOLS_INDEX {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/combined_mapping"

    input:
    path input_bam

    output:
    path '*'

    script:
    """
    ${params.tools.samtools} index -@ ${params.threads.base} ${input_bam}
    """
}