process BOWTIE2_MAPPING {
    debug params.debug
    label 'process_high'

    input:
    path fastq_1
    path fastq_2
    val index_name

    // output:
    // path '*'

    script:
    """
    {params.path_bowtie2} -p {params.threads} \
    --no-unal --end-to-end --fast \
    --un-conc {params.un} -x {params.ref_bowtie2} \
    -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
    """
}

process BOWTIE2_MAPPING_CONTAMINATION {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/contamination"

    input:
    path fastq_1
    path fastq_2
    
    output:
    path '*.1.fq', emit: fastq_1
    path '*.2.fq', emit: fastq_2
    path '*.sam', emit: sam
    path '*.report'

    script:
    { sampleId = fastq_1.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.bowtie2} -p ${params.threads.bowtie2} \
    --no-unal --end-to-end --fast \
    --un-conc ${sampleId}_contamination.fq -x ${params.references.contamination.bowtie2} \
    -1 ${fastq_1} -2 ${fastq_2} > ${sampleId}_contamination.sam 2> >(tee ${sampleId}_output.report >&2)
    """
}

process BOWTIE2_MAPPING_SPIKE {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/spike"

    input:
    path fastq_1
    path fastq_2

    output:
    path '*.1.fq', emit: fastq_1
    path '*.2.fq', emit: fastq_2
    path '*.sam', emit: sam
    path '*.report'

    script:
    { sampleId = fastq_1.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.bowtie2} -p ${params.threads.bowtie2} \
    --nofw --no-unal --end-to-end -L 16 -N 1 --mp 5 \
    --un-conc ${sampleId}_spike.fq -x ${params.references.spike_expand.bowtie2} \
    -1 ${fastq_1} -2 ${fastq_2} > ${sampleId}_spike.sam 2> >(tee ${sampleId}_output.report >&2)
    """
}

process BOWTIE2_MAPPING_SNCRNA {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/sncRNA"

    input:
    path fastq_1
    path fastq_2

    output:
    path '*.1.fq', emit: fastq_1
    path '*.2.fq', emit: fastq_2
    path '*.sam', emit: sam
    path '*.report'

    script:
    { sampleId = fastq_1.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.bowtie2} -p ${params.threads.bowtie2} \
    --nofw --all --no-unal --end-to-end -L 16 -N 1 --mp 5 \
    --un-conc ${sampleId}_sncRNA.fq -x ${params.references.sncRNA_human.bowtie2} \
    -1 ${fastq_1} -2 ${fastq_2} > ${sampleId}_sncRNA.sam 2> >(tee ${sampleId}_output.report >&2)
    """
}
