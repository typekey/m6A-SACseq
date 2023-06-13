process CUTADAPT_ADAPTER {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.TRIMMING}"

    input:
    tuple val(group), val(replicate), path(fastq_1), path(fastq_2)

    output:
    path "*_cut_R1.fastq"
    path "*_cut_R2.fastq"
    // path "*_trimming_adapter.report", optional:true

    script:
    """
    cutadapt -j ${params.threads.cutadapt} \
        -U 11 \
        --rename='{id}_{r1.cut_prefix} {comment}' \
        --max-n=0 -e 0.15 -q 20 --nextseq-trim=20 \
        -O 6 \
        --pair-filter=both \
        -a ${params.barcode.barcode3_r1} -A ${params.barcode.barcode3_r2} \
        -o ${group}_${replicate}_cut_R1.fastq -p${group}_${replicate}_cut_R2.fastq \
        ${fastq_1} ${fastq_2} >${group}_${replicate}_trimming_adapter.report 
    """
}


process CUTADAPT_PRIMER {
    debug params.debug
    label 'process_low'

    publishDir "${params.outdir.TRIMMING}"

    input:
    path(cut_R1)
    path(cut_R2)
    // path(trimming_adapter)

    output:
    path "*_cut_R1.fastq.gz"
    path "*_cut_R2.fastq.gz"
    // path "*_trimming_primer.report", optional:true
    // path "cut_adapter/${group}_cut_R2.fq.gz"
    // path "cut_adapter/${group}_short_R1.fq.gz"
    // path "cut_adapter/${group}_short_R2.fq.gz"
    // path "cut_adapter/${group}_step1_report"
    // path "cut_adapter/${group}_step2_report"

    script:
    { sampleId = cut_R1.name.toString().replace("_cut_R1.fastq", "") }
    """
    cutadapt -j ${params.threads.cutadapt} \
        -m 15 \
        -u -11 \
        -n 2 \
        -O 12 \
        -g ${params.barcode.primerF} -a ${params.barcode.primerR} \
        -G ${params.barcode.primerF} -A ${params.barcode.primerR} \
        --too-short-output=${sampleId}_short_R1.fastq --too-short-paired-output=${sampleId}_short_R2.fastq \
        -o ${sampleId}_cut_R1.fastq.gz -p ${sampleId}_cut_R2.fastq.gz \
        ${sampleId}_cut_R1.fastq ${sampleId}_cut_R2.fastq >${sampleId}_trimming_primer.report 
    """
}
