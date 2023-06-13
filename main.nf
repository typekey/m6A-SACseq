#!/usr/bin/env nextflow
/*
====================================================================
    m6A-SAC sequencing analysis pipeline
    Author: Lei Zheng
    Email : baimoc@163.com
    Github: https://github.com/qqdb/m6A-SACseq
====================================================================
*/

nextflow.enable.dsl = 2

include { FROM_CSV } from './modules/import_data/from_csv'
// include { FALCO_BEFORE } from './modules/qc/falco'
// include { FALCO_AFTER } from './modules/qc/falco'
include { FASTQC_BEFORE } from './modules/qc/fastqc'
include { FASTQC_AFTER } from './modules/qc/fastqc'
include { MULTIQC_REPORT } from './modules/report/multiqc_report'
include { CUTADAPT_ADAPTER } from './modules/trimming/cutadapt'
include { CUTADAPT_PRIMER } from './modules/trimming/cutadapt'


workflow QC{
    main:
    { println 'Start m6A-SAC-seq pipeline: '}

    fq_ch = FROM_CSV("./test_data/small_datasets/data.csv")
    // fq_ch
    //     .flatMap {[it[2], it[3]]}
    //     .view()
    FASTQC_BEFORE(fq_ch)
    // MULTIQC_REPORT(FALCO_BEFORE.out[0], "${params.publish_dir}/quality_control",'report_falco_before.html')
    CUTADAPT_ADAPTER(fq_ch) | CUTADAPT_PRIMER
    FASTQC_AFTER(CUTADAPT_PRIMER.out)
    emit:
    FASTQC_AFTER.out
}


workflow {
    // { println 'Start m6A-SAC-seq pipeline: '}

    fq_ch = FROM_CSV("./test_data/small_datasets/data.csv")
    FASTQC_BEFORE(fq_ch)
    CUTADAPT_ADAPTER(fq_ch) | CUTADAPT_PRIMER
    FASTQC_AFTER(CUTADAPT_PRIMER.out)
    cut_report = CUTADAPT_PRIMER.out[2] | collect

    MULTIQC_REPORT(FASTQC_BEFORE.out | concat(cut_report) | concat(FASTQC_AFTER.out) | collect, params.reprot_config)
}