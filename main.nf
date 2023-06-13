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
include { BOWTIE2_INDEX } from './modules/index/bowtie2.nf'
include { FASTQC_BEFORE } from './modules/qc/fastqc'
include { FASTQC_AFTER } from './modules/qc/fastqc'
include { BOWTIE2_MAPPING_CONTAMINATION } from './modules/mapping/bowtie2.nf'
include { BOWTIE2_MAPPING_SPIKE } from './modules/mapping/bowtie2.nf'
include { BOWTIE2_MAPPING_SNCRNA } from './modules/mapping/bowtie2.nf'
include { MULTIQC_REPORT } from './modules/report/multiqc_report'
include { CUTADAPT_ADAPTER } from './modules/trimming/cutadapt'
include { CUTADAPT_PRIMER } from './modules/trimming/cutadapt'

workflow {
    // { println 'Start m6A-SAC-seq pipeline: '}
    // 0. build index
    // Channel
    //     .fromPath([
    //         params.genome.contamination, 
    //         params.genome.genome_human, 
    //         params.genome.sncRNA_human, 
    //         params.genome.spike_expand])
    //     .view()
    //     .set{ genome_ch }
    // BOWTIE2_INDEX(genome_ch)

    // 1. QC
    fq_ch = FROM_CSV("./data/small_datasets/data.csv")
    // FASTQC_BEFORE(fq_ch)
    CUTADAPT_ADAPTER(fq_ch) \
    | CUTADAPT_PRIMER
    // FASTQC_AFTER(CUTADAPT_PRIMER.out)

    // 2. mapping
    BOWTIE2_MAPPING_CONTAMINATION(CUTADAPT_PRIMER.out)
    BOWTIE2_MAPPING_SPIKE(BOWTIE2_MAPPING_CONTAMINATION.out.fastq_1, BOWTIE2_MAPPING_CONTAMINATION.out.fastq_2)
    BOWTIE2_MAPPING_SNCRNA(BOWTIE2_MAPPING_SPIKE.out.fastq_1, BOWTIE2_MAPPING_SPIKE.out.fastq_2)
    // report
    // cut_report = CUTADAPT_PRIMER.out[2] | collect
    // MULTIQC_REPORT(FASTQC_BEFORE.out | concat(cut_report) | concat(FASTQC_AFTER.out) | collect, params.reprot_config)
    
}