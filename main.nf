#!/usr/bin/env nextflow
/*
====================================================================
    m6A-SAC sequencing analysis pipeline
    Author: Lei Zheng
    Email : type.zheng@gmail.com
    Github: https://github.com/qqdb/m6A-SACseq
====================================================================
*/

nextflow.enable.dsl = 2

include { PREPARE } from './subworkflow/prepare.nf'
include { FROM_CSV } from './modules/import_data/from_csv'
// include { FALCO_BEFORE } from './modules/qc/falco'
// include { FALCO_AFTER } from './modules/qc/falco'
// include { FALCO_UNMAPPED } from './modules/qc/falco'
include { BOWTIE2_INDEX } from './modules/index/bowtie2.nf'
include { STAR_INDEX } from './modules/index/star.nf'
include { FASTQC_BEFORE } from './modules/qc/fastqc'
include { FASTQC_AFTER } from './modules/qc/fastqc'
include { FASTQC_UNMAPPED } from './modules/qc/fastqc'
include { RNASEQC } from './modules/qc/rnaseqc.nf'
include { BOWTIE2_MAPPING_CONTAMINATION } from './modules/mapping/bowtie2.nf'
include { BOWTIE2_MAPPING_SPIKE } from './modules/mapping/bowtie2.nf'
include { BOWTIE2_MAPPING_SNCRNA } from './modules/mapping/bowtie2.nf'
include { STAR_MAPPING } from './modules/mapping/star.nf'
include { SAMTOOLS_SORT_FILTER; SAMTOOLS_MERGE; SAMTOOLS_MERGE_GENOME; SAMTOOLS_INDEX } from './modules/sort_filter/samtools.nf'
include { UMICOLLAPSE_DD } from './modules/drop_duplicates/umicollapse.nf'
include { MULTIQC_REPORT } from './modules/report/multiqc_report'
include { CUTADAPT_ADAPTER } from './modules/trimming/cutadapt'
include { CUTADAPT_PRIMER } from './modules/trimming/cutadapt'
include { BGZIP_CPMPRESS } from './modules/compress/bgzip.nf'
include { SPIKEIN_ALIGNMENT } from './subworkflow/spikin-alignment.nf'

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
    // STAR_INDEX(params.references.GRCh38.fasta, params.references.GRCh38.gtf)


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
    
    // 3. sort and filter
    SAMTOOLS_SORT_FILTER(BOWTIE2_MAPPING_CONTAMINATION.out.sam | concat(BOWTIE2_MAPPING_SPIKE.out.sam) | concat(BOWTIE2_MAPPING_SNCRNA.out.sam))

    // 4.mapping genome
    STAR_MAPPING(BOWTIE2_MAPPING_SNCRNA.out.fastq_1, BOWTIE2_MAPPING_SNCRNA.out.fastq_2)
    // BGZIP_CPMPRESS(STAR_MAPPING.out.unmapped)

    // 5. unmapped-QC
    // FASTQC_UNMAPPED(STAR_MAPPING.out.unmapped)
    
    // 6. merge run bam
    SAMTOOLS_SORT_FILTER.out
        .map {
            file ->
            def group = file.name.toString().split("_")[0]
            def reftype = file.name.toString().split("_")[2].replace(".bam","")
            return [group + "_" + reftype, file]
        }
        .groupTuple()
        .set { group_bam_ch }
    SAMTOOLS_MERGE(group_bam_ch)
    // separate genome
    SAMTOOLS_MERGE_GENOME(STAR_MAPPING.out.bam)
    // UMICOLLAPSE_DD(SAMTOOLS_MERGE.out)

    SAMTOOLS_INDEX(SAMTOOLS_MERGE.out)

    // BAM QC
    RNASEQC(SAMTOOLS_MERGE_GENOME.out)
    
    genome_unmapped_fq = STAR_MAPPING.out.unmapped
    SPIKEIN_ALIGNMENT(SAMTOOLS_MERGE.out, genome_unmapped_fq)

    // report
    // cut_report = CUTADAPT_PRIMER.out[2] | collect
    // MULTIQC_REPORT(FASTQC_BEFORE.out | concat(cut_report) | concat(FASTQC_AFTER.out) | concat(STAR_MAPPING.out.log) | collect)
    // MULTIQC_REPORT(FASTQC_UNMAPPED.out | collect)
    // MULTIQC_REPORT(RNASEQC.out.tsv)
}

// prepare pipeline
// workflow {
//     PREPARE()
// }