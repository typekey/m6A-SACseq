process STAR_MAPPING {
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/genome", pattern: "*.bam", saveAs: { filename -> "${sampleId}_genome.bam" }
    publishDir "${params.outdir.MAPPING}/genome/log", pattern: "*.out"
    publishDir "${params.outdir.MAPPING}/genome/SJ", pattern: "*SJ*"
    publishDir "${params.outdir.MAPPING}/genome/Unmapped", pattern: "*Unmapped*"

    input:
    path fastq_1
    path fastq_2

    output:
    path "*.bam", emit: bam
    path "*.out", emit: log
    path "*SJ*", emit: sj
    path "*Unmapped*", emit: unmapped

    script:
    { sampleId = fastq_1.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.star} \
      --runThreadN ${params.threads.bowtie2} \
      --genomeDir ${params.references.GRCh38.star} \
      --readFilesIn ${fastq_1} ${fastq_2} \
      --alignEndsType Local \
      --outFilterMatchNminOverLread 0.66 \
      --outFilterMatchNmin 15 \
      --outFilterMismatchNmax 5 \
      --outFilterMismatchNoverLmax 0.2 \
      --outFilterMultimapNmax 50 \
      --outSAMmultNmax -1 \
      --outReadsUnmapped Fastx \
      --outSAMattrRGline ID:${sampleId} SM:${sampleId} LB:RNA PL:Illumina PU:SE \
      --outSAMattributes NH HI AS nM NM MD jM jI MC \
      --limitBAMsortRAM 8000000000 \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${sampleId}_genome_
    """
}