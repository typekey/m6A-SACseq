
// 从 mapped spike bam 中抽取 fq1,fq2，并合并两个 fastq 文件
process MERGE_MAPPED_SPIKE {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.TMP}"

    input:
    path spike_bam

    output:
    path '*_spike.fq', emit: spike_fq

    script:
    { sample_id = spike_bam.name.toString().replace(".bam","") }
    """
    ${params.tools.samtools} fastq -@ ${params.threads.base} -1 ${sample_id}_R1.fq -2 ${sample_id}_R2.fq -0 /dev/null -s /dev/null -n ${spike_bam} 
    ${params.tools.fastp} -w ${params.threads.base} -j /dev/null -h /dev/null -i ${sample_id}_R1.fq -I ${sample_id}_R2.fq -m --merged_out ${sample_id}.fq
    """
}

// 直接合并 unmapped genome 的fq1，fq2，为一个 fastq 文件
process MERGE_UNMAPPED_GENOME {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.TMP}"

    input:
    tuple path(genome_unmapped_fq1), path(genome_unmapped_fq2)

    output:
    path '*'

    script:
    { sample_id = genome_unmapped_fq1.name.toString().split("_")[0..2].join("_") + "_unmapped" }
    """
    ${params.tools.fastp} -w ${params.threads.base} -j /dev/null -h /dev/null  -i ${genome_unmapped_fq1} -I ${genome_unmapped_fq2}  -m --stdout | grep --no-group-separator -A 2 -B 1  -P "CTAGAATTACACCA|TGGTGTAATTCTAG" || true; > ${sample_id}.fq
    """
}

// 合并重复样本的 fastq
process COMBINED_UNMAPPED {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.TMP}"

    input:
    tuple val(sample_id), path(genome_unmapped_fqs)

    output:
    path '*'

    script:
    """
    cat ${genome_unmapped_fqs} > ${sample_id}.fq
    """
}

//  
process MEGRE_SPIKE_FASTA {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.TMP}/spike_merge"

    input:
    path mapped_spike_fq
    path unmapped_genome_fq

    output:
    path '*.fq', emit: fq
    path '*.fa', emit: fa

    script:
    { sample_id = mapped_spike_fq.name.toString().replace(".fq", "") + "_merge"}
    """
    cat ${mapped_spike_fq} ${unmapped_genome_fq}  | tee ${sample_id}.fq | seqtk seq -a > ${sample_id}.fa
    """
}

process MAP_TO_SPIKIN {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.ALIGNMENT_SPIKE}"

    input:
    path spike_fa

    output:
    path '*.xml', emit: xml

    script:
    { sample_id = spike_fa.name.toString().replace(".fa", "") }
    """
    blastn -num_threads ${params.threads.base} -max_target_seqs 1 -db ${params.references.spike_degenerate.blastdb} -query ${spike_fa} -outfmt 5 > ${sample_id}.xml
    """
}

// TODO: blast2bam 需要编译，或拉取 conda 和 docker
process BLAST_TO_BAM {
    debug params.debug
    label 'process_low'
    
    input:
    path input_xml
    path input_fq

    output:
    path '*.bam'

    script:
    { output_bam = input_xml.name.toString().replace(".xml", ".bam") }
    """
    ${params.tools.blast2bam} ${input_xml} ${params.references.spike_degenerate.fasta} ${input_fq} | \
     samtools calmd -@ ${params.threads.base} --input-fmt-option 'filter=pos < 10 && pos + qlen > 33 && !flag.unmap' --output-fmt BAM - ${params.references.spike_degenerate.fasta} 2>/dev/null > ${output_bam}
    """
}

process SORT_BAM {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.ALIGNMENT_SPIKE}"

    input:
    path input_unsort_bam

    output:
    path sort_bam

    script:
    { output_bam =  = input_xml.name.toString().replace(".bam", "_sort.bam") }
    """
    ${params.tools.samtools} sort -@ {params.threads.base} -m 12G --write-index ${input_unsort_bam} -o {output_bam}
    """
}

process CALL_MUTATION_OF_SPIKE {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.ALIGNMENT_SPIKE}"

    input:
    path input_sort_bam

    output:
    path '*'

    script:
    { 
        header = "\t".join(["ref", "motif", "base", "count"])
        sample_id =  = input_xml.name.toString().replace(".bam", "")
    }
    """
    (
       echo ${header}
       {params.tools.call_spike_mutation} ${input_sort_bam} | awk '{{ t[$0]++ }} END{{ for (i in t) print t[i],i }}' | awk 'BEGIN{{OFS="\\t"}}{{print $2,$3,$4,$1}}'
    ) | bgzip -@ {params.threads.base} -l 9 > ${sample_id}.tsv.gz
    """
}



workflow SPIKEIN_ALIGNMENT{
    take: 
    merge_bam
    genome_unmapped_fq
    // path genome_unmapped_fq

    main:
    // 
    merge_bam.map {
            file ->
            def reftype = file.name.toString().split("_")[1].replace(".bam","")
            if (reftype == "spike") return file
        }
        .set{ spike_bam }

    MERGE_MAPPED_SPIKE(spike_bam)
    MERGE_UNMAPPED_GENOME(genome_unmapped_fq)

    MERGE_UNMAPPED_GENOME
        .out
        .map {
            fastq -> 
            def sample_id = fastq.name.toString().split("_")[0]
            return [sample_id + '_genome_unmapped', fastq]
        }
        .groupTuple()
        .set {unmapped_fq_ch}

    COMBINED_UNMAPPED(unmapped_fq_ch)

    MEGRE_SPIKE_FASTA(MERGE_MAPPED_SPIKE.out, COMBINED_UNMAPPED.out)
    MAP_TO_SPIKIN(MEGRE_SPIKE_FASTA.out.fa)

    BLAST_TO_BAM(MAP_TO_SPIKIN.out.xml, MEGRE_SPIKE_FASTA.out.fq)
}