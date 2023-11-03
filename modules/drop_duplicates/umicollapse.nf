process UMICOLLAPSE_DD {
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/drop_duplicates"

    input:
    path input_bam

    output:
    path '*'

    script:
    { group = input_bam.name.toString().split("_")[0..1].join("_") }
    """
    java -server -Xms4G -Xmx64G -Xss100M -Djava.io.tmpdir=/tmp -jar ${params.tools.umicollapse} bam \
            --two-pass -i ${input_bam} -o ${group}.bam  > ${group}.log
    """
}