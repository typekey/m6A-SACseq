process MAKE_BLASTDB {
    debug params.debug
    label 'process_low'
    publishDir "${params.references.blastdb_base}"

    input:
    path spike_fa

    output:
    path '*'

    script:
    { db_name = params.references.spike_degenerate.blastdb.toString().split("/").last()}
    """
    ${params.tools.makeblastdb} -in ${spike_fa} -dbtype 'nucl' -out ${db_name} 2>&1 >> ${db_name}.log
    """
}
