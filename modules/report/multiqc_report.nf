process MULTIQC_REPORT {
    debug params.debug
    label 'process_low'
    publishDir "${params.outdir.REPORT}"

    input:
    path zip

    output:
    path '*'
    

    script:
    """
    cp ${params.reprot_config}/* .
    echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
    multiqc $zip
    """
}
// # echo multiqc -f -m fastqc -n ${report_name} ${fastqc_data}
//     multiqc -f -m fastqc -n ${report_name} ${fastqc_data}
// cp $config/* .
//     echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
//     multiqc $zip