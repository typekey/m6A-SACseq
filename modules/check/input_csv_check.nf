process SAMPLES_CHECK {
    label 'process_single'

    input:
    path samplesheet

    output:
    path '*.csv'

    script:
    """
    echo TODO
    """
}