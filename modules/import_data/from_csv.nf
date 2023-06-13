workflow FROM_CSV{
    take:
    csv_path

    main:
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true)
        .map {
            row ->
            tuple(row.group, row.replicate, file(row.fastq_1), file(row.fastq_2))
        }
        .set { fq_ch }
    
    emit:
    data_ch = fq_ch
}