include { COLLAPSE_GTF } from '../modules/collapse_gtf/collapse_gtf.nf'
include { MAKE_BLASTDB } from '../modules/blast_db/blast.nf'

workflow PREPARE{
    // collapse annotation
    // COLLAPSE_GTF(params.references.GRCh38.gtf)

    MAKE_BLASTDB(params.references.spike_degenerate.fasta)
}