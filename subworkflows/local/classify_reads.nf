/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: CLASSIFY READS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLASSIFY_READS {

    take:
    ch_reads
    ch_kraken2_db

    main:

    ch_versions = Channel.empty()

    KRAKEN2_KRAKEN2(
        ch_reads,
        ch_kraken2_db,
        params.save_output_fastqs ?: false,
        params.save_reads_assignment ?: false
    )

    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions_kraken2)
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions_pigz)

    emit:
    kraken2_report      = KRAKEN2_KRAKEN2.out.report
    kraken2_assignment  = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    classified_fastq    = KRAKEN2_KRAKEN2.out.classified_reads_fastq
    unclassified_fastq  = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
    versions            = ch_versions
}
