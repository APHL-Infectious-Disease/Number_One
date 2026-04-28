/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQDL } from '../../modules/nf-core/fastqdl/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: FETCH READS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FETCH_READS {

    take:
    ch_accessions // [ meta, accession ]

    main:

    ch_versions = Channel.empty()

    FASTQDL(
        ch_accessions
    )

    ch_versions = ch_versions.mix(FASTQDL.out.versions_fastqdl)

    ch_reads = FASTQDL.out.fastq
        .map { meta, reads ->
            def read_list = reads instanceof List ? reads : [reads]
            meta.single_end = (read_list.size() == 1)
            tuple(meta, read_list)
        }

    emit:
    reads       = ch_reads
    runinfo     = FASTQDL.out.runinfo
    runmergers  = FASTQDL.out.runmergers
    versions    = ch_versions
}
