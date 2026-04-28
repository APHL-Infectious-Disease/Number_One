/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_KRAKEN_DB_UNPACK } from '../../modules/local/prepare_kraken_db_unpack/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE KRAKEN2 DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_KRAKEN_DB {

    main:

    if (params.kraken2_db_url) {
        ch_db_archive = Channel
            .fromPath(params.kraken2_db_url, checkIfExists: true)
            .ifEmpty {
                error "Could not resolve Kraken2 DB URL: ${params.kraken2_db_url}"
            }

        ch_unpack_input = ch_db_archive.map { archive ->
            tuple('kraken_db', archive)
        }

        PREPARE_KRAKEN_DB_UNPACK(
            ch_unpack_input
        )

        ch_db = PREPARE_KRAKEN_DB_UNPACK.out.db.first()

    } else if (params.kraken2_db) {
        ch_db = Channel
            .fromPath(params.kraken2_db, checkIfExists: true)
            .first()

    } else {
        error "Provide either --kraken2_db or --kraken2_db_url"
    }

    emit:
    db = ch_db
}
