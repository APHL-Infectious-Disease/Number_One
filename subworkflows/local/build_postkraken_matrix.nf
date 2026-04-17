/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { POSTKRAKEN_MATRIX } from '../../modules/local/postkraken_matrix/main'

workflow BUILD_POSTKRAKEN_MATRIX {

    take:
    ch_sra_meta_tsv
    ch_kraken_reports

    main:

    ch_versions = Channel.empty()

    // sra metadata is already a file channel; keep first file only
    ch_sra_meta_single = ch_sra_meta_tsv

    // kraken reports come in as tuple(meta, path) -> convert to plain paths
    ch_kraken_report_paths = ch_kraken_reports
        .map { item ->
            if (item instanceof List || item instanceof ArrayList) {
                return item[-1]
            }
            return item
        }
        .collect()

    POSTKRAKEN_MATRIX(
        ch_sra_meta_single,
        ch_kraken_report_paths
    )

    ch_versions = ch_versions.mix(POSTKRAKEN_MATRIX.out.versions)

    emit:
    matrix   = POSTKRAKEN_MATRIX.out.matrix
    versions = ch_versions
}