/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC } from '../../modules/nf-core/multiqc/main'

workflow REPORTING {

    take:
    ch_multiqc_files
    ch_versions_yml

    main:
    ch_versions = Channel.empty()

    // Convert incoming tuple(meta, path) or path into plain path objects
    ch_multiqc_paths = ch_multiqc_files.map { item ->
        if (item instanceof List || item instanceof ArrayList) {
            return item[-1]
        }
        return item
    }

    // Add the software versions YAML directly as another MultiQC input file
    ch_all_multiqc_files = ch_multiqc_paths.mix(ch_versions_yml)

    // MULTIQC module expects:
    // tuple val(meta), path(multiqc_files), path(multiqc_config), path(multiqc_logo), path(replace_names), path(sample_names)
    ch_multiqc_tuple = ch_all_multiqc_files
        .collect()
        .map { files ->
            tuple(
                [id: 'multiqc'],
                files,
                [],
                [],
                [],
                []
            )
        }

    MULTIQC(
        ch_multiqc_tuple
    )

    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    report   = MULTIQC.out.report
    data     = MULTIQC.out.data
    versions = ch_versions
}