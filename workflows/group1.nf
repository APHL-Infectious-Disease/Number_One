/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SRA_DISCOVERY }          from '../subworkflows/local/sra_discovery'
include { FETCH_READS }            from '../subworkflows/local/fetch_reads'
include { PREPROCESS_READS }       from '../subworkflows/local/preprocess_reads'
include { CLASSIFY_READS }         from '../subworkflows/local/classify_reads'
include { BUILD_POSTKRAKEN_MATRIX } from '../subworkflows/local/build_postkraken_matrix'
include { REPORTING }              from '../subworkflows/local/reporting'

workflow GROUP1 {

    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    if (params.mode == 'sra') {
        SRA_DISCOVERY()
        ch_input = SRA_DISCOVERY.out.reads

        FETCH_READS(
            ch_input
        )
        ch_reads = FETCH_READS.out.reads

    } else if (params.mode == 'samplesheet') {
        ch_reads = ch_samplesheet

    } else {
        error "Unsupported params.mode: ${params.mode}. Use 'sra' or 'samplesheet'."
    }

    PREPROCESS_READS(
        ch_reads
    )

    ch_kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists: true)

    CLASSIFY_READS(
        PREPROCESS_READS.out.reads,
        ch_kraken2_db
    )

    BUILD_POSTKRAKEN_MATRIX(
        params.mode == 'sra' ? SRA_DISCOVERY.out.tsv : Channel.empty(),
        CLASSIFY_READS.out.kraken2_report
    )

    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'group1_software_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    ch_multiqc_files = Channel.empty()
        .mix(PREPROCESS_READS.out.fastqc_raw)
        .mix(PREPROCESS_READS.out.fastqc_trimmed)
        .mix(PREPROCESS_READS.out.fastp_json)
        .mix(PREPROCESS_READS.out.seqkit_stats)
        .mix(CLASSIFY_READS.out.kraken2_report)

    REPORTING(
        ch_multiqc_files,
        ch_collated_versions
    )

    emit:
    reads                = PREPROCESS_READS.out.reads
    kraken2_report       = CLASSIFY_READS.out.kraken2_report
    metadata_postkraken  = BUILD_POSTKRAKEN_MATRIX.out.matrix
    multiqc_report       = REPORTING.out.report
    versions             = ch_versions
}
