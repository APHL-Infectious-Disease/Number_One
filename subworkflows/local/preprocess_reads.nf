/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_RAW }     from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../modules/nf-core/fastqc/main'
include { FASTP }                    from '../../modules/nf-core/fastp/main'
include { SEQKIT_STATS }             from '../../modules/nf-core/seqkit/stats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PREPROCESS READS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESS_READS {

    take:
    ch_reads // [ meta, [fastq1, fastq2?] ]

    main:

    ch_versions = Channel.empty()

    //
    // FASTQC on raw reads
    //
    ch_fastqc_raw_input = ch_reads
        .flatMap { meta, reads ->
            reads.collect { read -> tuple(meta, read) }
        }

    FASTQC_RAW(
        ch_fastqc_raw_input
    )

    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions_fastqc)

    //
    // FASTP input formatting
    // FASTP expects:
    // tuple val(meta), path(reads), path(adapter_fasta)
    // val discard_trimmed_pass
    // val save_trimmed_fail
    // val save_merged
    //
    ch_fastp_input = ch_reads.map { meta, reads ->
        tuple(meta, reads, [])
    }

    FASTP(
        ch_fastp_input,
        false,
        false,
        false
    )

    ch_versions = ch_versions.mix(FASTP.out.versions_fastp)

    ch_trimmed_reads = FASTP.out.reads

    //
    // FASTQC on trimmed reads
    //
    ch_fastqc_trimmed_input = ch_trimmed_reads
        .flatMap { meta, reads ->
            reads.collect { read -> tuple(meta, read) }
        }

    FASTQC_TRIMMED(
        ch_fastqc_trimmed_input
    )

    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions_fastqc)

    //
    // SeqKit stats on trimmed reads
    //
    ch_seqkit_input = ch_trimmed_reads
        .flatMap { meta, reads ->
            reads.collect { read -> tuple(meta, read) }
        }

    SEQKIT_STATS(
        ch_seqkit_input
    )

    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions_seqkit)

    emit:
    reads           = ch_trimmed_reads
    fastqc_raw      = FASTQC_RAW.out.zip
    fastqc_trimmed  = FASTQC_TRIMMED.out.zip
    fastp_html      = FASTP.out.html
    fastp_json      = FASTP.out.json
    seqkit_stats    = SEQKIT_STATS.out.stats
    versions        = ch_versions
}
