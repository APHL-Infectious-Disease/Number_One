/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ENTREZDIRECT_ESEARCH } from '../../modules/nf-core/entrezdirect/esearch'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBWORKFLOW: DISCOVER SRA RUNS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SRA_DISCOVERY {

    main:

    ch_versions = Channel.empty()

    ENTREZDIRECT_ESEARCH(
        tuple('group1_esearch', params.sra_query),
        "sra"
    )

    ch_versions = ch_versions.mix(ENTREZDIRECT_ESEARCH.out.versions_esearch)

    ch_runinfo_checked = ENTREZDIRECT_ESEARCH.out.runinfo
        .map { meta, csv ->
            if (!csv.exists() || csv.size() == 0) {
                error "ENTREZDIRECT_ESEARCH produced an empty runinfo CSV. Check params.sra_query: ${params.sra_query}"
            }
            tuple(meta, csv)
        }

    /*
     * Score candidate runs so the first max_runs selected are the most promising
     * for viruses likely to be in the Kraken database.
     */
    ch_accessions = ch_runinfo_checked
        .map { meta, csv -> csv }
        .splitCsv(header: true)
        .map { row ->
            def run       = row.Run?.toString()?.trim()
            def download  = row.download_path?.toString()?.trim()
            def strategy  = row.LibraryStrategy?.toString()?.trim()?.toUpperCase() ?: ''
            def source    = row.LibrarySource?.toString()?.trim()?.toUpperCase() ?: ''
            def sci       = row.ScientificName?.toString()?.toLowerCase() ?: ''
            def sample    = row.SampleName?.toString()?.toLowerCase() ?: ''
            def libname   = row.LibraryName?.toString()?.toLowerCase() ?: ''
            def text      = [sci, sample, libname].join(' ')

            if (!run || !run.startsWith('SRR') || !download) {
                return null
            }

            int score = 0

            // Prefer metagenomic libraries
            if (source == 'METAGENOMIC') score += 3

            // Prefer WGS or amplicon; both can be useful for your targets
            if (strategy == 'WGS') score += 2
            if (strategy == 'AMPLICON') score += 2

            // General viral/virome hints
            if (text.contains('virus')) score += 2
            if (text.contains('virome')) score += 3
            if (text.contains('viral')) score += 2

            // Specific target viruses / groups relevant to your DB
            if (text.contains('norovirus'))    score += 6
            if (text.contains('sapovirus'))    score += 6
            if (text.contains('astrovirus'))   score += 6
            if (text.contains('adenovirus'))   score += 6
            if (text.contains('rotavirus'))    score += 6
            if (text.contains('enterovirus'))  score += 6
            if (text.contains('poliovirus'))   score += 6
            if (text.contains('morbilli'))     score += 6
            if (text.contains('measles'))      score += 6
            if (text.contains('calicivirus'))  score += 4

            // Mild bonus for wastewater metagenome wording
            if (text.contains('wastewater metagenome')) score += 2
            if (text.contains('sewage')) score += 1

            // Penalize obviously non-viral libraries a bit
            if (source == 'GENOMIC' && !text.contains('virus') && !text.contains('virome')) score -= 3

            def meta = [
                id         : run,
                sample     : run,
                single_end : ((row.LibraryLayout ?: '').toString().toUpperCase() != 'PAIRED'),
                score      : score,
                strategy   : strategy,
                source     : source
            ]

            [
                score: score,
                run: run,
                meta: meta
            ]
        }
        .filter { it != null }
        .toSortedList { a, b ->
            b.score <=> a.score ?: a.run <=> b.run
        }
        .flatMap { rows ->
            rows.take(params.max_runs).collect { row ->
                tuple(row.meta, row.run)
            }
        }
        .ifEmpty {
            error "No suitable SRR accessions with download_path were parsed from runinfo CSV."
        }

    ch_meta_tsv = ch_runinfo_checked.map { meta, csv -> csv }

    emit:
    reads     = ch_accessions
    tsv       = ch_meta_tsv
    xml       = Channel.empty()
    versions  = ch_versions
}