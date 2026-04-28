process MAKE_ACCESSIONS_METADATA {
    tag "make_accessions_metadata"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    path accession_file

    output:
    path "accessions_metadata.csv", emit: csv

    script:
    """
    python ${moduleDir}/make_accessions_metadata.py \
        --accessions ${accession_file} \
        --out accessions_metadata.csv
    """
}