process ENTREZDIRECT_ESEARCH {
    tag "$meta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(term)
    val database

    output:
    tuple val(meta), path("*.csv"), emit: runinfo
    tuple val(meta), path("*.stderr.txt"), emit: stderr
    tuple val("${task.process}"), val('PYTHON_REQUESTS'), val('1'), emit: versions_esearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta}"

    """
    python ${moduleDir}/fetch_sra_runinfo.py \
        --db "${database}" \
        --query "${term}" \
        --out "${prefix}.csv" \
        2> "${prefix}.stderr.txt"

    test -s "${prefix}.csv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    touch "${prefix}.csv"
    touch "${prefix}.stderr.txt"
    """
}