process PREPARE_KRAKEN_DB_UNPACK {
    tag "$name"
    label 'process_single'

    input:
    tuple val(name), path(archive)

    output:
    path("kraken_db"), emit: db

    script:
    """
    mkdir -p kraken_db

    tar -xzf "${archive}" -C kraken_db

    if [ ! -f kraken_db/hash.k2d ]; then
        subdir=\$(find kraken_db -mindepth 1 -maxdepth 1 -type d | head -n 1)
        if [ -n "\$subdir" ] && [ -f "\$subdir/hash.k2d" ]; then
            mv "\$subdir"/* kraken_db/
            rmdir "\$subdir"
        fi
    fi

    test -f kraken_db/hash.k2d
    test -f kraken_db/opts.k2d
    test -f kraken_db/taxo.k2d
    """

    stub:
    """
    mkdir -p kraken_db
    touch kraken_db/hash.k2d
    touch kraken_db/opts.k2d
    touch kraken_db/taxo.k2d
    """
}