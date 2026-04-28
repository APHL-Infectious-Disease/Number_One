#!/usr/bin/env bash
set -euo pipefail

INPUT_XML="${1:-input.xml}"

if [[ ! -f "${INPUT_XML}" ]]; then
    echo "ERROR: Input XML not found: ${INPUT_XML}" >&2
    exit 1
fi

cp "${INPUT_XML}" sra.xml

cat > xtract_dynamic.cmd <<'EOF'
-pattern EXPERIMENT_PACKAGE
-first PRIMARY_ID
-element LIBRARY_STRATEGY
-element LIBRARY_SELECTION
-element TAXON_ID
-element SCIENTIFIC_NAME
-block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE
-if TAG -equals collection_date -element VALUE
-block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE
-if TAG -equals geo_loc_name -element VALUE
-block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE
-if TAG -equals lat_lon -element VALUE
EOF

cat > sra_meta.tsv <<'EOF'
run_accession	library_strategy	library_selection	taxon_id	scientific_name	collection_date	location	lat_lon
EOF

xtract -input sra.xml $(cat xtract_dynamic.cmd) > extracted.tsv || true

echo "DEBUG: extracted.tsv line count = $(wc -l < extracted.tsv 2>/dev/null || echo 0)" >&2
echo "DEBUG: first 5 lines of extracted.tsv:" >&2
head -5 extracted.tsv >&2 || true

if [[ -s extracted.tsv ]]; then
    grep "USA" extracted.tsv >> sra_meta.tsv || true
fi

if [[ $(wc -l < sra_meta.tsv) -le 1 ]] && [[ -s extracted.tsv ]]; then
    cat extracted.tsv >> sra_meta.tsv
fi

test -f sra_meta.tsv
