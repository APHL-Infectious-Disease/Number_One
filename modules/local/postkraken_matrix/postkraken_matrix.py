#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import re

COLUMN_MAP = {
    "Norovirus": ["norovirus"],
    "Rotavirus": ["rotavirus"],
    "Parvovirus": ["parvovirus"],
    "Enterovirus": ["enterovirus"],
    "Astrovirus": ["astrovirus"],
    "Poliovirus": ["poliovirus"],
    "Sapovirus": ["sapovirus"],
    "Adenovirus": ["adenovirus"],
    "Morbillivirus": ["morbillivirus", "measles virus"],
    "Unclassified": ["unclassified"],
}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sra-meta", required=True)
    parser.add_argument("--kraken-reports", nargs="+", required=True)
    return parser.parse_args()

def infer_sra_id(filepath: str) -> str:
    name = Path(filepath).name
    patterns = [
        r"(.+?)\.kraken2\.report\.txt$",
        r"(.+?)\.report\.txt$",
        r"(.+?)\.txt$",
    ]
    for pattern in patterns:
        m = re.match(pattern, name)
        if m:
            return m.group(1)
    return Path(filepath).stem

def normalize(text: str) -> str:
    return str(text).strip().lower()

def parse_kraken_report(report_file: str):
    sra_id = infer_sra_id(report_file)

    row = {
        "SRA accession": sra_id,
        "Norovirus": "0",
        "Rotavirus": "0",
        "Parvovirus": "0",
        "Enterovirus": "0",
        "Astrovirus": "0",
        "Poliovirus": "0",
        "Sapovirus": "0",
        "Adenovirus": "0",
        "Morbillivirus": "0",
        "Unclassified": "0",
    }

    with open(report_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            try:
                percent = float(parts[0].strip())
            except ValueError:
                continue

            taxid = parts[4].strip()
            organism = parts[5].strip()
            org_norm = normalize(organism)

            for col, keywords in COLUMN_MAP.items():
                if any(k in org_norm for k in keywords):
                    row[col] = f"{percent};{taxid};{organism}"
                    break

    return row

def load_metadata(meta_file: str) -> pd.DataFrame:
    # SRA runinfo is CSV, not TSV
    df = pd.read_csv(meta_file, dtype=str)

    # Map known SRA runinfo columns into your team's output columns
    rename_map = {}

    if "Run" in df.columns:
        rename_map["Run"] = "SRA accession"
    elif "run_accession" in df.columns:
        rename_map["run_accession"] = "SRA accession"

    # source
    if "ScientificName" in df.columns:
        rename_map["ScientificName"] = "source"
    elif "scientific_name" in df.columns:
        rename_map["scientific_name"] = "source"

    # collection date
    if "Collection_Date" in df.columns:
        rename_map["Collection_Date"] = "collection date"
    elif "collection_date" in df.columns:
        rename_map["collection_date"] = "collection date"

    # location
    if "geo_loc_name_country_calc" in df.columns:
        rename_map["geo_loc_name_country_calc"] = "location"
    elif "geo_loc_name" in df.columns:
        rename_map["geo_loc_name"] = "location"
    elif "location" in df.columns:
        rename_map["location"] = "location"

    df = df.rename(columns=rename_map)

    if "SRA accession" not in df.columns:
        raise ValueError("Could not find a Run/run_accession column in runinfo CSV")

    for col in ["source", "collection date", "location"]:
        if col not in df.columns:
            df[col] = ""

    return df[["SRA accession", "source", "collection date", "location"]].copy()

def main():
    args = parse_args()

    metadata_df = load_metadata(args.sra_meta)
    kraken_rows = [parse_kraken_report(fp) for fp in args.kraken_reports]
    kraken_df = pd.DataFrame(kraken_rows)

    final_df = metadata_df.merge(kraken_df, how="left", on="SRA accession")

    for col in [
        "Norovirus", "Rotavirus", "Parvovirus", "Enterovirus", "Astrovirus",
        "Poliovirus", "Sapovirus", "Adenovirus", "Morbillivirus", "Unclassified"
    ]:
        final_df[col] = final_df[col].fillna("0")

    final_df = final_df.sort_values("SRA accession")
    final_df.to_csv("metadata_postkraken.csv", index=False)

if __name__ == "__main__":
    main()
