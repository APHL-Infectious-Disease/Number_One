#!/usr/bin/env python3

import argparse
import csv
import sys
import time
import requests

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--query", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    session = requests.Session()
    session.headers.update({
        "User-Agent": "APHL-Infectious-Disease-group1/1.0"
    })

    esearch_params = {
        "db": args.db,
        "term": args.query,
        "retmax": 1000,
        "usehistory": "y",
        "retmode": "json",
    }

    resp = session.get(ESEARCH_URL, params=esearch_params, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    esearch = data.get("esearchresult", {})
    count = int(esearch.get("count", "0"))
    webenv = esearch.get("webenv", "")
    query_key = esearch.get("querykey", "")

    if count == 0 or not webenv or not query_key:
        raise RuntimeError(f"No SRA results found for query: {args.query}")

    efetch_params = {
        "db": args.db,
        "query_key": query_key,
        "WebEnv": webenv,
        "rettype": "runinfo",
        "retmode": "text",
    }

    fetch = session.get(EFETCH_URL, params=efetch_params, timeout=120)
    fetch.raise_for_status()
    text = fetch.text.strip()

    if not text:
        raise RuntimeError("efetch returned empty runinfo text")

    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        if not text.endswith("\n"):
            handle.write("\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
