#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import urllib.request
from pathlib import Path


BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output"
DEFAULT_OUT = Path("/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/opentargets_26_06")
DATASETS = [
    "disease",
    "target",
    "drug_molecule",
    "drug_mechanism_of_action",
    "clinical_indication",
    "association_by_datasource_direct",
]


def list_links(url: str) -> list[str]:
    with urllib.request.urlopen(url) as response:
        html = response.read().decode("utf-8", errors="replace")
    return re.findall(r'href="([^"]+)"', html)


def download_dataset(dataset: str, out_dir: Path, overwrite: bool) -> None:
    dataset_url = f"{BASE_URL}/{dataset}/"
    links = [link for link in list_links(dataset_url) if link.endswith(".parquet")]
    if not links:
        raise RuntimeError(f"No parquet files found at {dataset_url}")
    target_dir = out_dir / dataset
    target_dir.mkdir(parents=True, exist_ok=True)
    for link in links:
        url = dataset_url + link
        dest = target_dir / link
        if dest.exists() and not overwrite:
            print(f"skip existing {dest}")
            continue
        print(f"download {url} -> {dest}")
        urllib.request.urlretrieve(url, dest)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download selected Open Targets Platform parquet reference datasets.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--dataset", action="append", choices=DATASETS, help="Dataset to download; repeatable.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    datasets = args.dataset or DATASETS
    for dataset in datasets:
        download_dataset(dataset, args.out, args.overwrite)


if __name__ == "__main__":
    main()
