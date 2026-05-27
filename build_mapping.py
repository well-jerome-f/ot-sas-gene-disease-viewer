#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import pandas as pd


BASE_DIR = Path("/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet")

DEFAULT_VARIANTS = BASE_DIR / "gnomad.sas.enriched.investigate.parquet"
DEFAULT_OT_EVIDENCE = BASE_DIR / "opentargets_disease_target_credset_filtered.parquet"
DEFAULT_DISEASE = BASE_DIR / "disease" / "disease.parquet"
DEFAULT_OUTPUT = BASE_DIR / "ot_sas_disease_gene_variant_map.score_ge_0.25.parquet"

NON_DISEASE_THERAPEUTIC_AREAS = {
    "EFO_0000651",  # phenotype
    "EFO_0001444",  # measurement
    "EFO_0002571",  # medical procedure
    "EFO_0005932",  # animal disease
    "GO_0008150",  # biological_process
}

CHROM_ORDER = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}


def listify(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Iterable):
        return [str(v) for v in value]
    return []


def is_disease_trait(therapeutic_areas: object, excluded: set[str]) -> bool:
    areas = listify(therapeutic_areas)
    return bool(areas) and any(area not in excluded for area in areas)


def chrom_label(chrom: object) -> str:
    value = str(chrom)
    return value[3:] if value.lower().startswith("chr") else value


def chrom_sort_key(chrom: object) -> int:
    return CHROM_ORDER.get(chrom_label(chrom), 99)


def first_non_null(values: pd.Series) -> object:
    values = values.dropna()
    return values.iloc[0] if len(values) else None


def semicolon_join(values: pd.Series) -> str:
    seen: list[str] = []
    for value in values.dropna():
        for item in listify(value):
            item = str(item)
            if item and item not in seen:
                seen.append(item)
    return ";".join(seen)


def load_disease_terms(path: Path, excluded_areas: set[str]) -> pd.DataFrame:
    disease = pd.read_parquet(
        path,
        columns=["id", "name", "description", "therapeuticAreas"],
    )
    disease = disease.rename(
        columns={
            "id": "disease_id",
            "name": "disease_name_from_metadata",
            "description": "disease_description_from_metadata",
            "therapeuticAreas": "therapeutic_areas",
        }
    )
    disease["is_disease_trait"] = disease["therapeutic_areas"].map(
        lambda areas: is_disease_trait(areas, excluded_areas)
    )
    disease["therapeutic_areas"] = disease["therapeutic_areas"].map(lambda areas: ";".join(listify(areas)))
    return disease


def load_variants(path: Path) -> pd.DataFrame:
    columns = [
        "CHROM",
        "POS",
        "varid",
        "AF_sas",
        "AF_nfe",
        "sas_enrichment",
        "consequence",
        "impact",
        "lof",
        "symbol",
        "gene",
    ]
    variants = pd.read_parquet(path, columns=columns)
    variants = variants.rename(
        columns={
            "CHROM": "chrom",
            "POS": "pos",
            "gene": "ensembl_gene_id",
            "symbol": "gene_symbol",
        }
    )
    variants["chrom"] = variants["chrom"].map(chrom_label)
    variants["chrom_sort"] = variants["chrom"].map(chrom_sort_key)
    variants["variant_is_lof"] = variants["lof"].fillna("").astype(str).str.strip().ne("") | variants[
        "impact"
    ].eq("HIGH")
    variants["gene_has_lof"] = variants.groupby("ensembl_gene_id")["variant_is_lof"].transform("any")
    variants = variants.drop_duplicates(
        [
            "ensembl_gene_id",
            "chrom",
            "pos",
            "varid",
            "AF_sas",
            "AF_nfe",
            "sas_enrichment",
            "consequence",
            "variant_is_lof",
        ]
    )
    return variants.sort_values(["ensembl_gene_id", "chrom_sort", "pos", "varid"], kind="mergesort")


def load_and_aggregate_evidence(path: Path, disease_terms: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "targetId",
        "diseaseId",
        "diseaseFromSourceMappedId",
        "score",
        "resourceScore",
        "studyLocusId",
        "literature",
        "publicationDate",
        "evidenceDate",
        "name",
        "description",
    ]
    evidence = pd.read_parquet(path, columns=columns)
    evidence = evidence.rename(
        columns={
            "targetId": "ensembl_gene_id",
            "diseaseId": "disease_id",
            "diseaseFromSourceMappedId": "disease_from_source_mapped_id",
            "score": "ot_score",
            "resourceScore": "resource_score",
            "name": "disease_name",
            "description": "disease_description",
        }
    )
    evidence = evidence.merge(disease_terms, on="disease_id", how="left")
    evidence = evidence[evidence["is_disease_trait"].fillna(False)].copy()
    evidence["disease_name"] = evidence["disease_name"].fillna(evidence["disease_name_from_metadata"])
    evidence["disease_description"] = evidence["disease_description"].fillna(
        evidence["disease_description_from_metadata"]
    )

    best_idx = evidence.groupby(["ensembl_gene_id", "disease_id"], dropna=False)["ot_score"].idxmax()
    best = evidence.loc[
        best_idx,
        [
            "ensembl_gene_id",
            "disease_id",
            "disease_from_source_mapped_id",
            "disease_name",
            "disease_description",
            "ot_score",
            "resource_score",
            "studyLocusId",
            "publicationDate",
            "evidenceDate",
            "therapeutic_areas",
        ],
    ].rename(
        columns={
            "studyLocusId": "best_study_locus_id",
            "publicationDate": "best_publication_date",
            "evidenceDate": "best_evidence_date",
        }
    )

    counts = (
        evidence.groupby(["ensembl_gene_id", "disease_id"], dropna=False)
        .agg(
            evidence_count=("ot_score", "size"),
            study_locus_count=("studyLocusId", "nunique"),
            literature_pmids=("literature", semicolon_join),
        )
        .reset_index()
    )
    return best.merge(counts, on=["ensembl_gene_id", "disease_id"], how="left")


def build_mapping(
    variants_path: Path,
    evidence_path: Path,
    disease_path: Path,
    output_path: Path,
    excluded_areas: set[str],
    tsv_output_path: Path | None,
    min_score: float,
) -> pd.DataFrame:
    disease_terms = load_disease_terms(disease_path, excluded_areas)
    variants = load_variants(variants_path)
    disease_gene = load_and_aggregate_evidence(evidence_path, disease_terms)
    disease_gene = disease_gene[disease_gene["ot_score"] >= min_score].copy()

    mapped = disease_gene.merge(variants, on="ensembl_gene_id", how="inner")
    mapped = mapped[
        [
            "ensembl_gene_id",
            "gene_symbol",
            "disease_id",
            "disease_from_source_mapped_id",
            "disease_name",
            "disease_description",
            "ot_score",
            "resource_score",
            "evidence_count",
            "study_locus_count",
            "best_study_locus_id",
            "best_publication_date",
            "best_evidence_date",
            "literature_pmids",
            "therapeutic_areas",
            "chrom",
            "chrom_sort",
            "pos",
            "varid",
            "AF_sas",
            "AF_nfe",
            "sas_enrichment",
            "consequence",
            "impact",
            "lof",
            "variant_is_lof",
            "gene_has_lof",
        ]
    ].sort_values(
        ["gene_symbol", "ensembl_gene_id", "ot_score", "disease_name", "chrom_sort", "pos"],
        ascending=[True, True, False, True, True, True],
        kind="mergesort",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapped.to_parquet(output_path, index=False)
    stats = {
        "min_l2g_score": min_score,
        "mapped_rows": int(len(mapped)),
        "genes": int(mapped["ensembl_gene_id"].nunique()),
        "lof_positive_genes": int(mapped.loc[mapped["gene_has_lof"], "ensembl_gene_id"].nunique()),
        "disease_traits": int(mapped["disease_id"].nunique()),
        "disease_gene_pairs": int(mapped[["ensembl_gene_id", "disease_id"]].drop_duplicates().shape[0]),
        "variants": int(mapped["varid"].nunique()),
        "lof_variants": int(mapped.loc[mapped["variant_is_lof"], "varid"].nunique()),
        "median_l2g_score_flat_rows": float(mapped["ot_score"].median()) if len(mapped) else None,
        "median_l2g_score_disease_gene_pairs": float(
            mapped.drop_duplicates(["ensembl_gene_id", "disease_id"])["ot_score"].median()
        )
        if len(mapped)
        else None,
    }
    output_path.with_suffix(output_path.suffix + ".stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    if tsv_output_path:
        mapped.to_csv(tsv_output_path, sep="\t", index=False)
    return mapped


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a disease-gene-variant mapping for SAS-enriched coding variants and Open Targets GWAS credible-set evidence."
    )
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_OT_EVIDENCE)
    parser.add_argument("--disease", type=Path, default=DEFAULT_DISEASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--tsv-output", type=Path, default=None)
    parser.add_argument("--min-score", type=float, default=0.25, help="Minimum Open Targets/L2G score to retain.")
    parser.add_argument(
        "--exclude-therapeutic-area",
        action="append",
        default=sorted(NON_DISEASE_THERAPEUTIC_AREAS),
        help="Therapeutic area ID to exclude. Repeat to add multiple IDs.",
    )
    args = parser.parse_args()

    mapped = build_mapping(
        variants_path=args.variants,
        evidence_path=args.evidence,
        disease_path=args.disease,
        output_path=args.output,
        excluded_areas=set(args.exclude_therapeutic_area),
        tsv_output_path=args.tsv_output,
        min_score=args.min_score,
    )
    print(f"Wrote {len(mapped):,} disease-gene-variant rows to {args.output}")
    print(f"Genes: {mapped['ensembl_gene_id'].nunique():,}")
    print(f"Diseases: {mapped['disease_id'].nunique():,}")
    print(f"Variants: {mapped['varid'].nunique():,}")
    if args.tsv_output:
        print(f"Wrote TSV copy to {args.tsv_output}")


if __name__ == "__main__":
    main()
