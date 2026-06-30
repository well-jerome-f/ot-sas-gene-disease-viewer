#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd
import polars as pl


APP_DIR = Path(__file__).resolve().parent
BASE_DIR = Path("/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet")

DEFAULT_VEP = BASE_DIR / "vep.annotated.variant_gene_consequence.parquet"
DEFAULT_GNOMAD = BASE_DIR / "gnomad.india_or_sas.af_ge_0.001.enriched_vs_nfe_fin.parquet"
DEFAULT_GENOME_INDIA = BASE_DIR / "genome.india.af_ge_0.001.not_in_gnomad.parquet"
DEFAULT_EVIDENCE = BASE_DIR / "opentargets_disease_target_credset_filtered.parquet"
DEFAULT_DISEASE = BASE_DIR / "disease" / "disease.parquet"
DEFAULT_UNMET_NEEDS = APP_DIR / "data" / "unmet-needs-index.csv"

DEFAULT_ANNOTATED = BASE_DIR / "vep115.coding.with_sas_india_enrichment.parquet"
DEFAULT_HIGH_MODERATE = BASE_DIR / "vep115.plof_deleterious_missense.gene_variants.with_global_maf.parquet"
DEFAULT_EUR_GENES = BASE_DIR / "vep115.plof_deleterious_missense.eur_rare_excluded_genes.parquet"
DEFAULT_OUTPUT = BASE_DIR / "ot_sas_gi_vep115_plof_deleterious_missense_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet"

CHROM_ORDER = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}

LOF_CONSEQUENCE_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
}

ALPHAMISSENSE_DELETERIOUS_CLASSES = {
    "likely_pathogenic",
    "likely pathogenic",
    "pathogenic",
}

ESM1B_CANDIDATE_COLUMNS = [
    "ESM1b",
    "ESM1b_score",
    "ESM1b_rankscore",
    "esm1b",
    "esm1b_score",
    "esm1b_rankscore",
]


def import_build_mapping_helpers():
    sys.path.insert(0, str(APP_DIR))
    from build_mapping import (  # type: ignore
        NON_DISEASE_THERAPEUTIC_AREAS,
        load_and_aggregate_evidence,
        load_disease_terms,
        load_unmet_needs,
    )

    return NON_DISEASE_THERAPEUTIC_AREAS, load_and_aggregate_evidence, load_disease_terms, load_unmet_needs


def chrom_label(chrom: object) -> str:
    value = str(chrom)
    return value[3:] if value.lower().startswith("chr") else value


def chrom_sort_key(chrom: object) -> int:
    return CHROM_ORDER.get(chrom_label(chrom), 99)


def cumulative_frequency(values: pd.Series) -> float:
    freqs = pd.to_numeric(values, errors="coerce").dropna().clip(lower=0.0, upper=1.0)
    if freqs.empty:
        return 0.0
    if (freqs >= 1.0).any():
        return 1.0
    return 1.0 - math.exp(float((1.0 - freqs).map(math.log).sum()))


def optional_numeric_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            return col
    return None


def parse_float_expr(name: str) -> pl.Expr:
    if name not in pl.scan_parquet(DEFAULT_VEP).collect_schema():
        return pl.lit(None, dtype=pl.Float64)
    return (
        pl.col(name)
        .cast(pl.Utf8, strict=False)
        .replace({"": None, ".": None, "invalid_field": None})
        .cast(pl.Float64, strict=False)
    )


def annotate_vep_with_frequency(
    vep_path: Path,
    gnomad_path: Path,
    genome_india_path: Path,
    output_path: Path,
) -> pl.DataFrame:
    vep_schema = pl.scan_parquet(vep_path).collect_schema()

    def maybe_float(name: str) -> pl.Expr:
        if name not in vep_schema:
            return pl.lit(None, dtype=pl.Float64)
        return (
            pl.col(name)
            .cast(pl.Utf8, strict=False)
            .replace({"": None, ".": None, "invalid_field": None})
            .cast(pl.Float64, strict=False)
        )

    gnomad = (
        pl.scan_parquet(gnomad_path)
        .select(
            "varid",
            pl.col("AF_sas").cast(pl.Float64).alias("AF_sas_source"),
            pl.col("AF_nfe").cast(pl.Float64).alias("AF_nfe_source"),
            pl.col("AF_fin").cast(pl.Float64).alias("AF_fin_source"),
            pl.col("genome_india_af").cast(pl.Float64).alias("genome_india_af_gnomad_match"),
            pl.col("sas_vs_nfe_enrichment").cast(pl.Float64),
            pl.col("sas_vs_fin_enrichment").cast(pl.Float64),
            pl.col("genome_india_vs_nfe_enrichment").cast(pl.Float64),
            pl.col("genome_india_vs_fin_enrichment").cast(pl.Float64),
            pl.col("genome_india_match_type"),
            pl.col("genome_india_varid"),
        )
        .unique("varid")
    )
    genome_india = (
        pl.scan_parquet(genome_india_path)
        .select(
            "varid",
            pl.col("AF").cast(pl.Float64).alias("genome_india_af_private"),
            pl.col("present_gnomad"),
        )
        .unique("varid")
    )

    annotated = (
        pl.scan_parquet(vep_path)
        .join(gnomad, on="varid", how="left")
        .join(genome_india, on="varid", how="left")
        .with_columns(
            [
                pl.coalesce([pl.col("AF_sas_source"), maybe_float("gnomAD4_1_joint_SAS_AF")]).alias("AF_sas"),
                pl.coalesce([pl.col("AF_nfe_source"), maybe_float("gnomAD4_1_joint_NFE_AF")]).alias("AF_nfe"),
                pl.coalesce([pl.col("AF_fin_source"), maybe_float("gnomAD4_1_joint_FIN_AF")]).alias("AF_fin"),
                pl.coalesce([pl.col("genome_india_af_gnomad_match"), pl.col("genome_india_af_private")]).alias(
                    "genome_india_af"
                ),
                maybe_float("gnomAD4_1_joint_AF").alias("gnomad_joint_af"),
                maybe_float("gnomAD4_1_joint_POPMAX_AF").alias("gnomad_joint_popmax_af"),
                maybe_float("gnomAD4_1_joint_SAS_AF").alias("gnomad_joint_sas_af"),
                maybe_float("AllOfUs_gvs_sas_af").alias("allofus_sas_af"),
                maybe_float("AllOfUs_gvs_eur_af").alias("allofus_eur_af"),
                pl.col("CADD_PHRED").cast(pl.Float64, strict=False).alias("cadd_phred"),
                pl.col("CADD_RAW").cast(pl.Float64, strict=False).alias("cadd_raw"),
            ]
        )
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    annotated.collect().write_parquet(output_path, compression="zstd")
    return pl.read_parquet(output_path)


def is_lof_consequence(consequence: object, impact: object) -> bool:
    if str(impact) == "HIGH":
        return True
    terms = {term.strip() for term in str(consequence or "").split("&")}
    return bool(terms & LOF_CONSEQUENCE_TERMS)


def prepare_high_moderate_variants(
    annotated_path: Path,
    output_path: Path,
    max_plof_sas_af: float,
    cadd_missense_threshold: float,
    alphamissense_threshold: float,
    esm1b_threshold: float | None,
) -> pd.DataFrame:
    df = pd.read_parquet(annotated_path)
    df = df[df["IMPACT"].isin(["HIGH", "MODERATE"])].copy()
    df = df[df["gene_id"].notna() & df["varid"].notna()].copy()

    df["chrom"] = df["CHROM"].map(chrom_label)
    df["chrom_sort"] = df["chrom"].map(chrom_sort_key)
    df["pos"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    df["ensembl_gene_id"] = df["gene_id"]
    df["gene_symbol"] = df["SYMBOL"].fillna("")
    df["consequence"] = df["Consequence"]
    df["impact"] = df["IMPACT"]
    df["lof"] = [
        consequence if is_lof_consequence(consequence, impact) else ""
        for consequence, impact in zip(df["consequence"], df["impact"])
    ]
    df["variant_is_lof"] = df["lof"].astype(str).str.len().gt(0)
    df["variant_is_missense"] = df["consequence"].fillna("").astype(str).str.contains("missense_variant", regex=False)

    for col in [
        "AF_sas",
        "AF_nfe",
        "AF_fin",
        "genome_india_af",
        "cadd_phred",
        "sas_vs_nfe_enrichment",
        "sas_vs_fin_enrichment",
        "genome_india_vs_nfe_enrichment",
        "genome_india_vs_fin_enrichment",
        "gnomad_joint_sas_af",
        "allofus_sas_af",
        "allofus_eur_af",
        "LOEUF",
        "am_pathogenicity",
        "CADD_PHRED",
        "CADD_RAW",
        "REVEL",
        "ClinPred",
        "EVE_SCORE",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    esm1b_col = optional_numeric_column(df, ESM1B_CANDIDATE_COLUMNS)

    df["plof_is_not_common_sas"] = df["variant_is_lof"] & df["AF_sas"].fillna(0.0).le(max_plof_sas_af)
    alpha_class = df.get("am_class", pd.Series(index=df.index, dtype=object)).fillna("").astype(str).str.lower()
    alpha_score = df.get("am_pathogenicity", pd.Series(index=df.index, dtype=float))
    df["missense_deleterious_by_cadd"] = df["variant_is_missense"] & df["cadd_phred"].ge(cadd_missense_threshold)
    df["missense_deleterious_by_alphamissense"] = df["variant_is_missense"] & (
        alpha_class.isin(ALPHAMISSENSE_DELETERIOUS_CLASSES) | pd.to_numeric(alpha_score, errors="coerce").ge(alphamissense_threshold)
    )
    if esm1b_col and esm1b_threshold is not None:
        df["ESM1b_score"] = df[esm1b_col]
        df["missense_deleterious_by_esm1b"] = df["variant_is_missense"] & df["ESM1b_score"].ge(esm1b_threshold)
    else:
        df["ESM1b_score"] = pd.NA
        df["missense_deleterious_by_esm1b"] = False

    df["variant_is_deleterious_missense"] = (
        df["missense_deleterious_by_cadd"]
        | df["missense_deleterious_by_alphamissense"]
        | df["missense_deleterious_by_esm1b"]
    )
    df["variant_selected_for_global_maf"] = df["plof_is_not_common_sas"] | df["variant_is_deleterious_missense"]
    df = df[df["variant_selected_for_global_maf"]].copy()

    df["sas_enrichment"] = df[["sas_vs_nfe_enrichment", "sas_vs_fin_enrichment"]].max(axis=1, skipna=True)
    df["genome_india_enrichment"] = df[
        ["genome_india_vs_nfe_enrichment", "genome_india_vs_fin_enrichment"]
    ].max(axis=1, skipna=True)
    df["gene_has_lof"] = df.groupby("ensembl_gene_id")["plof_is_not_common_sas"].transform("any")
    df["gene_has_deleterious_missense"] = df.groupby("ensembl_gene_id")["variant_is_deleterious_missense"].transform("any")

    unique_variant_classes = df.drop_duplicates(
        [
            "ensembl_gene_id",
            "varid",
            "plof_is_not_common_sas",
            "variant_is_deleterious_missense",
            "variant_selected_for_global_maf",
        ]
    )
    for source_col, prefix in [("AF_sas", "cum_af_sas"), ("genome_india_af", "cum_af_genome_india")]:
        lof_caf = (
            unique_variant_classes[unique_variant_classes["plof_is_not_common_sas"]]
            .groupby("ensembl_gene_id")[source_col]
            .apply(cumulative_frequency)
            .rename(f"{prefix}_lof")
        )
        missense_caf = (
            unique_variant_classes[unique_variant_classes["variant_is_deleterious_missense"]]
            .groupby("ensembl_gene_id")[source_col]
            .apply(cumulative_frequency)
            .rename(f"{prefix}_missense")
        )
        global_maf = (
            unique_variant_classes[unique_variant_classes["variant_selected_for_global_maf"]]
            .groupby("ensembl_gene_id")[source_col]
            .apply(cumulative_frequency)
            .rename(f"global_maf_{source_col.lower()}")
        )
        df = df.join(lof_caf, on="ensembl_gene_id").join(missense_caf, on="ensembl_gene_id")
        df = df.join(global_maf, on="ensembl_gene_id")
        df[f"{prefix}_lof"] = df[f"{prefix}_lof"].fillna(0.0)
        df[f"{prefix}_missense"] = df[f"{prefix}_missense"].fillna(0.0)
        df[f"global_maf_{source_col.lower()}"] = df[f"global_maf_{source_col.lower()}"].fillna(0.0)

    keep = [
        "ensembl_gene_id",
        "gene_symbol",
        "chrom",
        "chrom_sort",
        "pos",
        "varid",
        "REF",
        "ALT",
        "Allele",
        "Feature",
        "BIOTYPE",
        "consequence",
        "impact",
        "lof",
        "variant_is_lof",
        "variant_is_missense",
        "plof_is_not_common_sas",
        "variant_is_deleterious_missense",
        "variant_selected_for_global_maf",
        "missense_deleterious_by_cadd",
        "missense_deleterious_by_alphamissense",
        "missense_deleterious_by_esm1b",
        "gene_has_lof",
        "gene_has_deleterious_missense",
        "AF_sas",
        "AF_nfe",
        "AF_fin",
        "genome_india_af",
        "sas_enrichment",
        "genome_india_enrichment",
        "sas_vs_nfe_enrichment",
        "sas_vs_fin_enrichment",
        "genome_india_vs_nfe_enrichment",
        "genome_india_vs_fin_enrichment",
        "gnomad_joint_sas_af",
        "allofus_sas_af",
        "allofus_eur_af",
        "cadd_phred",
        "CADD_RAW",
        "REVEL",
        "ESM1b_score",
        "SIFT",
        "PolyPhen",
        "am_class",
        "am_pathogenicity",
        "LOEUF",
        "ClinPred",
        "EVE_CLASS",
        "EVE_SCORE",
        "cum_af_sas_lof",
        "cum_af_sas_missense",
        "cum_af_genome_india_lof",
        "cum_af_genome_india_missense",
        "global_maf_af_sas",
        "global_maf_genome_india_af",
    ]
    keep = [col for col in keep if col in df.columns]
    df = df[keep].drop_duplicates()
    df = df.sort_values(["ensembl_gene_id", "chrom_sort", "pos", "varid"], kind="mergesort")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False)
    return df


def write_eur_excluded_genes(variants: pd.DataFrame, output_path: Path, eur_rare_af_max: float) -> set[str]:
    eur_af = variants[["AF_nfe", "AF_fin"]].apply(pd.to_numeric, errors="coerce").fillna(0.0).max(axis=1)
    eur_rare_mask = eur_af.gt(0) & eur_af.le(eur_rare_af_max)
    excluded = (
        variants.loc[eur_rare_mask, ["ensembl_gene_id", "gene_symbol"]]
        .drop_duplicates()
        .sort_values(["gene_symbol", "ensembl_gene_id"], kind="mergesort")
    )
    excluded.to_parquet(output_path, index=False)
    return set(excluded["ensembl_gene_id"].dropna())


def semicolon_join(values: pd.Series) -> str:
    seen: list[str] = []
    for value in values.dropna():
        for item in str(value).split(";"):
            if item and item not in seen:
                seen.append(item)
    return ";".join(seen)


def build_disease_mapping(
    variants: pd.DataFrame,
    evidence_path: Path,
    disease_path: Path,
    unmet_needs_path: Path,
    output_path: Path,
    eur_genes: set[str],
    min_score: float,
    eur_rare_af_max: float,
    max_plof_sas_af: float,
    cadd_missense_threshold: float,
    alphamissense_threshold: float,
    esm1b_threshold: float | None,
) -> pd.DataFrame:
    (
        NON_DISEASE_THERAPEUTIC_AREAS,
        load_and_aggregate_evidence,
        load_disease_terms,
        load_unmet_needs,
    ) = import_build_mapping_helpers()

    variants = variants[~variants["ensembl_gene_id"].isin(eur_genes)].copy()
    disease_terms = load_disease_terms(disease_path, NON_DISEASE_THERAPEUTIC_AREAS)
    unmet_needs = load_unmet_needs(unmet_needs_path)
    disease_gene = load_and_aggregate_evidence(evidence_path, disease_terms, unmet_needs)
    disease_gene = disease_gene[disease_gene["ot_score"] >= min_score].copy()
    disease_gene = disease_gene[~disease_gene["ensembl_gene_id"].isin(eur_genes)].copy()

    mapped = disease_gene.merge(variants, on="ensembl_gene_id", how="inner")
    mapped["gene_symbol"] = mapped["gene_symbol"].fillna("")
    mapped = mapped.sort_values(
        ["gene_symbol", "ensembl_gene_id", "ot_score", "disease_name", "chrom_sort", "pos"],
        ascending=[True, True, False, True, True, True],
        kind="mergesort",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapped.to_parquet(output_path, index=False)

    unique_disease_gene = mapped.drop_duplicates(["ensembl_gene_id", "disease_id"])
    stats = {
        "source_annotation": "VEP 115.1 coding-only VCF",
        "min_l2g_score": min_score,
        "mapped_rows": int(len(mapped)),
        "genes": int(mapped["ensembl_gene_id"].nunique()),
        "lof_positive_genes": int(mapped.loc[mapped["gene_has_lof"], "ensembl_gene_id"].nunique()),
        "disease_traits": int(mapped["disease_id"].nunique()),
        "disease_gene_pairs": int(unique_disease_gene.shape[0]),
        "variants": int(mapped["varid"].nunique()),
        "plof_variants_not_common_sas": int(mapped.loc[mapped["plof_is_not_common_sas"], "varid"].nunique()),
        "deleterious_missense_variants": int(
            mapped.loc[mapped["variant_is_deleterious_missense"], "varid"].nunique()
        ),
        "global_maf_variants": int(mapped.loc[mapped["variant_selected_for_global_maf"], "varid"].nunique()),
        "lof_variants": int(mapped.loc[mapped["plof_is_not_common_sas"], "varid"].nunique()),
        "genes_with_deleterious_missense": int(
            mapped.loc[mapped["gene_has_deleterious_missense"], "ensembl_gene_id"].nunique()
        ),
        "max_plof_sas_af": max_plof_sas_af,
        "cadd_missense_threshold": cadd_missense_threshold,
        "alphamissense_threshold": alphamissense_threshold,
        "esm1b_threshold": esm1b_threshold,
        "esm1b_available": bool(mapped["ESM1b_score"].notna().any()) if "ESM1b_score" in mapped.columns else False,
        "median_l2g_score_disease_gene_pairs": float(unique_disease_gene["ot_score"].median()) if len(mapped) else None,
        "exclude_eur_rare_genes": True,
        "eur_rare_af_max": eur_rare_af_max,
        "eur_rare_populations": ["AF_nfe", "AF_fin"],
        "eur_rare_excluded_genes": len(eur_genes),
        "disease_traits_with_unmet_needs": int(mapped.loc[mapped["unmet_need_index"].notna(), "disease_id"].nunique()),
        "disease_gene_pairs_with_unmet_needs": int(
            mapped.loc[mapped["unmet_need_index"].notna(), ["ensembl_gene_id", "disease_id"]]
            .drop_duplicates()
            .shape[0]
        ),
    }
    output_path.with_suffix(output_path.suffix + ".stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    return mapped


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a VEP 115 disease-gene-variant mapping for the static viewer.")
    parser.add_argument("--vep", type=Path, default=DEFAULT_VEP)
    parser.add_argument("--gnomad", type=Path, default=DEFAULT_GNOMAD)
    parser.add_argument("--genome-india", type=Path, default=DEFAULT_GENOME_INDIA)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--disease", type=Path, default=DEFAULT_DISEASE)
    parser.add_argument("--unmet-needs", type=Path, default=DEFAULT_UNMET_NEEDS)
    parser.add_argument("--annotated-output", type=Path, default=DEFAULT_ANNOTATED)
    parser.add_argument("--high-moderate-output", type=Path, default=DEFAULT_HIGH_MODERATE)
    parser.add_argument("--eur-genes-output", type=Path, default=DEFAULT_EUR_GENES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--min-score", type=float, default=0.15)
    parser.add_argument("--eur-rare-af-max", type=float, default=0.001)
    parser.add_argument(
        "--max-plof-sas-af",
        type=float,
        default=0.3,
        help="Retain pLoF only when SAS AF is less than or equal to this value.",
    )
    parser.add_argument(
        "--cadd-missense-threshold",
        type=float,
        default=20.0,
        help="CADD PHRED threshold for deleterious missense classification.",
    )
    parser.add_argument(
        "--alphamissense-threshold",
        type=float,
        default=0.564,
        help="AlphaMissense pathogenicity threshold for deleterious missense classification.",
    )
    parser.add_argument(
        "--esm1b-threshold",
        type=float,
        default=None,
        help="Optional ESM1b threshold if an ESM1b score column is available.",
    )
    args = parser.parse_args()

    annotated = annotate_vep_with_frequency(args.vep, args.gnomad, args.genome_india, args.annotated_output)
    print(f"Wrote annotated VEP/frequency table: {args.annotated_output} ({annotated.height:,} rows)")

    variants = prepare_high_moderate_variants(
        args.annotated_output,
        args.high_moderate_output,
        max_plof_sas_af=args.max_plof_sas_af,
        cadd_missense_threshold=args.cadd_missense_threshold,
        alphamissense_threshold=args.alphamissense_threshold,
        esm1b_threshold=args.esm1b_threshold,
    )
    print(f"Wrote selected pLoF/deleterious missense variants: {args.high_moderate_output} ({len(variants):,} rows)")

    eur_genes = write_eur_excluded_genes(variants, args.eur_genes_output, args.eur_rare_af_max)
    print(f"Wrote EUR rare excluded genes: {args.eur_genes_output} ({len(eur_genes):,} genes)")

    mapped = build_disease_mapping(
        variants=variants,
        evidence_path=args.evidence,
        disease_path=args.disease,
        unmet_needs_path=args.unmet_needs,
        output_path=args.output,
        eur_genes=eur_genes,
        min_score=args.min_score,
        eur_rare_af_max=args.eur_rare_af_max,
        max_plof_sas_af=args.max_plof_sas_af,
        cadd_missense_threshold=args.cadd_missense_threshold,
        alphamissense_threshold=args.alphamissense_threshold,
        esm1b_threshold=args.esm1b_threshold,
    )
    print(f"Wrote final disease-gene-variant map: {args.output} ({len(mapped):,} rows)")
    print(f"Genes: {mapped['ensembl_gene_id'].nunique():,}")
    print(f"Diseases: {mapped['disease_id'].nunique():,}")
    print(f"Variants: {mapped['varid'].nunique():,}")


if __name__ == "__main__":
    main()
