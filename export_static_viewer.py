#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import pandas as pd


APP_DIR = Path(__file__).resolve().parent
DEFAULT_PARQUET = APP_DIR / "data" / "ot_sas_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet"
DEFAULT_STATIC_DIR = APP_DIR / "static"
DEFAULT_DB = DEFAULT_STATIC_DIR / "data" / "ot_sas_viewer.sqlite"


def write_table(conn: sqlite3.Connection, name: str, df: pd.DataFrame) -> None:
    df.to_sql(name, conn, if_exists="replace", index=False)


def build_sqlite(parquet_path: Path, db_path: Path) -> dict[str, int | float | str | bool | None]:
    df = pd.read_parquet(parquet_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()

    unique_disease_gene = df.drop_duplicates(["ensembl_gene_id", "disease_id"]).copy()

    genes = (
        df.groupby(["ensembl_gene_id", "gene_symbol"], dropna=False)
        .agg(
            disease_count=("disease_id", "nunique"),
            variant_count=("varid", "nunique"),
            max_score=("ot_score", "max"),
            gene_has_lof=("gene_has_lof", "max"),
            cum_af_sas_lof=("cum_af_sas_lof", "max"),
            cum_af_sas_missense=("cum_af_sas_missense", "max"),
        )
        .reset_index()
    )
    genes["gene_has_lof"] = genes["gene_has_lof"].astype(int)
    genes["gene_symbol"] = genes["gene_symbol"].fillna("")
    genes["symbol_missing"] = genes["gene_symbol"].eq("").astype(int)
    genes = genes.sort_values(["symbol_missing", "gene_symbol", "ensembl_gene_id"], kind="mergesort").drop(
        columns=["symbol_missing"]
    )

    disease_gene = unique_disease_gene[
        [
            "ensembl_gene_id",
            "disease_id",
            "disease_name",
            "disease_description",
            "ot_score",
            "unmet_need_index",
            "unmet_need_category",
            "unmet_need_prevalence",
            "unmet_need_key_needs",
            "evidence_count",
            "study_locus_count",
        ]
    ].sort_values(["ensembl_gene_id", "ot_score", "disease_name"], ascending=[True, False, True], kind="mergesort")

    variants = (
        df[
            [
                "ensembl_gene_id",
                "chrom",
                "chrom_sort",
                "pos",
                "varid",
                "AF_sas",
                "AF_nfe",
                "cadd_phred",
                "sas_enrichment",
                "consequence",
            ]
        ]
        .drop_duplicates(
            ["ensembl_gene_id", "chrom", "pos", "varid", "AF_sas", "AF_nfe", "cadd_phred", "sas_enrichment", "consequence"]
        )
        .sort_values(["ensembl_gene_id", "chrom_sort", "pos", "varid"], kind="mergesort")
    )

    stats = {
        "mapped_rows": int(len(df)),
        "genes": int(df["ensembl_gene_id"].nunique()),
        "lof_positive_genes": int(df.loc[df["gene_has_lof"], "ensembl_gene_id"].nunique()),
        "disease_traits": int(df["disease_id"].nunique()),
        "disease_gene_pairs": int(unique_disease_gene.shape[0]),
        "variants": int(df["varid"].nunique()),
        "lof_variants": int(df.loc[df["variant_is_lof"], "varid"].nunique()),
        "median_l2g_score_disease_gene_pairs": float(unique_disease_gene["ot_score"].median()),
        "min_l2g_score": float(df["ot_score"].min()),
        "disease_traits_with_unmet_needs": int(df.loc[df["unmet_need_index"].notna(), "disease_id"].nunique()),
        "disease_gene_pairs_with_unmet_needs": int(
            df.loc[df["unmet_need_index"].notna(), ["ensembl_gene_id", "disease_id"]].drop_duplicates().shape[0]
        ),
    }
    stats_json = parquet_path.with_suffix(parquet_path.suffix + ".stats.json")
    if stats_json.exists():
        stats.update(json.loads(stats_json.read_text()))

    conn = sqlite3.connect(db_path)
    try:
        write_table(conn, "genes", genes)
        write_table(conn, "disease_gene", disease_gene)
        write_table(conn, "variants", variants)
        stats_df = pd.DataFrame([{"key": key, "value": json.dumps(value)} for key, value in stats.items()])
        write_table(conn, "stats", stats_df)
        conn.executescript(
            """
            CREATE INDEX idx_genes_symbol ON genes(gene_symbol, ensembl_gene_id);
            CREATE INDEX idx_disease_gene_gene ON disease_gene(ensembl_gene_id, ot_score DESC);
            CREATE INDEX idx_variants_gene ON variants(ensembl_gene_id, chrom_sort, pos, varid);
            VACUUM;
            """
        )
    finally:
        conn.close()
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description="Export a static SQLite database for the browser-only viewer.")
    parser.add_argument("--parquet", type=Path, default=DEFAULT_PARQUET)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    args = parser.parse_args()

    stats = build_sqlite(args.parquet, args.db)
    print(f"Wrote SQLite database: {args.db}")
    print(f"Genes: {stats['genes']:,}")
    print(f"Disease-gene pairs: {stats['disease_gene_pairs']:,}")
    print(f"Variants: {stats['variants']:,}")


if __name__ == "__main__":
    main()
