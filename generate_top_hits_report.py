#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
from pathlib import Path

import pandas as pd
import polars as pl


APP_DIR = Path(__file__).resolve().parent
DEFAULT_PARQUET = (
    APP_DIR
    / "data"
    / "ot_sas_gi_vep115_plof_deleterious_missense_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet"
)
DEFAULT_HTML = APP_DIR / "reports" / "top_100_hits.html"
DEFAULT_OT_REF = Path("/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/opentargets_26_06")


def esc(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return html.escape(str(value))


def fmt(value: object, digits: int = 4) -> str:
    if value is None or pd.isna(value):
        return ""
    try:
        return f"{float(value):.{digits}g}"
    except Exception:
        return esc(value)


def first_non_empty(series: pd.Series) -> str:
    for value in series.dropna():
        text = str(value).strip()
        if text:
            return text
    return ""


def variant_table(df: pd.DataFrame) -> str:
    keep = [
        "varid",
        "consequence",
        "missense_interpretation",
        "AF_sas",
        "genome_india_af",
        "AF_nfe",
        "AF_fin",
        "cadd_phred",
        "am_pathogenicity",
        "EVE_SCORE",
        "REVEL",
        "ClinPred",
        "global_maf_af_sas",
        "global_maf_genome_india_af",
    ]
    rows = []
    for _, row in df.drop_duplicates(["varid", "consequence", "missense_interpretation"]).sort_values(
        ["chrom_sort", "pos", "varid"], kind="mergesort"
    ).iterrows():
        rows.append(
            "<tr>"
            + "".join(
                f"<td>{fmt(row[col]) if col not in ['varid', 'consequence', 'missense_interpretation'] else esc(row[col])}</td>"
                for col in keep
                if col in row
            )
            + "</tr>"
        )
    headers = "".join(f"<th>{esc(col)}</th>" for col in keep if col in df.columns)
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def load_target_reference(ref_dir: Path, target_ids: set[str]) -> dict[str, dict[str, str]]:
    target = (
        pl.scan_parquet(str(ref_dir / "target" / "*.parquet"))
        .filter(pl.col("id").is_in(sorted(target_ids)))
        .select("id", "approvedName", "functionDescriptions", "pathways")
        .collect()
        .to_pandas()
    )
    result = {}
    for row in target.itertuples(index=False):
        funcs = row.functionDescriptions if isinstance(row.functionDescriptions, list) else []
        pathways = row.pathways if isinstance(row.pathways, list) else []
        pathway_names = []
        for p in pathways[:8]:
            if isinstance(p, dict) and p.get("pathway"):
                pathway_names.append(str(p["pathway"]))
        result[row.id] = {
            "approved_name": row.approvedName or "",
            "function": " ".join(str(x) for x in funcs[:3]) if funcs else "",
            "pathways": "; ".join(pathway_names),
        }
    return result


def load_drug_reference(ref_dir: Path, target_ids: set[str], disease_ids: set[str]) -> dict[tuple[str, str], str]:
    moa = (
        pl.scan_parquet(str(ref_dir / "drug_mechanism_of_action" / "*.parquet"))
        .select("actionType", "mechanismOfAction", "chemblIds", "targetName", "targets")
        .explode("targets")
        .filter(pl.col("targets").is_in(sorted(target_ids)))
        .explode("chemblIds")
        .rename({"targets": "targetId", "chemblIds": "drugId"})
    )
    indications = (
        pl.scan_parquet(str(ref_dir / "clinical_indication" / "*.parquet"))
        .filter(pl.col("diseaseId").is_in(sorted(disease_ids)))
        .select("drugId", "diseaseId", "maxClinicalStage")
    )
    molecules = pl.scan_parquet(str(ref_dir / "drug_molecule" / "*.parquet")).select(
        "id", "name", "maximumClinicalStage", "description"
    )
    joined = (
        moa.join(indications, on="drugId", how="inner")
        .join(molecules, left_on="drugId", right_on="id", how="left")
        .collect()
        .to_pandas()
    )
    result: dict[tuple[str, str], list[str]] = {}
    for row in joined.itertuples(index=False):
        key = (row.targetId, row.diseaseId)
        drug = row.name or row.drugId
        stage = row.maxClinicalStage or row.maximumClinicalStage or ""
        mechanism = row.mechanismOfAction or ""
        item = f"{drug} ({stage}; {mechanism})"
        result.setdefault(key, [])
        if item not in result[key]:
            result[key].append(item)
    return {key: "; ".join(values[:8]) for key, values in result.items()}


def load_biology_reference(ref_dir: Path, target_ids: set[str], disease_ids: set[str]) -> dict[tuple[str, str], str]:
    assoc = (
        pl.scan_parquet(str(ref_dir / "association_by_datasource_direct" / "*.parquet"))
        .filter(pl.col("targetId").is_in(sorted(target_ids)) & pl.col("diseaseId").is_in(sorted(disease_ids)))
        .select("diseaseId", "targetId", "aggregationValue", "associationScore", "evidenceCount")
        .collect()
        .to_pandas()
    )
    if assoc.empty:
        return {}
    assoc = assoc.sort_values(["targetId", "diseaseId", "associationScore"], ascending=[True, True, False])
    result: dict[tuple[str, str], list[str]] = {}
    for row in assoc.itertuples(index=False):
        key = (row.targetId, row.diseaseId)
        item = f"{row.aggregationValue}: score {fmt(row.associationScore, 3)}, evidence {fmt(row.evidenceCount, 0)}"
        result.setdefault(key, [])
        if item not in result[key]:
            result[key].append(item)
    return {key: "; ".join(values[:8]) for key, values in result.items()}


def build_report(parquet_path: Path, html_path: Path, limit: int, ot_ref_dir: Path) -> None:
    df = pd.read_parquet(parquet_path)
    df["unmet_rank"] = pd.to_numeric(df.get("unmet_need_index"), errors="coerce").fillna(-1)
    df["rank_maf"] = pd.to_numeric(df.get("global_maf_af_sas"), errors="coerce").fillna(0)
    pair_rank = (
        df.drop_duplicates(["ensembl_gene_id", "disease_id"])
        .sort_values(
            ["unmet_rank", "ot_score", "rank_maf", "gene_symbol", "disease_name"],
            ascending=[False, False, False, True, True],
            kind="mergesort",
        )
        .head(limit)
    )
    target_ids = set(pair_rank["ensembl_gene_id"].dropna())
    disease_ids = set(pair_rank["disease_id"].dropna())
    target_ref = load_target_reference(ot_ref_dir, target_ids)
    drug_ref = load_drug_reference(ot_ref_dir, target_ids, disease_ids)
    biology_ref = load_biology_reference(ot_ref_dir, target_ids, disease_ids)

    css = """
    @page { size: letter; margin: 0.65in; }
    body { font-family: Arial, sans-serif; color: #17202a; font-size: 10.5pt; line-height: 1.35; }
    h1 { font-size: 22pt; margin: 0 0 12px; }
    h2 { font-size: 16pt; margin: 0 0 8px; }
    h3 { font-size: 12pt; margin: 12px 0 5px; }
    .muted { color: #52616f; }
    .toc a { color: #075985; text-decoration: none; }
    .toc li { margin: 4px 0; }
    .page { break-before: page; page-break-before: always; }
    .grid { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 8px; margin: 8px 0 10px; }
    .metric { border: 1px solid #d7dce2; border-radius: 4px; padding: 6px; }
    .metric span { display: block; color: #52616f; font-size: 8.5pt; }
    .metric strong { display: block; font-size: 11pt; }
    table { width: 100%; border-collapse: collapse; table-layout: fixed; margin-top: 6px; }
    th, td { border: 1px solid #d7dce2; padding: 4px; vertical-align: top; overflow-wrap: anywhere; }
    th { background: #f1f5f9; font-size: 8pt; }
    td { font-size: 8pt; }
    .note { border-left: 3px solid #94a3b8; padding-left: 8px; color: #475569; }
    """

    toc_items = []
    pages = []
    for i, pair in enumerate(pair_rank.itertuples(index=False), start=1):
        anchor = f"hit-{i}"
        pair_df = df[(df["ensembl_gene_id"] == pair.ensembl_gene_id) & (df["disease_id"] == pair.disease_id)]
        title = f"{pair.gene_symbol or pair.ensembl_gene_id} - {pair.disease_name}"
        key = (pair.ensembl_gene_id, pair.disease_id)
        toc_items.append(
            f'<li><a href="#{anchor}">{i}. {esc(title)}</a> '
            f'<span class="muted">L2G {fmt(pair.ot_score, 3)}; unmet {fmt(getattr(pair, "unmet_need_index", None), 3)}</span></li>'
        )
        target_info = target_ref.get(pair.ensembl_gene_id, {})
        gene_summary = target_info.get("function") or target_info.get("approved_name") or (
            "No Open Targets target function description available."
        )
        if target_info.get("pathways"):
            gene_summary += f" Pathways: {target_info['pathways']}."
        drugs = drug_ref.get(key) or (
            "No same target-disease clinical indication found in the downloaded Open Targets drug tables."
        )
        biology = biology_ref.get(key)
        if biology:
            biology = (
                f"Datasource-level direct associations: {biology}. GWAS credible-set evidence in this analysis has "
                f"L2G {fmt(pair.ot_score, 3)}, {fmt(getattr(pair, 'evidence_count', None), 0)} evidence rows, and "
                f"{fmt(getattr(pair, 'study_locus_count', None), 0)} study loci."
            )
        else:
            biology = (
                f"GWAS credible-set evidence supports this gene-disease pair with L2G score {fmt(pair.ot_score, 3)} "
                f"across {fmt(getattr(pair, 'evidence_count', None), 0)} evidence rows and "
                f"{fmt(getattr(pair, 'study_locus_count', None), 0)} study loci. No datasource-level direct association "
                f"summary was found in the downloaded association table."
            )
        prevalence = first_non_empty(pair_df.get("unmet_need_prevalence", pd.Series(dtype=object)))
        if not prevalence:
            prevalence = "No SAS/Indian prevalence estimate is available in the current local cache."

        pages.append(
            f"""
            <section class="page" id="{anchor}">
              <h2>{i}. {esc(title)}</h2>
              <div class="muted">{esc(pair.ensembl_gene_id)} | {esc(pair.disease_id)}</div>
              <div class="grid">
                <div class="metric"><span>L2G score</span><strong>{fmt(pair.ot_score, 3)}</strong></div>
                <div class="metric"><span>Unmet need index</span><strong>{fmt(getattr(pair, 'unmet_need_index', None), 3)}</strong></div>
                <div class="metric"><span>SAS global MAF</span><strong>{fmt(pair_df['global_maf_af_sas'].max(), 4)}</strong></div>
                <div class="metric"><span>Genome India global MAF</span><strong>{fmt(pair_df['global_maf_genome_india_af'].max(), 4)}</strong></div>
                <div class="metric"><span>Variants</span><strong>{pair_df['varid'].nunique()}</strong></div>
                <div class="metric"><span>Interpretations</span><strong>{esc('; '.join(sorted(set(pair_df['missense_interpretation'].dropna()))))}</strong></div>
              </div>
              <h3>Disease Description</h3>
              <p>{esc(first_non_empty(pair_df['disease_description']))}</p>
              <h3>Gene Summary</h3>
              <p class="note">{esc(gene_summary)}</p>
              <h3>Known Drugs / Clinical Investigation</h3>
              <p class="note">{esc(drugs)}</p>
              <h3>Biological Evidence</h3>
              <p>{esc(biology)}</p>
              <h3>SAS / Indian Prevalence</h3>
              <p>{esc(prevalence)}</p>
              <h3>Alleles and Scores</h3>
              {variant_table(pair_df)}
            </section>
            """
        )

    html_text = f"""
    <!doctype html>
    <html>
      <head>
        <meta charset="utf-8" />
        <title>Top 100 SAS Variant-Disease Hits</title>
        <style>{css}</style>
      </head>
      <body>
        <h1>Top {limit} SAS Variant-Disease Hits</h1>
        <p class="muted">Ranked by unmet-need index, L2G score, and SAS global MAF. Generated from {esc(parquet_path.name)}.</p>
        <h2>Table of Contents</h2>
        <ol class="toc">{''.join(toc_items)}</ol>
        {''.join(pages)}
      </body>
    </html>
    """
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(html_text)
    print(f"Wrote HTML report: {html_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate top-hit HTML report for PDF printing.")
    parser.add_argument("--parquet", type=Path, default=DEFAULT_PARQUET)
    parser.add_argument("--html", type=Path, default=DEFAULT_HTML)
    parser.add_argument("--ot-ref-dir", type=Path, default=DEFAULT_OT_REF)
    parser.add_argument("--limit", type=int, default=100)
    args = parser.parse_args()
    build_report(args.parquet, args.html, args.limit, args.ot_ref_dir)


if __name__ == "__main__":
    main()
