#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import dash
import pandas as pd
from dash import Input, Output, State, dash_table, dcc, html
from flask import send_file


APP_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_PATH = APP_DIR / "data" / "ot_sas_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet"
DATA_PATH = Path(os.environ.get("OT_SAS_MAPPING_PATH", DEFAULT_DATA_PATH))


def load_mapping() -> pd.DataFrame:
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"Mapping file not found: {DATA_PATH}. Run build_mapping.py before starting the app."
        )
    df = pd.read_parquet(DATA_PATH)
    text_cols = ["gene_symbol", "disease_name", "disease_description", "consequence"]
    for col in text_cols:
        df[col] = df[col].fillna("")
    return df


def make_gene_table(df: pd.DataFrame) -> pd.DataFrame:
    gene = (
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
    gene["symbol_missing"] = gene["gene_symbol"].fillna("").eq("")
    gene = gene.sort_values(["symbol_missing", "gene_symbol", "ensembl_gene_id"], kind="mergesort").drop(
        columns=["symbol_missing"]
    )
    return gene


def round_float(value: object, digits: int = 6) -> object:
    if pd.isna(value):
        return None
    return round(float(value), digits)


mapping_df = load_mapping()
gene_df = make_gene_table(mapping_df)
unique_disease_gene_df = mapping_df.drop_duplicates(["ensembl_gene_id", "disease_id"])
APP_STATS = {
    "min_l2g_score": float(mapping_df["ot_score"].min()) if len(mapping_df) else None,
    "genes": int(mapping_df["ensembl_gene_id"].nunique()),
    "lof_positive_genes": int(mapping_df.loc[mapping_df["gene_has_lof"], "ensembl_gene_id"].nunique()),
    "disease_traits": int(mapping_df["disease_id"].nunique()),
    "disease_gene_pairs": int(unique_disease_gene_df.shape[0]),
    "variants": int(mapping_df["varid"].nunique()),
    "lof_variants": int(mapping_df.loc[mapping_df["variant_is_lof"], "varid"].nunique()),
    "mapped_rows": int(len(mapping_df)),
    "median_l2g_score": float(unique_disease_gene_df["ot_score"].median()) if len(unique_disease_gene_df) else None,
}


def gene_records(df: pd.DataFrame) -> list[dict]:
    records = df.to_dict("records")
    for record in records:
        record["id"] = record["ensembl_gene_id"]
        record["gene_has_lof"] = "yes" if record.get("gene_has_lof") else ""
    return records


app = dash.Dash(__name__, title="Open Targets SAS Variant Browser")
server = app.server


@server.route("/download/mapping.parquet")
def download_mapping():
    return send_file(DATA_PATH, as_attachment=True, download_name=DATA_PATH.name)

BASE_TABLE_STYLE = {
    "fontFamily": "Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
    "fontSize": "13px",
    "border": "1px solid #d7dce2",
}

CELL_STYLE = {
    "padding": "7px 9px",
    "whiteSpace": "normal",
    "height": "auto",
    "lineHeight": "1.25",
    "textAlign": "left",
    "border": "1px solid #e2e6ea",
}

HEADER_STYLE = {
    "backgroundColor": "#eef2f5",
    "fontWeight": "700",
    "border": "1px solid #ccd3db",
}


def data_table(table_id: str, columns: list[dict], page_size: int, **kwargs) -> dash_table.DataTable:
    data = kwargs.pop("data", [])
    extra_style_data_conditional = kwargs.pop("style_data_conditional", [])
    return dash_table.DataTable(
        id=table_id,
        columns=columns,
        data=data,
        page_size=page_size,
        sort_action="native",
        filter_action="native",
        row_selectable=False,
        style_table={"overflowX": "auto", "border": "1px solid #d7dce2"},
        style_cell=CELL_STYLE,
        style_header=HEADER_STYLE,
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#f8fafc"},
            {"if": {"state": "active"}, "backgroundColor": "#dbeafe", "border": "1px solid #2563eb"},
        ]
        + extra_style_data_conditional,
        css=[{"selector": ".dash-spreadsheet-container .dash-spreadsheet-inner table", "rule": "table-layout: fixed;"}],
        **kwargs,
    )


app.layout = html.Div(
    [
        html.Div(
            [
                html.H1("Open Targets SAS Variant Browser"),
                html.Div(
                    [
                        html.Div(f"{len(gene_df):,} genes"),
                        html.Div(f"{mapping_df['disease_id'].nunique():,} disease traits"),
                        html.Div(f"{mapping_df['varid'].nunique():,} variants"),
                        html.Div(f"L2G >= {APP_STATS['min_l2g_score']:.2f}"),
                    ],
                    className="metrics",
                ),
            ],
            className="topbar",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Label("Gene search"),
                        dcc.Input(
                            id="gene-search",
                            type="text",
                            debounce=True,
                            placeholder="Symbol or Ensembl ID",
                            className="search",
                        ),
                        data_table(
                            "gene-table",
                            [
                                {"name": "Ensembl gene id", "id": "ensembl_gene_id"},
                                {"name": "symbol", "id": "gene_symbol"},
                                {"name": "LoF", "id": "gene_has_lof"},
                                {
                                    "name": "cum AF LoF",
                                    "id": "cum_af_sas_lof",
                                    "type": "numeric",
                                    "format": {"specifier": ".4g"},
                                },
                                {
                                    "name": "cum AF missense",
                                    "id": "cum_af_sas_missense",
                                    "type": "numeric",
                                    "format": {"specifier": ".4g"},
                                },
                                {"name": "diseases", "id": "disease_count", "type": "numeric"},
                                {"name": "variants", "id": "variant_count", "type": "numeric"},
                                {
                                    "name": "max score",
                                    "id": "max_score",
                                    "type": "numeric",
                                    "format": {"specifier": ".3f"},
                                },
                            ],
                            page_size=18,
                            data=gene_records(gene_df),
                            style_cell_conditional=[
                                {"if": {"column_id": "ensembl_gene_id"}, "width": "180px"},
                                {"if": {"column_id": "gene_symbol"}, "width": "100px"},
                                {"if": {"column_id": "gene_has_lof"}, "width": "52px"},
                                {"if": {"column_id": "cum_af_sas_lof"}, "width": "95px", "textAlign": "right"},
                                {"if": {"column_id": "cum_af_sas_missense"}, "width": "120px", "textAlign": "right"},
                                {"if": {"column_id": "disease_count"}, "width": "80px", "textAlign": "right"},
                                {"if": {"column_id": "variant_count"}, "width": "80px", "textAlign": "right"},
                                {"if": {"column_id": "max_score"}, "width": "90px", "textAlign": "right"},
                            ],
                            style_data_conditional=[
                                {
                                    "if": {"filter_query": '{gene_has_lof} = "yes"'},
                                    "backgroundColor": "#dbeafe",
                                    "color": "#0f172a",
                                    "fontWeight": "700",
                                }
                            ],
                        ),
                    ],
                    className="left-pane",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Span("Disease-gene pairs"),
                                        html.Strong(f"{APP_STATS['disease_gene_pairs']:,}"),
                                    ]
                                ),
                                html.Div(
                                    [
                                        html.Span("LoF-positive genes"),
                                        html.Strong(f"{APP_STATS['lof_positive_genes']:,}"),
                                    ]
                                ),
                                html.Div(
                                    [
                                        html.Span("LoF variants"),
                                        html.Strong(f"{APP_STATS['lof_variants']:,}"),
                                    ]
                                ),
                                html.Div(
                                    [
                                        html.Span("Median L2G"),
                                        html.Strong(f"{APP_STATS['median_l2g_score']:.3f}"),
                                    ]
                                ),
                                html.Div(
                                    [
                                        html.Span("Mapped rows"),
                                        html.Strong(f"{APP_STATS['mapped_rows']:,}"),
                                    ]
                                ),
                                html.A("Download mapped parquet", href="/download/mapping.parquet", className="download"),
                            ],
                            className="stats",
                        ),
                        html.Div(id="selected-gene", className="selected"),
                        html.H2("Associated diseases"),
                        data_table(
                            "disease-table",
                            [
                                {"name": "disease", "id": "disease_name"},
                                {"name": "disease id", "id": "disease_id"},
                                {
                                    "name": "score",
                                    "id": "ot_score",
                                    "type": "numeric",
                                    "format": {"specifier": ".3f"},
                                },
                                {"name": "evidence rows", "id": "evidence_count", "type": "numeric"},
                                {"name": "study loci", "id": "study_locus_count", "type": "numeric"},
                                {"name": "description", "id": "disease_description"},
                            ],
                            page_size=10,
                            style_cell_conditional=[
                                {"if": {"column_id": "disease_name"}, "width": "210px"},
                                {"if": {"column_id": "disease_id"}, "width": "130px"},
                                {"if": {"column_id": "ot_score"}, "width": "80px", "textAlign": "right"},
                                {"if": {"column_id": "evidence_count"}, "width": "90px", "textAlign": "right"},
                                {"if": {"column_id": "study_locus_count"}, "width": "80px", "textAlign": "right"},
                                {"if": {"column_id": "disease_description"}, "width": "420px"},
                            ],
                        ),
                        html.H2("Mapped coding variants"),
                        data_table(
                            "variant-table",
                            [
                                {"name": "chrom", "id": "chrom"},
                                {"name": "pos", "id": "pos", "type": "numeric"},
                                {"name": "varid", "id": "varid"},
                                {"name": "AF_sas", "id": "AF_sas", "type": "numeric", "format": {"specifier": ".6g"}},
                                {"name": "AF_nfe", "id": "AF_nfe", "type": "numeric", "format": {"specifier": ".6g"}},
                                {"name": "CADD", "id": "cadd_phred", "type": "numeric", "format": {"specifier": ".3g"}},
                                {
                                    "name": "sas_enrichment",
                                    "id": "sas_enrichment",
                                    "type": "numeric",
                                    "format": {"specifier": ".3f"},
                                },
                                {"name": "consequence", "id": "consequence"},
                            ],
                            page_size=12,
                            style_cell_conditional=[
                                {"if": {"column_id": "chrom"}, "width": "60px"},
                                {"if": {"column_id": "pos"}, "width": "100px", "textAlign": "right"},
                                {"if": {"column_id": "varid"}, "width": "150px"},
                                {"if": {"column_id": "AF_sas"}, "width": "90px", "textAlign": "right"},
                                {"if": {"column_id": "AF_nfe"}, "width": "90px", "textAlign": "right"},
                                {"if": {"column_id": "cadd_phred"}, "width": "80px", "textAlign": "right"},
                                {"if": {"column_id": "sas_enrichment"}, "width": "110px", "textAlign": "right"},
                                {"if": {"column_id": "consequence"}, "width": "190px"},
                            ],
                        ),
                    ],
                    className="right-pane",
                ),
            ],
            className="workspace",
        ),
    ]
)

app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            body {
                margin: 0;
                background: #f5f7fa;
                color: #17202a;
                font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            }
            .topbar {
                display: flex;
                justify-content: space-between;
                align-items: center;
                gap: 24px;
                padding: 18px 22px;
                background: #ffffff;
                border-bottom: 1px solid #d7dce2;
            }
            h1 {
                margin: 0;
                font-size: 22px;
                font-weight: 750;
                letter-spacing: 0;
            }
            h2 {
                margin: 18px 0 8px;
                font-size: 15px;
                font-weight: 750;
                letter-spacing: 0;
            }
            label {
                display: block;
                margin: 0 0 8px;
                font-size: 13px;
                font-weight: 700;
            }
            .metrics {
                display: flex;
                gap: 10px;
                flex-wrap: wrap;
                font-size: 13px;
                color: #34495e;
            }
            .metrics div {
                padding: 6px 9px;
                border: 1px solid #d7dce2;
                background: #f8fafc;
                border-radius: 6px;
            }
            .workspace {
                display: grid;
                grid-template-columns: minmax(380px, 34vw) minmax(0, 1fr);
                gap: 18px;
                padding: 18px;
                height: calc(100vh - 82px);
                box-sizing: border-box;
            }
            .left-pane, .right-pane {
                min-width: 0;
                overflow: auto;
            }
            .search {
                width: 100%;
                box-sizing: border-box;
                border: 1px solid #b9c2cc;
                border-radius: 6px;
                padding: 9px 10px;
                margin-bottom: 12px;
                font-size: 14px;
                background: #ffffff;
            }
            .selected {
                min-height: 24px;
                font-size: 14px;
                font-weight: 700;
                color: #0f172a;
            }
            .stats {
                display: grid;
                grid-template-columns: repeat(5, minmax(120px, 1fr)) minmax(180px, auto);
                gap: 10px;
                align-items: stretch;
                margin-bottom: 14px;
            }
            .stats div, .download {
                background: #ffffff;
                border: 1px solid #d7dce2;
                border-radius: 6px;
                padding: 9px 10px;
                min-height: 44px;
                box-sizing: border-box;
            }
            .stats span {
                display: block;
                font-size: 12px;
                color: #52616f;
                margin-bottom: 4px;
            }
            .stats strong {
                display: block;
                font-size: 16px;
            }
            .download {
                display: flex;
                align-items: center;
                justify-content: center;
                color: #075985;
                font-size: 13px;
                font-weight: 700;
                text-decoration: none;
            }
            @media (max-width: 900px) {
                .topbar {
                    display: block;
                }
                .metrics {
                    margin-top: 12px;
                }
                .workspace {
                    display: block;
                    height: auto;
                }
                .right-pane {
                    margin-top: 18px;
                }
                .stats {
                    grid-template-columns: 1fr 1fr;
                }
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""


@app.callback(
    Output("gene-table", "data"),
    Input("gene-search", "value"),
)
def filter_genes(query: str | None) -> list[dict]:
    if not query:
        data = gene_df
    else:
        q = query.strip().casefold()
        data = gene_df[
            gene_df["ensembl_gene_id"].str.casefold().str.contains(q, regex=False)
            | gene_df["gene_symbol"].fillna("").str.casefold().str.contains(q, regex=False)
        ]
    return gene_records(data)


@app.callback(
    Output("selected-gene", "children"),
    Output("disease-table", "data"),
    Output("variant-table", "data"),
    Input("gene-table", "active_cell"),
    State("gene-table", "derived_virtual_data"),
    State("gene-table", "data"),
)
def update_gene_details(active_cell: dict | None, virtual_rows: list[dict] | None, table_rows: list[dict] | None):
    rows = virtual_rows or table_rows or []
    if not rows:
        return "No gene selected", [], []
    if active_cell and active_cell.get("row_id"):
        ensembl_gene_id = active_cell["row_id"]
        selected = next((row for row in rows if row.get("ensembl_gene_id") == ensembl_gene_id), rows[0])
    else:
        row_index = active_cell["row"] if active_cell else 0
        if row_index >= len(rows):
            row_index = 0
        selected = rows[row_index]
        ensembl_gene_id = selected["ensembl_gene_id"]
    gene_symbol = selected.get("gene_symbol") or ""

    gene_rows = mapping_df[mapping_df["ensembl_gene_id"] == ensembl_gene_id]
    diseases = (
        gene_rows[
            [
                "disease_id",
                "disease_name",
                "disease_description",
                "ot_score",
                "evidence_count",
                "study_locus_count",
            ]
        ]
        .drop_duplicates("disease_id")
        .sort_values(["ot_score", "disease_name"], ascending=[False, True], kind="mergesort")
    )
    variants = (
        gene_rows[
            ["chrom", "chrom_sort", "pos", "varid", "AF_sas", "AF_nfe", "cadd_phred", "sas_enrichment", "consequence"]
        ]
        .drop_duplicates(["chrom", "pos", "varid", "AF_sas", "AF_nfe", "cadd_phred", "sas_enrichment", "consequence"])
        .sort_values(["chrom_sort", "pos", "varid"], kind="mergesort")
        .drop(columns=["chrom_sort"])
    )
    variants["AF_sas"] = variants["AF_sas"].map(lambda value: round_float(value, 8))
    variants["AF_nfe"] = variants["AF_nfe"].map(lambda value: round_float(value, 8))
    variants["cadd_phred"] = variants["cadd_phred"].map(lambda value: round_float(value, 3))
    variants["sas_enrichment"] = variants["sas_enrichment"].map(lambda value: round_float(value, 3))

    label = f"{gene_symbol} ({ensembl_gene_id})" if gene_symbol else ensembl_gene_id
    return label, diseases.to_dict("records"), variants.to_dict("records")


if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=8050)
