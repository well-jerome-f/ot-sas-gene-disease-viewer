# Open Targets SAS Gene-Disease Variant Viewer

This folder contains a small local workflow for browsing SAS-enriched coding variants against Open Targets GWAS credible-set target-disease evidence.

The default build keeps disease-gene associations with Open Targets/L2G score `>= 0.25`. Genes with at least one mapped LoF variant are highlighted in blue in the viewer.

## Build the mapped file

```bash
python3 /Users/jerome/work_data/codex_workspace/ot_sas_gene_disease_viewer/build_mapping.py
```

Default output:

```text
/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.score_ge_0.25.parquet
```

The disease filter keeps terms assigned to disease therapeutic areas and excludes broad non-disease areas:

```text
EFO_0000651 phenotype
EFO_0001444 measurement
EFO_0002571 medical procedure
EFO_0005932 animal disease
GO_0008150 biological_process
```

To also write a TSV copy:

```bash
python3 /Users/jerome/work_data/codex_workspace/ot_sas_gene_disease_viewer/build_mapping.py \
  --tsv-output /Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.tsv
```

To change the score cutoff:

```bash
python3 /Users/jerome/work_data/codex_workspace/ot_sas_gene_disease_viewer/build_mapping.py --min-score 0.5
```

## Start the viewer

```bash
python3 app.py
```

Then open:

```text
http://127.0.0.1:8050
```

Clicking a gene row updates the associated disease table, sorted by score, and the mapped coding variant table, sorted by chromosome and position.

## Data Access

While the viewer is running, the mapped parquet is available through:

```text
http://127.0.0.1:8050/download/mapping.parquet
```

The local file is:

```text
/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.score_ge_0.25.parquet
```

Build statistics are written beside the parquet as:

```text
/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.score_ge_0.25.parquet.stats.json
```

## Deploy

This is a Dash/Flask app, so it needs a Python web host. GitHub Pages cannot run it directly because GitHub Pages only serves static files.

The repo is ready for a service such as Render, Railway, Fly.io, or a VM. The production start command is:

```bash
gunicorn app:server --bind 0.0.0.0:$PORT
```

The app reads data from:

```text
data/ot_sas_disease_gene_variant_map.score_ge_0.25.parquet
```

You can override that path with:

```bash
OT_SAS_MAPPING_PATH=/path/to/mapping.parquet
```

### Render Example

1. Push this repository to GitHub.
2. In Render, create a new Web Service from the GitHub repository.
3. Use Python as the runtime.
4. Use this build command:

```bash
pip install -r requirements.txt
```

5. Use this start command:

```bash
gunicorn app:server --bind 0.0.0.0:$PORT
```

6. Deploy and share the Render URL.
