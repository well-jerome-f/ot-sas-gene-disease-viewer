# Open Targets SAS Gene-Disease Variant Viewer

This folder contains a small local workflow for browsing SAS-enriched coding variants against Open Targets GWAS credible-set target-disease evidence.

The default build keeps disease-gene associations with Open Targets/L2G score `>= 0.15`.
It also excludes genes where NFE/FIN European populations carry rare coding variants with `0 < AF <= 0.001`, which creates the India-exclusive opportunity set.
Genes with at least one mapped LoF variant are highlighted in blue in the viewer.
Gene-level cumulative AF columns use `AF_sas` and are calculated as `1 - product(1 - AF_i)` separately for LoF and missense variants.

## Build the mapped file

```bash
python3 /Users/jerome/work_data/codex_workspace/ot_sas_gene_disease_viewer/build_mapping.py
```

Default output:

```text
/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet
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

To rebuild without the EUR rare-gene exclusion:

```bash
python3 /Users/jerome/work_data/codex_workspace/ot_sas_gene_disease_viewer/build_mapping.py --disable-eur-rare-gene-exclusion
```

## Static Viewer

The static viewer exports the filtered data to SQLite and loads it in-browser with sql.js. It can be hosted on GitHub Pages without a Python server.

Build the SQLite file:

```bash
python3 export_static_viewer.py
```

Open locally from this folder:

```bash
python3 -m http.server 8080 --directory static
```

Then open:

```text
http://127.0.0.1:8080
```

GitHub Pages should publish the `docs/` folder. GitHub's branch-source Pages UI generally supports only `/` or `/docs`, so `docs/` mirrors the static viewer files.

```text
docs/index.html
```

The browser-loaded database is:

```text
docs/data/ot_sas_viewer.sqlite
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
/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet
```

Build statistics are written beside the parquet as:

```text
/Users/jerome/work_data/codex_workspace/gnomad_sas_parquet/ot_sas_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet.stats.json
```

## Deploy

Use the static viewer for GitHub Pages. It avoids Render cold starts and does not need Python at runtime.

### GitHub Pages

1. Commit and push the repository to GitHub.
2. In GitHub, open the repository settings.
3. Go to **Pages**.
4. Set the source to deploy from a branch.
5. Select branch `main` and folder `/docs`.
6. Save. GitHub will publish a URL like:

```text
https://<user>.github.io/<repo>/
```

The static page loads:

```text
data/ot_sas_viewer.sqlite
```

### Render

The Dash/Flask app can still run on Render, Railway, Fly.io, or a VM. The production start command is:

```bash
gunicorn app:server --bind 0.0.0.0:$PORT
```

The app reads data from:

```text
data/ot_sas_disease_gene_variant_map.india_exclusive.score_ge_0.15.parquet
```

You can override that path with:

```bash
OT_SAS_MAPPING_PATH=/path/to/mapping.parquet
```

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
