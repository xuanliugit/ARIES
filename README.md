# ARIES Lookup: Connect enzyme function to executable biochemical transformation

Static client-side lookup page for EC-number reaction templates from:

`/home/n-z/xliu254/ezspecificity/analysis/2026-07-19-aries-db-build-methods/results/ec_number_reaction_template_rows.csv`

Generated data contains 2,464 EC buckets and 6,333 deduplicated EC/template SMARTS entries from 18,171 source rows.
Each template includes model-ready examples from `analysis/2026-07-16-clean-brenda-hal-model1-split/data/all_model_ready_examples.csv`; the page shows one example by default and expands the rest on demand.

## Data Sources

- Full EC recommended names: https://ftp.expasy.org/databases/enzyme/enzyme.dat
- EC class/subclass labels: https://enzyme.expasy.org/enzyme-byclass.html
- Browser rendering: vendored RDKit.js assets from `@rdkit/rdkit@2025.3.4-1.0.0`.
- Example sequence/source metadata: `analysis/2026-07-16-clean-brenda-hal-model1-split/data/all_model_ready_examples.csv`

## Rebuild

From this branch worktree:

```bash
conda activate ezsp
python scripts/build_data.py
```

To use another CSV:

```bash
python scripts/build_data.py --csv /path/to/ec_number_reaction_template_rows.csv
```

## Serve Locally

```bash
python -m http.server 4173
```

Then open `http://localhost:4173/`.

## Validate

The data/static-server checks require only Python:

```bash
conda activate ezsp
python -m py_compile scripts/build_data.py scripts/validate_static.py
python scripts/validate_static.py
```

`scripts/validate_page.py` contains an optional Playwright browser check for EC drilldown, SMARTS search, copy buttons, and RDKit SVG rendering.
