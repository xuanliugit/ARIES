# ARIES Lookup: Connect enzyme function to executable biochemical transformation

Static client-side lookup page for EC-number ARIES FG reaction templates from
the July 29 BRENDA/halogenase rebuild:

- `/home/n-z/xliu254/ezspecificity/data/aries_brenda/examples.csv`
- `/home/n-z/xliu254/ezspecificity/data/aries_halogenase/examples.csv`

Full mapped reactions come from the matching `reaction_details.csv` in each
dataset directory. These locations and the output contract are documented in
`analysis/2026-07-29-aries-rxn-fg-brenda-halogenase-data/README.md`.

Generated data contains 2,343 EC buckets and 6,696 deduplicated EC/template
SMARTS entries from 18,588 accepted enzyme–reaction examples: 18,009 BRENDA
and 579 halogenase. They cover 12,432 source reactions, 5,593 source enzyme
records, and 5,573 unique protein sequences. Every example includes its updated
EC number, enzyme sequence/identifiers, and full mapped reaction.

## Data Sources

- Full EC recommended names: https://ftp.expasy.org/databases/enzyme/enzyme.dat
- EC class/subclass labels: https://enzyme.expasy.org/enzyme-byclass.html
- Browser rendering: vendored RDKit.js assets from `@rdkit/rdkit@2025.3.4-1.0.0`.
- Example sequence/source metadata: `data/aries_brenda/examples.csv` and
  `data/aries_halogenase/examples.csv`
- Full mapped reactions: `data/aries_brenda/reaction_details.csv` and
  `data/aries_halogenase/reaction_details.csv`

## Rebuild

From this branch worktree:

```bash
conda activate ezsp
python scripts/build_data.py
```

If ExPASy snapshots have already been downloaded, the build can run without
network access:

```bash
python scripts/build_data.py \
  --enzyme-dat /path/to/enzyme.dat \
  --byclass-html /path/to/enzyme-byclass.html
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
