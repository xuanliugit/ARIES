# ARIES Lookup: Connect enzyme function to executable biochemical transformation

Static client-side lookup page for EC-number ARIES FG reaction templates from
the current BRENDA/halogenase exports, verified on October 2, 2026:

- `/home/n-z/xliu254/ezspecificity/data/aries_brenda/examples.csv`
- `/home/n-z/xliu254/ezspecificity/data/aries_halogenase/examples.csv`

Full mapped reactions come from the matching `reaction_details.csv` in each
dataset directory. These locations and the output contract are documented in
`analysis/2026-07-29-aries-rxn-fg-brenda-halogenase-data/README.md`.

Generated data contains 22,920 accepted enzyme–reaction examples: 22,326 BRENDA
and 594 halogenase, with 4,409 distinct template IDs. They cover 15,178 source
reactions, 6,701 source enzyme records, and 6,680 unique protein sequences.
Every example includes its source annotations, enzyme sequence/identifiers,
and full mapped reaction.

The manuscript count accepts four dot-separated EC fields: a numeric first
field, then numeric or `-` fields. Under that same rule, 22,039 examples carry
2,690 distinct EC annotations and 7,819 distinct `(EC, template_id)`
associations, covering 4,338 of the 4,409 templates. This is a syntax rule,
not a claim that every annotation is independently validated.

The lookup also preserves 84 B-suffixed source annotations and `Other.`.
Including these gives 2,775 source buckets and 8,042 lookup entries, with
881 examples outside the standard-format subset (447 B-suffixed and 434
unclassified). The old header called these 8,042 entries "templates"; they
are associations, not distinct transformations. The header now reports the
distinct-template total and the manuscript-compatible EC coverage separately.
Source file SHA256 hashes are embedded in the generated metadata.

The current example-file hashes exactly match the manuscript audit in
`analysis/2026-09-14-manuscript-results-revision/ec_association_check.json`.
A fresh audit and count comparison are saved in
`analysis/2026-10-02-lookup-ec-count-audit/` in the main worktree.
The [published count audit](data/count_audit.json) includes the verified counts,
source hashes, and breakdown by dataset.

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
