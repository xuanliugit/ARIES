# ARIES Lookup: BRENDA-only ARIES-DB

The lookup serves only the canonical BRENDA collection, verified November 7, 2026.
The default builder reads `data/BRENDA/examples.csv` and the matching
`reaction_details.csv` from the ARIES repository. Halogenase and the other family
collections remain available for specialist evaluation but are outside ARIES-DB.

The page checks BRENDA source labels, example IDs and record counts before
showing lookup results. A mixed or incomplete bundle is rejected with a refresh
message. Rebuilding the export updates the data and application asset versions
in `index.html` so a refreshed page loads the matching BRENDA bundle.

The database contains **22,326 examples, 14,891 source reactions, 9,630 canonical
substrates, 6,392 exact enzyme sequences and 4,377 distinct template IDs**. Of these,
7,584 examples are issue cases and 14,742 are nonissue cases. Every exported example
retains its original ID, sequence, template SMARTS, source identifiers and full
mapped reaction.

Standard-format EC coverage comprises 21,879 examples, 2,681 EC annotations and
7,793 EC-template associations, covering 4,327 templates. Standard-format means
four dot-separated fields, with numeric first field and numeric or `-` subsequent
fields; it does not imply independent verification of EC function. Preserving
84 provisional B-suffixed source annotations gives 2,765 source buckets and 7,979
lookup associations in total. A lookup association is not a unique template.

`data/count_audit.json` records the canonical dataset release, exact input hashes,
and comparison to the independent Fig. 2 audit. Model embeddings are not used in
this symbolic lookup. The underlying ARIES model combines frozen ESM-C residue
and GROVER atom/bond features with template-context attention and an atom-level
prediction head; training results are separate from this lookup release.

Rebuild on a CPU compute node in the `ezsp` environment:

```bash
python scripts/build_data.py --project-root /home/n-z/xliu254/ezspecificity \
  --enzyme-dat /home/n-z/xliu254/ezspecificity/analysis/2026-10-02-lookup-ec-count-audit/data/enzyme.dat \
  --byclass-html /home/n-z/xliu254/ezspecificity/analysis/2026-10-02-lookup-ec-count-audit/data/enzyme-byclass.html
python scripts/validate_static.py
```

Source names use [ENZYME](https://ftp.expasy.org/databases/enzyme/enzyme.dat)
and the [EC browser](https://enzyme.expasy.org/enzyme-byclass.html). Browser
chemistry rendering uses vendored RDKit.js. Serve with `python -m http.server 4173`.
The optional `scripts/run_validate_with_server.sh` checks EC drilldown, SMARTS
search, clipboard copying and RDKit SVG rendering with Playwright.

Full provenance and Fig. 2 outputs:
`analysis/2026-11-02-brenda-only-aries-db/` in the main ARIES repository.

Validated in Chromium on an mmli compute node using the existing QIIME2 shared
libraries: `LD_LIBRARY_PATH=/home/apps/software/QIIME2/2025.4/lib:${LD_LIBRARY_PATH:-} bash scripts/run_validate_with_server.sh`. No packages were installed.
