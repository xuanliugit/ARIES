# ARIES RXN

ARIES RXN extracts compact reaction graph rewrites and applies them with an
RDKit-based engine. The visible template omits default hydrogen and neutral
charge constraints; the applicator restores chemically valid states.

Three extraction modes are available:

- `core`: reacting atoms and their attachment-degree constraints.
- `context`: core plus nearby heteroatom context.
- `fg`: core plus complete AccFG functional groups containing reacting atoms,
  using `AccFG(lite=True, exclude_fgs=["rings"])`.

```python
from aries_rxn import apply_aries_template, extract_aries_template

template = extract_aries_template(
    "[CH3:1][C:2](=[O:3])[OH:4]",
    "[CH3:1][C:2](=[O:3])[Cl:4]",
    mode="fg",
)
result = apply_aries_template(template, "CCC(=O)O")
print(template.compact_smarts)
print(sorted(result.products))
print(sorted(result.predicted_centers))
```

Run the test suite from this directory:

```bash
python -m pytest
```

