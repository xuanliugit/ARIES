# ARIES

<div align="center">

[![Code](https://img.shields.io/badge/Code-xuanliugit/ARIES-blue?logoWidth=40)](https://github.com/xuanliugit/ARIES)
[![Dataset](https://img.shields.io/badge/Dataset-Hugging_Face-yellow?logoWidth=40)](https://huggingface.co/datasets/xuan-liu/ARIES-data/)
[![Model Weights](https://img.shields.io/badge/Model_Weights-Hugging_Face-yellow?logoWidth=40)](https://huggingface.co/xuan-liu/ARIES)

</div>

ARIES (atom-resolved regioselectivity inference from enzymes and symbolic rules) predicts enzyme-conditioned reaction sites on substrate atoms.

## Installation

Create and activate the dedicated conda environment, then install ARIES from
the repository root:

```bash
conda create -n aries python=3.12 -y
conda activate aries
python -m pip install --upgrade pip
python -m pip install -e ./aries
```

To convert raw reaction SMILES into ARIES inputs, install the reaction
preprocessing components in the same `aries` environment:

```bash
python -m pip install -e ./aries_rxn
python -m pip install -e "./aries[process]"
```

For development and testing:

```bash
conda activate aries
python -m pip install -e "./aries[test]"
pytest aries/tests
```

## Model weights

[Download ARIES model weights](https://huggingface.co/xuan-liu/ARIES) into
`model_weights/` and retain the directory names.

```bash
hf download xuan-liu/ARIES --local-dir model_weights
```

See [`model_weights/README.md`](model_weights/README.md) for the layout and
base/adapter compatibility rules.

## Data

Download the compact inference and evaluation datasets, exact GROVER LMDBs,
ESM-C 6B residue LMDBs, ground-truth atom labels, and canonical split files
from [`xuan-liu/ARIES-data`](https://huggingface.co/datasets/xuan-liu/ARIES-data):

```bash
hf download xuan-liu/ARIES-data \
  --repo-type dataset \
  --local-dir ARIES-data
```

BRENDA and Halogenase provide six protocols with seeds 42-44: random, unknown
enzyme, unknown substrate, unknown transformation, unknown transformation plus
substrate, and unknown transformation plus substrate plus enzyme. P450,
glycosyltransferase, methyltransferase, fatty-acid-desaturase, and
prenyltransferase provide random, unknown-enzyme, unknown-substrate, and joint
unknown-enzyme-plus-substrate splits with seeds 42-44. Released adapters use
the random splits.

Each released directory also contains an independent `reaction.csv` with the
original reaction SMILES, enzyme name/accession/sequence, and EC number. The EC
number is metadata and is not used as a model feature.

## Reaction preprocessing

Convert raw enzymatic reactions into model-ready substrate/template rows and
atom labels with:

```bash
conda activate aries
aries-data-process \
  --input-csv reaction.csv \
  --output-dir processed \
  --device cuda
```

`reaction.csv` uses the columns `example_id`, `reaction_smiles`, `enzyme_name`,
`enzyme_accession`, `enzyme_sequence`, and `ec_number`; multiple EC annotations
are separated by semicolons. The command writes
`inputs.csv` and `atom_labels.csv`. It performs reaction mapping, substrate and
product selection, ARIES template extraction, and atom-label generation; it
does not compute GROVER or ESM-C embeddings. Run
`aries-precompute-embeddings` on the resulting `inputs.csv` as the next step.

## Inference

ARIES inference reads dataset-style LMDB embeddings. For new examples, prepare
a label-free CSV with these required fields:

| Column | Description |
|---|---|
| `example_id` | Unique, non-empty example identifier. |
| `enzyme_sequence` | Protein sequence; whitespace is removed and letters are normalized to uppercase. |
| `canonical_substrate_smiles` | Canonical isomeric substrate SMILES. |
| `mapped_substrate_smiles` | Exact mapped SMILES whose atom order defines prediction indices. |
| `site_template_smarts` | Reaction-site template used by ARIES. |

`eligible_atom_indices` and arbitrary provenance columns are optional. Generate
the canonical LMDB files with an explicit local GROVER checkpoint and source
checkout; the GROVER encoder checkpoint is not distributed with ARIES:

```bash
conda activate aries
aries-precompute-embeddings \
  --input-csv inputs.csv \
  --output-dir embeddings \
  --grover-checkpoint /path/to/grover_large.pt \
  --grover-code-dir /path/to/grover \
  --esm-model biohub/ESMC-6B \
  --device cuda
```

This command deduplicates substrates and sequences, then runs GROVER and ESM-C
sequentially so both encoders do not occupy GPU memory together. It writes only
`grover_exact_aries_fg.lmdb`, `esmc6b_residue_embeddings.lmdb`,
`esmc6b_embedding_index.csv`, and `embedding_summary.json`. Use `--resume` to
validate and skip complete entries, and `--local-files-only` when the ESM-C
model is already downloaded.

The bundled 21-case experimental example already contains the canonical LMDB
and index files. Run base-model inference directly with:

```bash
conda activate aries
aries-inference \
  --input-csv data/example_data/inputs.csv \
  --embedding-dir data/example_data \
  --weights-dir model_weights/aries-xl-brenda \
  --output-csv predictions_brenda.csv \
  --device auto
```

Apply the Halogenase adapter to its required ARIES-XL-Brenda base:

```bash
aries-inference \
  --input-csv data/example_data/inputs.csv \
  --embedding-dir data/example_data \
  --weights-dir model_weights/aries-xl-brenda \
  --lora-weights-dir model_weights/aries-xl-hal-lora \
  --output-csv predictions_hal_specialist.csv \
  --device auto
```

Apply a specialist-hub adapter to its required ARIES-XL random base, for
example P450:

```bash
aries-inference \
  --input-csv data/example_data/inputs.csv \
  --embedding-dir data/example_data \
  --weights-dir model_weights/aries-random-xl \
  --lora-weights-dir model_weights/specialist_hub/p450-random \
  --output-csv predictions_p450_specialist.csv \
  --device auto
```

There is no model-name argument: `--weights-dir` selects the base model, and
`--lora-weights-dir` optionally selects a compatible adapter. ARIES loads
seeds 42-44 sequentially, averages their atom probabilities, and writes both
all-atom and eligible-atom ranks. Use `--device cpu`, `--device cuda`, or
`--device auto`.

The example input and its accepted legacy column aliases are documented in
[`data/example_data/README.md`](data/example_data/README.md). GROVER and ESM-C
LMDB embeddings must be precomputed before inference.

The Hugging Face release uses compact CSVs without provenance, audit, or
selection-reason columns. Ground truth is stored separately in
`atom_labels.csv`, so the same `examples.csv` can be passed directly to
`aries-inference`.

