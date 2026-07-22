#!/usr/bin/env python3
"""Build the static EC reaction-template lookup dataset."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen


EXPASY_ENZYME_DAT_URL = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
EXPASY_BYCLASS_URL = "https://enzyme.expasy.org/enzyme-byclass.html"

SOURCE_CSV = (
    Path(__file__).resolve().parents[2]
    / "ezspecificity"
    / "analysis"
    / "2026-07-19-aries-db-build-methods"
    / "results"
    / "ec_number_reaction_template_rows.csv"
)

REFERENCE_EXAMPLES_CSV = (
    Path(__file__).resolve().parents[2]
    / "ezspecificity"
    / "analysis"
    / "2026-07-16-clean-brenda-hal-model1-split"
    / "data"
    / "all_model_ready_examples.csv"
)

REACTION_SOURCE_CSVS = [
    (
        "brenda_pool",
        Path(__file__).resolve().parents[2]
        / "ezspecificity"
        / "data"
        / "brenda"
        / "site_selectivity_candidate_pool.csv",
    ),
    (
        "mixed_aries",
        Path(__file__).resolve().parents[2]
        / "ezspecificity"
        / "data"
        / "mixed_brenda_halogenase_corrected_fold0_20260703"
        / "aries_examples.csv",
    ),
]


def fetch_text(url: str) -> str:
    request = Request(
        url,
        headers={"User-Agent": "EC reaction template lookup data builder"},
    )
    with urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="replace")


def clean_name(value: str) -> str:
    value = re.sub(r"\s+", " ", value).strip()
    return value[:-1] if value.endswith(".") else value


def parse_enzyme_dat(text: str) -> dict[str, str]:
    names: dict[str, str] = {}
    current_id: str | None = None
    de_lines: list[str] = []

    for line in text.splitlines():
        if line.startswith("ID   "):
            current_id = line[5:].strip()
            de_lines = []
        elif line.startswith("DE   ") and current_id:
            de_lines.append(line[5:].strip())
        elif line == "//" and current_id:
            if de_lines:
                names[current_id] = clean_name(" ".join(de_lines))
            current_id = None
            de_lines = []

    return names


def normalize_class_code(raw_code: str) -> str:
    code = re.sub(r"\s+", "", raw_code.strip())
    code = code.replace("-.-", "-.-").replace(".-", ".-")
    parts = [part for part in code.split(".") if part != ""]
    parts = [part for part in parts if part != "-"]
    return ".".join(parts)


def parse_byclass_names(text: str) -> dict[str, str]:
    names: dict[str, str] = {}
    span_pattern = re.compile(
        r"HREF\s*=\s*\"/EC/[^\"]+\">"
        r"\s*<span[^>]*schema:name[^>]*>(?P<code>.*?)</span>"
        r"\s*</A>\s*<span[^>]*schema:description[^>]*>(?P<name>.*?)</span>",
        re.IGNORECASE | re.DOTALL,
    )
    for match in span_pattern.finditer(text):
        code = normalize_class_code(html.unescape(re.sub(r"<[^>]+>", "", match.group("code"))))
        name = clean_name(html.unescape(re.sub(r"<[^>]+>", "", match.group("name"))))
        if code and name:
            names[code] = name
    return names


def ec_sort_key(ec: str) -> tuple[int, tuple[int, ...], str]:
    if not re.match(r"^\d+(?:\.\d+){0,3}$", ec):
        return (1, (), ec)
    return (0, tuple(int(part) for part in ec.split(".")), ec)


def first_int(values: list[str]) -> int | None:
    parsed: list[int] = []
    for value in values:
        if value == "":
            continue
        try:
            parsed.append(int(value))
        except ValueError:
            continue
    return parsed[0] if parsed else None


def json_list(value: str) -> list:
    if not value:
        return []
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return []
    return parsed if isinstance(parsed, list) else [parsed]


def first_json_string(value: str) -> str:
    values = json_list(value)
    for item in values:
        if item is not None and str(item) != "":
            return str(item)
    return ""


def id_list(value: str) -> list[str]:
    return [str(item) for item in json_list(value) if str(item) != ""]


def clean_id(value: str) -> str:
    text = str(value or "")
    return text[:-2] if text.endswith(".0") else text


def protein_label(dataset: str, enzyme_ids: list[str], uniprots: list[str]) -> str:
    if uniprots:
        return f"UniProt {uniprots[0]}"
    enzyme_id = enzyme_ids[0] if enzyme_ids else ""
    dataset_label = {"brenda": "BRENDA", "halogenase": "Halogenase"}.get(dataset.lower(), dataset.title())
    if enzyme_id:
        return f"{dataset_label} enzyme {enzyme_id}"
    return f"{dataset_label} protein"


def parse_mapped_reaction_json(value: str, reaction_id: str = "") -> str:
    for item in json_list(value):
        if not isinstance(item, dict):
            continue
        mapped = item.get("mapped_reaction_smiles") or item.get("reaction_smiles") or ""
        if not mapped:
            continue
        if not reaction_id or str(item.get("positive_reaction", "")) == str(reaction_id):
            return mapped
    return ""


def put_reaction(reaction_map: dict[tuple[str, str, str, str], tuple[int, str]], key: tuple[str, str, str, str], value: str, score: int) -> None:
    if not value:
        return
    old = reaction_map.get(key)
    if old is None or score > old[0]:
        reaction_map[key] = (score, value)


def load_full_reaction_map(paths: list[tuple[str, Path]]) -> dict[tuple[str, str, str, str], str]:
    reaction_map: dict[tuple[str, str, str, str], tuple[int, str]] = {}
    for source_name, path in paths:
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if source_name == "brenda_pool":
                    key = (
                        "brenda",
                        row.get("source_pair_id", ""),
                        row.get("reaction", ""),
                        row.get("enzyme", ""),
                    )
                    observed = row.get("is_observed_product") == "1"
                    positive = parse_mapped_reaction_json(row.get("positive_mapped_reaction_smiles", ""), row.get("reaction", ""))
                    generated = parse_mapped_reaction_json(row.get("generation_mapped_reaction_smiles", ""), row.get("reaction", ""))
                    put_reaction(reaction_map, key, positive, 30 if observed else 20)
                    put_reaction(reaction_map, key, generated, 10)
                elif source_name == "mixed_aries":
                    dataset = row.get("dataset") or row.get("source_dataset") or ""
                    enzyme_id = clean_id(row.get("enzyme_id", ""))
                    reaction_ids = [
                        clean_id(row.get("row_reaction_id", "")),
                        clean_id(row.get("reaction_id", "")),
                    ]
                    for reaction_id in reaction_ids:
                        if not reaction_id:
                            continue
                        key = (dataset, row.get("source_pair_id", ""), reaction_id, enzyme_id)
                        put_reaction(reaction_map, key, row.get("mapped_reaction_smiles", ""), 25)
                        put_reaction(reaction_map, key, row.get("source_reaction_smiles", ""), 5)
    return {key: value for key, (_, value) in reaction_map.items()}


def load_reference_examples(path: Path, reaction_map: dict[tuple[str, str, str, str], str]) -> dict[str, dict]:
    examples: dict[str, dict] = {}
    if not path.exists():
        return examples

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            example_id = row.get("example_id", "")
            if not example_id:
                continue
            pair_ids = id_list(row.get("source_pair_ids", "")) or [row.get("source_pair_id", "")]
            reaction_ids = id_list(row.get("source_reaction_ids", ""))
            enzyme_ids = id_list(row.get("source_enzyme_ids", ""))
            uniprots = id_list(row.get("source_uniprots", ""))
            dataset = row.get("source_dataset") or first_json_string(row.get("dataset_memberships", "")) or "source"

            full_reaction = ""
            for pair_id in pair_ids:
                for reaction_id in reaction_ids:
                    for enzyme_id in enzyme_ids:
                        full_reaction = reaction_map.get((dataset, pair_id, reaction_id, enzyme_id), "")
                        if full_reaction:
                            break
                    if full_reaction:
                        break
                if full_reaction:
                    break

            canonical_substrate = row.get("canonical_substrate_smiles", "")
            mapped_substrate = row.get("mapped_substrate_smiles", "")
            if mapped_substrate == canonical_substrate:
                mapped_substrate = ""

            examples[example_id] = {
                "id": example_id,
                "proteinName": protein_label(dataset, enzyme_ids, uniprots),
                "proteinSequence": row.get("enzyme_sequence", ""),
                "dataset": dataset,
                "sourcePairIds": pair_ids,
                "sourceReactionIds": reaction_ids,
                "sourceEnzymeIds": enzyme_ids,
                "mappedSubstrateSmiles": mapped_substrate,
                "canonicalSubstrateSmiles": canonical_substrate,
                "fullReactionSmiles": full_reaction,
            }
    return examples


def numeric_parent_key(ec: str) -> str:
    parts: list[str] = []
    for part in ec.strip(".").split("."):
        if part.isdigit():
            parts.append(part)
        else:
            break
    return ".".join(parts)


def ec_name(ec: str, full_names: dict[str, str], class_names: dict[str, str]) -> tuple[str, str]:
    if ec in full_names:
        return full_names[ec], "expasy-entry"
    if ec in class_names:
        return class_names[ec], "expasy-class"
    parent = numeric_parent_key(ec)
    if parent in class_names:
        return f"Parent class: {class_names[parent]}", "expasy-parent-class"
    if parent in full_names:
        return f"Parent entry: {full_names[parent]}", "expasy-parent-entry"
    if ec == "Other.":
        return "Other / unclassified EC bucket", "dataset-bucket"
    return "Unlisted EC bucket", "dataset-bucket"


def build_dataset(csv_path: Path, full_names: dict[str, str], class_names: dict[str, str], example_details: dict[str, dict]) -> dict:
    rows = list(csv.DictReader(csv_path.open(newline="")))
    grouped: dict[tuple[str, str, str, str], dict] = {}

    for row in rows:
        key = (
            row["ec_number"],
            row["template_id"],
            row["radius0_template_smarts"],
            row["site_template_smarts"],
        )
        entry = grouped.setdefault(
            key,
            {
                "ec": row["ec_number"],
                "templateId": row["template_id"],
                "radius0Smarts": row["radius0_template_smarts"],
                "siteSmarts": row["site_template_smarts"],
                "sourceDatasets": Counter(),
                "exampleIds": [],
                "sourceReactionIds": [],
                "sourceEnzymeIds": [],
                "sourcePairIds": [],
                "examples": [],
                "selectivityIssueCount": 0,
                "rowCount": 0,
                "legalSiteCounts": [],
                "numAtoms": [],
                "positiveAtomCounts": [],
                "legalAtomCounts": [],
            },
        )
        entry["rowCount"] += 1
        entry["sourceDatasets"][row["source_dataset"]] += 1
        if row["selectivity_issue"] == "1":
            entry["selectivityIssueCount"] += 1
        for field, target in [
            ("example_id", "exampleIds"),
            ("source_reaction_id", "sourceReactionIds"),
            ("source_enzyme_id", "sourceEnzymeIds"),
            ("source_pair_id", "sourcePairIds"),
        ]:
            value = row[field]
            if value and value not in entry[target] and len(entry[target]) < 8:
                entry[target].append(value)
        example = example_details.get(row["example_id"])
        if example and all(existing["id"] != example["id"] for existing in entry["examples"]):
            entry["examples"].append(example)
        for field, target in [
            ("legal_site_count", "legalSiteCounts"),
            ("num_atoms", "numAtoms"),
            ("positive_atom_count", "positiveAtomCounts"),
            ("legal_atom_count", "legalAtomCounts"),
        ]:
            entry[target].append(row[field])

    templates = []
    for entry in grouped.values():
        entry["sourceDatasets"] = dict(sorted(entry["sourceDatasets"].items()))
        entry["legalSiteCount"] = first_int(entry.pop("legalSiteCounts"))
        entry["numAtoms"] = first_int(entry.pop("numAtoms"))
        entry["positiveAtomCount"] = first_int(entry.pop("positiveAtomCounts"))
        entry["legalAtomCount"] = first_int(entry.pop("legalAtomCounts"))
        entry["examples"].sort(key=lambda example: (0 if example.get("fullReactionSmiles") else 1, example["id"]))
        templates.append(entry)

    templates.sort(key=lambda item: (ec_sort_key(item["ec"]), item["templateId"]))

    ec_counts = Counter(row["ec_number"] for row in rows)
    template_counts = Counter(item["ec"] for item in templates)
    ec_numbers = sorted(ec_counts, key=ec_sort_key)

    ec_entries = []
    prefix_set: set[str] = set()
    for ec in ec_numbers:
        parts = ec.split(".")
        numeric_parts: list[str] = []
        for part in parts:
            if not part.isdigit():
                break
            numeric_parts.append(part)
        for i in range(1, len(numeric_parts) + 1):
            prefix_set.add(".".join(numeric_parts[:i]))
        name, name_source = ec_name(ec, full_names, class_names)
        ec_entries.append(
            {
                "ec": ec,
                "name": name,
                "nameSource": name_source,
                "rowCount": ec_counts[ec],
                "templateCount": template_counts[ec],
            }
        )

    prefixes = {}
    for prefix in sorted(prefix_set, key=ec_sort_key):
        full_matches = [ec for ec in ec_numbers if ec.startswith(prefix + ".") or ec == prefix]
        prefixes[prefix] = {
            "name": class_names.get(prefix) or full_names.get(prefix) or "",
            "rowCount": sum(ec_counts[ec] for ec in full_matches),
            "templateCount": sum(template_counts[ec] for ec in full_matches),
            "ecCount": len(full_matches),
        }

    return {
        "metadata": {
            "builtAt": datetime.now(timezone.utc).isoformat(),
            "sourceCsv": str(csv_path),
            "sourceCsvRows": len(rows),
            "uniqueEcCount": len(ec_numbers),
            "uniqueTemplateCount": len(templates),
            "referenceExamplesCsv": str(REFERENCE_EXAMPLES_CSV),
            "referenceExampleCount": len(example_details),
            "templatesWithExamples": sum(1 for item in templates if item["examples"]),
            "examplesWithFullReactionCount": sum(
                1 for item in templates for example in item["examples"] if example.get("fullReactionSmiles")
            ),
            "expasyEnzymeDatUrl": EXPASY_ENZYME_DAT_URL,
            "expasyByClassUrl": EXPASY_BYCLASS_URL,
        },
        "ecEntries": ec_entries,
        "prefixes": prefixes,
        "templates": templates,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, default=SOURCE_CSV)
    parser.add_argument("--reference-examples", type=Path, default=REFERENCE_EXAMPLES_CSV)
    parser.add_argument("--output", type=Path, default=Path("data/ec_templates.js"))
    args = parser.parse_args()

    enzyme_dat = fetch_text(EXPASY_ENZYME_DAT_URL)
    byclass_html = fetch_text(EXPASY_BYCLASS_URL)
    reaction_map = load_full_reaction_map(REACTION_SOURCE_CSVS)
    example_details = load_reference_examples(args.reference_examples, reaction_map)
    dataset = build_dataset(
        args.csv,
        full_names=parse_enzyme_dat(enzyme_dat),
        class_names=parse_byclass_names(byclass_html),
        example_details=example_details,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(dataset, ensure_ascii=False, separators=(",", ":"))
    args.output.write_text(
        "window.EC_TEMPLATE_DATA = " + payload + ";\n",
        encoding="utf-8",
    )
    print(json.dumps(dataset["metadata"], indent=2))


if __name__ == "__main__":
    main()
