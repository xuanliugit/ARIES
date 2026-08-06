#!/usr/bin/env python3
"""Build the static lookup dataset from the current ARIES FG exports."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen


EXPASY_ENZYME_DAT_URL = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
EXPASY_BYCLASS_URL = "https://enzyme.expasy.org/enzyme-byclass.html"

PROJECT_ROOT = Path(__file__).resolve().parents[2] / "ezspecificity"
DEFAULT_SOURCES = (
    (
        "brenda",
        PROJECT_ROOT / "data" / "aries_brenda" / "examples.csv",
        PROJECT_ROOT / "data" / "aries_brenda" / "reaction_details.csv",
    ),
    (
        "halogenase",
        PROJECT_ROOT / "data" / "aries_halogenase" / "examples.csv",
        PROJECT_ROOT / "data" / "aries_halogenase" / "reaction_details.csv",
    ),
)


def fetch_text(url: str) -> str:
    request = Request(url, headers={"User-Agent": "ARIES EC lookup data builder"})
    with urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="replace")


def load_text(path: Path | None, url: str) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path else fetch_text(url)


def raise_csv_field_limit() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


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
    parts = [part for part in re.sub(r"\s+", "", raw_code.strip()).split(".") if part not in {"", "-"}]
    return ".".join(parts)


def parse_byclass_names(text: str) -> dict[str, str]:
    names: dict[str, str] = {}
    pattern = re.compile(
        r'HREF\s*=\s*"/EC/[^"]+">\s*<span[^>]*schema:name[^>]*>(?P<code>.*?)</span>'
        r"\s*</A>\s*<span[^>]*schema:description[^>]*>(?P<name>.*?)</span>",
        re.IGNORECASE | re.DOTALL,
    )
    for match in pattern.finditer(text):
        code = normalize_class_code(html.unescape(re.sub(r"<[^>]+>", "", match.group("code"))))
        name = clean_name(html.unescape(re.sub(r"<[^>]+>", "", match.group("name"))))
        if code and name:
            names[code] = name
    return names


def json_list(value: str) -> list:
    if not value:
        return []
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return []
    return parsed if isinstance(parsed, list) else [parsed]


def id_list(value: str) -> list[str]:
    values = json_list(value)
    if not values and value.strip().lower() not in {"", "[]", "null", "none", "nan"}:
        values = re.split(r"[,;|\s]+", value.strip())
    result: list[str] = []
    for item in values:
        text = str(item).strip()
        if text.endswith(".0"):
            text = text[:-2]
        if text and text.lower() not in {"[]", "nan", "none", "null"} and text not in result:
            result.append(text)
    return result


def maybe_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def ec_sort_key(ec: str) -> tuple[int, tuple[int, ...], str]:
    if not re.match(r"^\d+(?:\.\d+){0,3}$", ec):
        return (1, (), ec)
    return (0, tuple(int(part) for part in ec.split(".")), ec)


def numeric_parent_key(ec: str) -> str:
    parts: list[str] = []
    for part in ec.strip(".").split("."):
        if not part.isdigit():
            break
        parts.append(part)
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


def load_reactions(dataset: str, path: Path) -> dict[tuple[str, str], str]:
    reactions: dict[tuple[str, str], str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("funnel_stage") != "11_aries_fg_chemical_funnel_accepted":
                continue
            reaction_id = row.get("reaction_index", "")
            mapped_reaction = row.get("mapped_rxn", "")
            if reaction_id and mapped_reaction:
                reactions[(dataset, reaction_id)] = mapped_reaction
    return reactions


def load_examples(sources: tuple[tuple[str, Path, Path], ...]) -> tuple[list[dict], dict]:
    reaction_map: dict[tuple[str, str], str] = {}
    for dataset, _, reaction_path in sources:
        reaction_map.update(load_reactions(dataset, reaction_path))

    rows: list[dict] = []
    source_example_counts: Counter[str] = Counter()
    for dataset, examples_path, _ in sources:
        with examples_path.open(newline="") as handle:
            for raw in csv.DictReader(handle):
                if raw.get("molecule_label") != "1" or raw.get("is_negative_pair") != "0":
                    continue
                source_example_counts[dataset] += 1
                reaction_id = str(raw.get("source_reaction_id", ""))
                ec_numbers = id_list(raw.get("source_ec_numbers", "")) or ["Other."]
                uniprots = id_list(raw.get("source_uniprots", ""))
                canonical_substrate = raw.get("canonical_substrate_smiles", "")
                mapped_substrate = raw.get("mapped_substrate_smiles", "")
                if mapped_substrate == canonical_substrate:
                    mapped_substrate = ""
                example = {
                    "id": raw.get("example_id", ""),
                    "uniprotIds": uniprots,
                    "proteinSequence": raw.get("enzyme_sequence", ""),
                    "dataset": dataset,
                    "sourcePairIds": [str(raw.get("source_pair_id", ""))],
                    "sourceReactionIds": [reaction_id],
                    "sourceEnzymeIds": [str(raw.get("source_enzyme_id", ""))],
                    "selectivityIssue": raw.get("selectivity_issue") == "1",
                    "potentialSiteCount": maybe_int(raw.get("potential_site_count", "")),
                    "positiveSiteCount": maybe_int(raw.get("positive_site_count", "")),
                    "numAtoms": maybe_int(raw.get("num_atoms", "")),
                    "positiveAtomCount": maybe_int(raw.get("positive_atom_count", "")),
                    "potentialAtomCount": maybe_int(raw.get("potential_atom_count", "")),
                    "observedGroundTruthAtomCount": maybe_int(raw.get("observed_ground_truth_atom_count", "")),
                    "mappedSubstrateSmiles": mapped_substrate,
                    "canonicalSubstrateSmiles": canonical_substrate,
                    "fullReactionSmiles": reaction_map.get((dataset, reaction_id), ""),
                }
                for ec in ec_numbers:
                    rows.append(
                        {
                            "ec": ec,
                            "templateId": raw.get("template_id", ""),
                            "siteSmarts": raw.get("aries_fg_template_smarts") or raw.get("site_template_smarts", ""),
                            "dataset": dataset,
                            "example": example,
                        }
                    )

    metadata = {
        "sourceExampleCounts": dict(sorted(source_example_counts.items())),
        "positiveExampleCount": sum(source_example_counts.values()),
        "uniqueReactionCount": len({(row["dataset"], row["example"]["sourceReactionIds"][0]) for row in rows}),
        "uniqueEnzymeCount": len({(row["dataset"], row["example"]["sourceEnzymeIds"][0]) for row in rows}),
        "uniqueProteinSequenceCount": len({row["example"]["proteinSequence"] for row in rows}),
    }
    return rows, metadata


def build_dataset(
    sources: tuple[tuple[str, Path, Path], ...],
    full_names: dict[str, str],
    class_names: dict[str, str],
) -> dict:
    rows, source_metadata = load_examples(sources)
    grouped: dict[tuple[str, str, str], dict] = {}
    for row in rows:
        key = (row["ec"], row["templateId"], row["siteSmarts"])
        entry = grouped.setdefault(
            key,
            {
                "ec": row["ec"],
                "templateId": row["templateId"],
                "siteSmarts": row["siteSmarts"],
                "sourceDatasets": Counter(),
                "exampleIds": [],
                "sourceReactionIds": [],
                "sourceEnzymeIds": [],
                "sourcePairIds": [],
                "examples": [],
                "selectivityIssueCount": 0,
                "rowCount": 0,
            },
        )
        example = row["example"]
        entry["rowCount"] += 1
        entry["sourceDatasets"][row["dataset"]] += 1
        entry["selectivityIssueCount"] += int(example["selectivityIssue"])
        for value, target in (
            (example["id"], "exampleIds"),
            (example["sourceReactionIds"][0], "sourceReactionIds"),
            (example["sourceEnzymeIds"][0], "sourceEnzymeIds"),
            (example["sourcePairIds"][0], "sourcePairIds"),
        ):
            if value and value not in entry[target] and len(entry[target]) < 8:
                entry[target].append(value)
        if all(old["id"] != example["id"] for old in entry["examples"]):
            entry["examples"].append(example)

    templates = list(grouped.values())
    for entry in templates:
        entry["sourceDatasets"] = dict(sorted(entry["sourceDatasets"].items()))
        entry["examples"].sort(key=lambda example: example["id"])
    templates.sort(key=lambda item: (ec_sort_key(item["ec"]), item["templateId"]))

    ec_counts = Counter(row["ec"] for row in rows)
    template_counts = Counter(item["ec"] for item in templates)
    ec_numbers = sorted(ec_counts, key=ec_sort_key)
    prefix_set: set[str] = set()
    ec_entries = []
    for ec in ec_numbers:
        numeric_parts: list[str] = []
        for part in ec.split("."):
            if not part.isdigit():
                break
            numeric_parts.append(part)
        for index in range(1, len(numeric_parts) + 1):
            prefix_set.add(".".join(numeric_parts[:index]))
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
        matches = [ec for ec in ec_numbers if ec == prefix or ec.startswith(prefix + ".")]
        prefixes[prefix] = {
            "name": class_names.get(prefix) or full_names.get(prefix) or "",
            "rowCount": sum(ec_counts[ec] for ec in matches),
            "templateCount": sum(template_counts[ec] for ec in matches),
            "ecCount": len(matches),
        }

    metadata = {
        "builtAt": datetime.now(timezone.utc).isoformat(),
        "sourceExamplesCsvs": [str(examples_path) for _, examples_path, _ in sources],
        "reactionDetailsCsvs": [str(reaction_path) for _, _, reaction_path in sources],
        "sourceCsvRows": len(rows),
        "uniqueEcCount": len(ec_numbers),
        "uniqueTemplateCount": len(templates),
        **source_metadata,
        "templatesWithExamples": sum(bool(item["examples"]) for item in templates),
        "examplesWithUniprotCount": sum(
            bool(example["uniprotIds"]) for item in templates for example in item["examples"]
        ),
        "examplesWithFullReactionCount": sum(
            bool(example["fullReactionSmiles"]) for item in templates for example in item["examples"]
        ),
        "expasyEnzymeDatUrl": EXPASY_ENZYME_DAT_URL,
        "expasyByClassUrl": EXPASY_BYCLASS_URL,
    }
    return {"metadata": metadata, "ecEntries": ec_entries, "prefixes": prefixes, "templates": templates}


def main() -> None:
    raise_csv_field_limit()
    parser = argparse.ArgumentParser()
    parser.add_argument("--enzyme-dat", type=Path, help="Local ENZYME enzyme.dat; downloaded when omitted")
    parser.add_argument("--byclass-html", type=Path, help="Local ExPASy EC class page; downloaded when omitted")
    parser.add_argument("--output", type=Path, default=Path("data/ec_templates.js"))
    args = parser.parse_args()

    dataset = build_dataset(
        DEFAULT_SOURCES,
        full_names=parse_enzyme_dat(load_text(args.enzyme_dat, EXPASY_ENZYME_DAT_URL)),
        class_names=parse_byclass_names(load_text(args.byclass_html, EXPASY_BYCLASS_URL)),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(dataset, ensure_ascii=False, separators=(",", ":"))
    args.output.write_text("window.EC_TEMPLATE_DATA = " + payload + ";\n", encoding="utf-8")
    print(json.dumps(dataset["metadata"], indent=2))


if __name__ == "__main__":
    main()
