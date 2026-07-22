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


def build_dataset(csv_path: Path, full_names: dict[str, str], class_names: dict[str, str]) -> dict:
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
    parser.add_argument("--output", type=Path, default=Path("data/ec_templates.js"))
    args = parser.parse_args()

    enzyme_dat = fetch_text(EXPASY_ENZYME_DAT_URL)
    byclass_html = fetch_text(EXPASY_BYCLASS_URL)
    dataset = build_dataset(
        args.csv,
        full_names=parse_enzyme_dat(enzyme_dat),
        class_names=parse_byclass_names(byclass_html),
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
