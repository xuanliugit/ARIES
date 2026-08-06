#!/usr/bin/env python3
"""Validate static data and local assets without a browser runtime."""

from __future__ import annotations

import json
from pathlib import Path


def load_data() -> dict:
    text = Path("data/ec_templates.js").read_text(encoding="utf-8")
    prefix = "window.EC_TEMPLATE_DATA = "
    assert text.startswith(prefix), "data/ec_templates.js does not define EC_TEMPLATE_DATA"
    return json.loads(text.removeprefix(prefix).rstrip(";\n"))


def main() -> None:
    data = load_data()
    assert data["metadata"]["sourceCsvRows"] == 21131
    assert data["metadata"]["positiveExampleCount"] == 21131
    assert data["metadata"]["sourceExampleCounts"] == {"brenda": 19891, "halogenase": 1240}
    assert data["metadata"]["uniqueReactionCount"] == 13532
    assert data["metadata"]["uniqueEnzymeCount"] == 6180
    assert data["metadata"]["uniqueProteinSequenceCount"] == 6159
    assert data["metadata"]["uniqueEcCount"] == len(data["ecEntries"]) == 2578
    assert data["metadata"]["uniqueTemplateCount"] == len(data["templates"]) == 7313
    assert data["metadata"]["templatesWithExamples"] == 7313
    assert data["metadata"]["examplesWithUniprotCount"] == 21131
    assert data["metadata"]["examplesWithFullReactionCount"] == 21131
    assert all(entry["name"] for entry in data["ecEntries"])
    assert data["prefixes"]["1"]["name"] == "Oxidoreductases"
    assert data["prefixes"]["1.1"]["name"].startswith("Acting on the CH-OH")
    assert data["prefixes"]["1.1.1"]["name"].startswith("With NAD")
    assert any("[c;D2:1]>>[c;D3:1]-[Cl;D1:2]" == item["siteSmarts"] for item in data["templates"])
    assert any(" || " in item["siteSmarts"] for item in data["templates"])
    assert any("alcohol dehydrogenase" == item["name"] for item in data["ecEntries"])
    assert all(item["examples"] for item in data["templates"])
    assert any(item["examples"][0].get("fullReactionSmiles") for item in data["templates"])
    assert all(
        "uniprotIds" in example and example.get("proteinSequence")
        for item in data["templates"]
        for example in item["examples"]
    )
    assert all(
        key not in item
        for item in data["templates"]
        for key in ["potentialSiteCount", "numAtoms", "positiveAtomCount", "potentialAtomCount"]
    )
    assert all(
        all(
            key in example
            for key in [
                "potentialSiteCount",
                "positiveSiteCount",
                "numAtoms",
                "positiveAtomCount",
                "potentialAtomCount",
                "observedGroundTruthAtomCount",
            ]
        )
        for item in data["templates"]
        for example in item["examples"]
    )
    examples_by_id = {
        example["id"]: example
        for item in data["templates"]
        for example in item["examples"]
    }
    template_selectivity_issue_count = sum(item.get("selectivityIssueCount", 0) for item in data["templates"])
    example_selectivity_issue_count = sum(
        1
        for item in data["templates"]
        for example in item["examples"]
        if example.get("selectivityIssue")
    )
    assert template_selectivity_issue_count == example_selectivity_issue_count > 0
    assert len(examples_by_id) == 21131
    brenda_example = examples_by_id["aries-brenda-19581e27de7ced00ff1c"]
    assert brenda_example["sourceReactionIds"] == ["39253"]
    assert brenda_example["sourceEnzymeIds"] == ["25201"]
    assert brenda_example["uniprotIds"] == ["A0A3Q0KMZ9"]
    assert brenda_example["fullReactionSmiles"]
    halogenase_example = examples_by_id["aries-halogenase-4e07408562bedb8b60ce"]
    assert halogenase_example["sourceReactionIds"] == ["288"]
    assert halogenase_example["sourceEnzymeIds"] == ["201"]
    assert halogenase_example["uniprotIds"] == ["A0A0Q0FAF2"]
    assert halogenase_example["fullReactionSmiles"]

    for path in [
        "index.html",
        "styles.css",
        "app.js",
        "data/ec_templates.js",
        "vendor/RDKit_minimal.js",
        "vendor/RDKit_minimal.wasm",
    ]:
        asset = Path(path)
        assert asset.exists() and asset.stat().st_size > 0, f"missing asset: {path}"

    print(
        "static ok:",
        f"{len(data['ecEntries'])} EC buckets,",
        f"{len(data['templates'])} templates",
    )


if __name__ == "__main__":
    main()
