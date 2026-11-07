#!/usr/bin/env python3
"""Validate static data and local assets without a browser runtime."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


def load_data() -> dict:
    text = Path("data/ec_templates.js").read_text(encoding="utf-8")
    prefix = "window.EC_TEMPLATE_DATA = "
    assert text.startswith(prefix), "data/ec_templates.js does not define EC_TEMPLATE_DATA"
    return json.loads(text.removeprefix(prefix).rstrip(";\n"))


def main() -> None:
    data = load_data()
    revision = hashlib.sha256(Path("data/ec_templates.js").read_bytes()).hexdigest()[:16]
    html_text = Path("index.html").read_text(encoding="utf-8")
    assert f'src="data/ec_templates.js?v={revision}"' in html_text
    assert f'src="app.js?v={revision}"' in html_text
    assert data["metadata"]["sourceCsvRows"] == 22326
    assert data["metadata"]["positiveExampleCount"] == 22326
    assert data["metadata"]["sourceExampleCounts"] == {"brenda": 22326}
    assert data["metadata"]["uniqueReactionCount"] == 14891
    assert data["metadata"]["uniqueEnzymeCount"] == 6398
    assert data["metadata"]["uniqueProteinSequenceCount"] == 6392
    assert data["metadata"]["uniqueEcCount"] == len(data["ecEntries"]) == 2765
    assert data["metadata"]["uniqueTemplateCount"] == 4377
    assert data["metadata"]["lookupEntryCount"] == len(data["templates"]) == 7979
    assert data["metadata"]["ecAnnotationCount"] == 2681
    assert data["metadata"]["ecTemplateAssociationCount"] == 7793
    assert data["metadata"]["examplesWithEcCount"] == 21879
    assert data["metadata"]["templatesWithEcCount"] == 4327
    assert data["metadata"]["templatesWithExamples"] == 7979
    assert data["metadata"]["examplesWithUniprotCount"] == 22326
    assert data["metadata"]["examplesWithFullReactionCount"] == 22326
    assert all(entry["name"] for entry in data["ecEntries"])
    assert all(entry["ec"] != "[]" for entry in data["ecEntries"])
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
        "[]" not in example["uniprotIds"]
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
    assert len(examples_by_id) == 22326
    standard_entries = [
        item for item in data["templates"]
        if re.fullmatch(data["metadata"]["standardEcPattern"], item["ec"])
    ]
    assert len({item["templateId"] for item in data["templates"]}) == 4377
    assert len({item["ec"] for item in standard_entries}) == 2681
    assert len({(item["ec"], item["templateId"]) for item in standard_entries}) == 7793
    assert len({item["templateId"] for item in standard_entries}) == 4327
    assert len({example["id"] for item in standard_entries for example in item["examples"]}) == 21879
    brenda_example = examples_by_id["aries-brenda-0005a738cd26e819ab8b"]
    assert brenda_example["sourceReactionIds"] == ["15890"]
    assert brenda_example["sourceEnzymeIds"] == ["9106"]
    assert brenda_example["uniprotIds"] == ["B3VA58"]
    assert brenda_example["fullReactionSmiles"]
    assert all(example["dataset"] == "brenda" for example in examples_by_id.values())
    assert all(example_id.startswith("aries-brenda-") for example_id in examples_by_id)
    assert data["metadata"]["datasetScope"] == "BRENDA only"
    assert data["metadata"]["releaseId"] == "aries-db-2026-10-31-grover-fixed4-v2"

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
        f"{data['metadata']['uniqueTemplateCount']} distinct templates,",
        f"{len(data['templates'])} lookup entries",
    )


if __name__ == "__main__":
    main()
