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
    assert data["metadata"]["sourceCsvRows"] == 18171
    assert data["metadata"]["uniqueEcCount"] == len(data["ecEntries"]) == 2464
    assert data["metadata"]["uniqueTemplateCount"] == len(data["templates"]) == 6333
    assert data["metadata"]["referenceExampleCount"] == 18171
    assert data["metadata"]["templatesWithExamples"] == 6333
    assert data["metadata"]["examplesWithFullReactionCount"] >= 11000
    assert all(entry["name"] for entry in data["ecEntries"])
    assert data["prefixes"]["1"]["name"] == "Oxidoreductases"
    assert data["prefixes"]["1.1"]["name"].startswith("Acting on the CH-OH")
    assert data["prefixes"]["1.1.1"]["name"].startswith("With NAD")
    assert any("[O&H1&+0&D1:1]>>" in item["siteSmarts"] for item in data["templates"])
    assert any("alcohol dehydrogenase" == item["name"] for item in data["ecEntries"])
    assert all(item["examples"] for item in data["templates"])
    assert any(item["examples"][0].get("fullReactionSmiles") for item in data["templates"])
    assert all(example.get("proteinName") and example.get("proteinSequence") for item in data["templates"] for example in item["examples"])

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
