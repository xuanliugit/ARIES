from __future__ import annotations

import json

import pytest

from aries_rxn import (
    AriesTemplate,
    apply_aries_template,
    canonical_smiles,
    extract_aries_template,
)


def assert_product(template: AriesTemplate, substrate: str, expected: str) -> None:
    result = apply_aries_template(template, substrate, max_products=100)
    assert canonical_smiles(expected) in result.products, result


def test_repeat_pattern_enumerates_partial_and_complete_rewrites() -> None:
    template = extract_aries_template(
        "[NH2:1][CH2:2][CH2:3][NH2:4]",
        "[OH:1][CH2:2][CH2:3][OH:4]",
        mode="core",
    )
    assert template.status == "ok"
    assert len(template.rules) == 1
    assert template.rules[0].observed_repeats == 2
    result = apply_aries_template(template, "NCCN")
    assert result.products == {canonical_smiles("NCCO"), canonical_smiles("OCCO")}


def test_alcohol_attachment_generalizes_but_carboxylic_acid_is_blocked() -> None:
    template = extract_aries_template("[CH3:1][OH:2]", "[CH4:1]", mode="core")
    assert_product(template, "CCO", "CC")
    assert_product(template, "CC(C)(C)O", "CC(C)C")
    assert not apply_aries_template(template, "C(=O)O").products


def test_fg_keeps_complete_guanidine_group() -> None:
    template = extract_aries_template(
        "[NH2:1][C:2](=[NH:3])[NH2:4]",
        "[NH2:1][C:2](=[O:3])[NH2:4]",
        mode="fg",
    )
    assert template.status == "ok"
    assert template.rules
    reactant = template.rules[0].reactant_smarts
    assert sum(reactant.count(f":{local_id}]") for local_id in range(1, 5)) == 4


def test_fg_excludes_ring_group_expansion_by_default() -> None:
    template = extract_aries_template(
        "[cH:1]1[cH:2][cH:3][c:4]([OH:7])[cH:5][cH:6]1",
        "[cH:1]1[cH:2][cH:3][c:4]([Cl:7])[cH:5][cH:6]1",
        mode="fg",
    )
    assert template.status == "ok"
    assert template.rules
    assert len(template.rules[0].context_ids) < 6


@pytest.mark.parametrize("mode", ["core", "context", "fg"])
def test_unchanged_tetrahedral_stereo_is_not_serialized(mode: str) -> None:
    template = extract_aries_template(
        "[C@@H:1]([OH:2])([CH3:3])[NH2:4]",
        "[C@@H:1]([Cl:2])([CH3:3])[NH2:4]",
        mode=mode,
    )
    assert template.status == "ok"
    assert "@@" not in template.compact_smarts
    assert "@" not in template.compact_smarts


def test_changed_tetrahedral_stereo_is_an_explicit_action() -> None:
    template = extract_aries_template(
        "[C@@H:1]([OH:2])([CH3:3])[NH2:4]",
        "[C@H:1]([Cl:2])([CH3:3])[NH2:4]",
        mode="core",
    )
    assert template.status == "ok"
    assert any(rule.tetra_actions for rule in template.rules)


def test_json_round_trip_is_stable() -> None:
    template = extract_aries_template(
        "[CH3:1][C:2](=[O:3])[OH:4]",
        "[CH3:1][C:2](=[O:3])[Cl:4]",
        mode="fg",
    )
    encoded = template.to_json()
    restored = AriesTemplate.from_json(encoded)
    assert restored.to_json() == encoded
    assert json.loads(encoded)["mode"] == "fg"


def test_predicted_centers_use_input_atom_maps() -> None:
    template = extract_aries_template("[CH3:1][OH:2]", "[CH4:1]", mode="core")
    result = apply_aries_template(template, "[CH3:11][CH2:12][OH:13]")
    assert (12, 13) in result.predicted_centers
    assert result.product_centers[canonical_smiles("CC")] == {(12, 13)}


def test_disconnected_reactants_allow_intramolecular_coupling() -> None:
    template = extract_aries_template(
        "[NH2:1].[CH3:2][Cl:3]",
        "[NH:1][CH3:2]",
        mode="core",
    )
    result = apply_aries_template(template, "NCCCl")
    assert result.products
