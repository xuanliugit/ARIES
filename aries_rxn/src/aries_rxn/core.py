#!/usr/bin/env python
"""ARIES RXN v3: compact reaction graph rewrites with explicit application semantics.

The visible rule language deliberately omits default hydrogen and neutral charge
constraints.  The applicator owns those defaults, while non-default states such
as [O-] and [nH] remain explicit.

Core and context modes are self contained apart from RDKit and NetworkX.
Functional-group mode uses AccFG's lite catalog. RDChiral is used as a
behavioral reference, not as the application engine.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import re
from dataclasses import asdict, dataclass, field, replace
from typing import Any, Iterable, Iterator, Literal, Sequence

import networkx as nx
from rdkit import Chem, RDLogger


RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")

VERSION = "ARIES_RXN/3"
Mode = Literal["core", "context", "fg"]
NodeKey = tuple[str, int]
_ACCFG_LITE = None
_INPUT_TETRA_MARKER_PROP = "_ariesInputTetraMarker"


@dataclass(frozen=True)
class AtomState:
    atomic_num: int
    aromatic: bool
    degree: int
    formal_charge: int
    radical_electrons: int
    special_h: int | None

    @classmethod
    def from_atom(cls, atom: Chem.Atom) -> "AtomState":
        return cls(
            atomic_num=int(atom.GetAtomicNum()),
            aromatic=bool(atom.GetIsAromatic()),
            degree=int(atom.GetDegree()),
            formal_charge=int(atom.GetFormalCharge()),
            radical_electrons=int(atom.GetNumRadicalElectrons()),
            special_h=special_h_count(atom),
        )


@dataclass(frozen=True)
class BondState:
    bond_type: str
    aromatic: bool
    stereo: str

    @classmethod
    def from_bond(cls, bond: Chem.Bond) -> "BondState":
        return cls(
            bond_type=str(bond.GetBondType()),
            aromatic=bool(bond.GetIsAromatic()),
            stereo=normalized_bond_stereo(bond),
        )


@dataclass(frozen=True)
class DoubleBondStereoSpec:
    begin_id: int
    end_id: int
    begin_neighbor_id: int
    end_neighbor_id: int
    stereo: str

    def text(self) -> str:
        return f"{self.begin_id}={self.end_id}:{self.stereo}"


@dataclass(frozen=True)
class SourceTetraSpec:
    atom_map: int
    action: str
    product_tag: str
    product_neighbor_maps: tuple[int, ...]


@dataclass
class AriesRule:
    reactant_smarts: str
    product_smarts: str
    center_ids: tuple[int, ...]
    attachment_ids: tuple[int, ...]
    changed_bond_ids: tuple[tuple[int, int], ...]
    reactant_double_stereo: tuple[DoubleBondStereoSpec, ...]
    product_double_stereo: tuple[DoubleBondStereoSpec, ...]
    tetra_actions: dict[int, str]
    context_ids: tuple[int, ...]
    stereo_support_ids: tuple[int, ...]
    observed_repeats: int
    source_occurrences: tuple[tuple[int, ...], ...]
    source_assignments: tuple[tuple[tuple[int, int], ...], ...]
    graph_hash: str
    reactant_states: dict[int, AtomState] = field(repr=False)
    product_states: dict[int, AtomState] = field(repr=False)

    @property
    def smarts(self) -> str:
        return f"{self.reactant_smarts}>>{self.product_smarts}"

    @property
    def repeat_text(self) -> str:
        if self.observed_repeats <= 1:
            return "1"
        return f"1..{self.observed_repeats}"

    def text(self, index: int = 1) -> str:
        centers = ",".join(str(value) for value in self.center_ids)
        stereo_parts = [
            f"r:{spec.text()}" for spec in self.reactant_double_stereo
        ] + [f"p:{spec.text()}" for spec in self.product_double_stereo]
        stereo_parts.extend(
            f"t:{local_id}:{action}"
            for local_id, action in sorted(self.tetra_actions.items())
        )
        stereo = f"; stereo={','.join(stereo_parts)}" if stereo_parts else ""
        return (
            f"R{index} {self.smarts} "
            f"{{repeat={self.repeat_text}; observed={self.observed_repeats}; "
            f"center={centers}{stereo}}}"
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "reactant_smarts": self.reactant_smarts,
            "product_smarts": self.product_smarts,
            "center_ids": list(self.center_ids),
            "attachment_ids": list(self.attachment_ids),
            "changed_bond_ids": [list(values) for values in self.changed_bond_ids],
            "reactant_double_stereo": [
                asdict(value) for value in self.reactant_double_stereo
            ],
            "product_double_stereo": [
                asdict(value) for value in self.product_double_stereo
            ],
            "tetra_actions": {
                str(key): value for key, value in self.tetra_actions.items()
            },
            "context_ids": list(self.context_ids),
            "stereo_support_ids": list(self.stereo_support_ids),
            "observed_repeats": self.observed_repeats,
            "source_occurrences": [list(values) for values in self.source_occurrences],
            "source_assignments": [
                [list(pair) for pair in assignment]
                for assignment in self.source_assignments
            ],
            "graph_hash": self.graph_hash,
            "reactant_states": {
                str(key): asdict(value) for key, value in self.reactant_states.items()
            },
            "product_states": {
                str(key): asdict(value) for key, value in self.product_states.items()
            },
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "AriesRule":
        return cls(
            reactant_smarts=str(payload["reactant_smarts"]),
            product_smarts=str(payload["product_smarts"]),
            center_ids=tuple(int(value) for value in payload["center_ids"]),
            attachment_ids=tuple(
                int(value) for value in payload.get("attachment_ids", [])
            ),
            changed_bond_ids=tuple(
                tuple(int(value) for value in pair)
                for pair in payload.get("changed_bond_ids", [])
            ),
            reactant_double_stereo=tuple(
                DoubleBondStereoSpec(**value)
                for value in payload.get("reactant_double_stereo", [])
            ),
            product_double_stereo=tuple(
                DoubleBondStereoSpec(**value)
                for value in payload.get("product_double_stereo", [])
            ),
            tetra_actions={
                int(key): str(value)
                for key, value in payload.get("tetra_actions", {}).items()
            },
            context_ids=tuple(int(value) for value in payload.get("context_ids", [])),
            stereo_support_ids=tuple(
                int(value) for value in payload.get("stereo_support_ids", [])
            ),
            observed_repeats=int(payload.get("observed_repeats", 1)),
            source_occurrences=tuple(
                tuple(int(value) for value in occurrence)
                for occurrence in payload.get("source_occurrences", [])
            ),
            source_assignments=tuple(
                tuple(
                    (int(pair[0]), int(pair[1]))
                    for pair in assignment
                )
                for assignment in payload.get("source_assignments", [])
            ),
            graph_hash=str(payload.get("graph_hash", "")),
            reactant_states={
                int(key): AtomState(**value)
                for key, value in payload.get("reactant_states", {}).items()
            },
            product_states={
                int(key): AtomState(**value)
                for key, value in payload.get("product_states", {}).items()
            },
        )


@dataclass
class AriesTemplate:
    mode: Mode
    rules: tuple[AriesRule, ...]
    status: str = "ok"
    error: str = ""
    version: str = VERSION
    source_mapped_substrate: str = ""
    source_remote_tetra_specs: tuple[SourceTetraSpec, ...] = ()

    @property
    def template_id(self) -> str:
        # Counts and identity concern only the graph rewrite itself. Display
        # metadata such as repeat, observed, center, and stereo actions is not
        # part of the compact template identity.
        payload = self.compact_smarts if self.status == "ok" else self.text
        return stable_id("arxn3", payload)

    @property
    def text(self) -> str:
        if self.status != "ok":
            return f"{self.version} {self.mode} ERROR: {self.error}"
        lines = [f"{self.version} {self.mode}"]
        lines.extend(rule.text(index) for index, rule in enumerate(self.rules, start=1))
        return "\n".join(lines)

    @property
    def compact_smarts(self) -> str:
        return " || ".join(
            canonicalize_local_atom_maps(rule.smarts)
            for rule in self.rules
        )

    @property
    def observed_rewrite_count(self) -> int:
        return sum(rule.observed_repeats for rule in self.rules)

    def as_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "mode": self.mode,
            "status": self.status,
            "error": self.error,
            "source_mapped_substrate": self.source_mapped_substrate,
            "source_remote_tetra_specs": [
                asdict(spec) for spec in self.source_remote_tetra_specs
            ],
            "rules": [rule.as_dict() for rule in self.rules],
        }

    def to_json(self, indent: int | None = None) -> str:
        return json.dumps(
            self.as_dict(),
            sort_keys=True,
            separators=(",", ":") if indent is None else None,
            indent=indent,
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "AriesTemplate":
        return cls(
            version=str(payload.get("version", VERSION)),
            mode=str(payload["mode"]),
            status=str(payload.get("status", "ok")),
            error=str(payload.get("error", "")),
            source_mapped_substrate=str(
                payload.get("source_mapped_substrate", "")
            ),
            source_remote_tetra_specs=tuple(
                SourceTetraSpec(
                    atom_map=int(value["atom_map"]),
                    action=str(value["action"]),
                    product_tag=str(value["product_tag"]),
                    product_neighbor_maps=tuple(
                        int(item)
                        for item in value.get("product_neighbor_maps", [])
                    ),
                )
                for value in payload.get("source_remote_tetra_specs", [])
            ),
            rules=tuple(AriesRule.from_dict(rule) for rule in payload.get("rules", [])),
        )

    @classmethod
    def from_json(cls, text: str) -> "AriesTemplate":
        return cls.from_dict(json.loads(text))


@dataclass(frozen=True)
class RuleMatch:
    rule_index: int
    query_to_substrate: tuple[int, ...]
    center_atom_indices: tuple[int, ...]
    center_atom_maps: tuple[int, ...]


@dataclass
class ApplicationResult:
    products: set[str] = field(default_factory=set)
    intramolecular_products: set[str] = field(default_factory=set)
    intermolecular_products: set[str] = field(default_factory=set)
    mapped_products: set[str] = field(default_factory=set)
    product_centers: dict[str, set[tuple[int, ...]]] = field(default_factory=dict)
    predicted_centers: set[tuple[int, ...]] = field(default_factory=set)
    match_count: int = 0
    attempted_combinations: int = 0
    invalid_product_count: int = 0
    invalid_product_errors: dict[str, int] = field(default_factory=dict)
    truncated: bool = False
    error: str = ""


@dataclass
class _CandidateRule:
    graph: nx.Graph
    selected_keys: set[NodeKey]
    center_keys: set[NodeKey]
    attachment_keys: set[NodeKey]
    changed_bond_keys: set[frozenset[NodeKey]]
    reactant_double_stereo: tuple[DoubleBondStereoSpec, ...]
    product_double_stereo: tuple[DoubleBondStereoSpec, ...]
    tetra_actions: dict[int, str]
    context_keys: set[NodeKey]
    stereo_keys: set[NodeKey]
    occurrence_maps: tuple[int, ...]
    source_assignment: tuple[tuple[int, int], ...]
    local_ids: dict[NodeKey, int]
    reactant_smarts: str
    product_smarts: str
    reactant_states: dict[int, AtomState]
    product_states: dict[int, AtomState]
    graph_hash: str


class _DisjointSet:
    def __init__(self) -> None:
        self.parent: dict[NodeKey, NodeKey] = {}

    def add(self, value: NodeKey) -> None:
        self.parent.setdefault(value, value)

    def find(self, value: NodeKey) -> NodeKey:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: NodeKey, right: NodeKey) -> None:
        self.add(left)
        self.add(right)
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root

    def components(self) -> list[set[NodeKey]]:
        grouped: dict[NodeKey, set[NodeKey]] = {}
        for value in self.parent:
            grouped.setdefault(self.find(value), set()).add(value)
        return list(grouped.values())


def stable_id(prefix: str, text: str) -> str:
    digest = hashlib.sha1(text.encode("utf-8")).hexdigest()[:16]
    return f"{prefix}_{digest}"


def canonicalize_local_atom_maps(smarts: str) -> str:
    """Renumber rule-local maps by first appearance on left, then right."""

    mapping: dict[int, int] = {}

    def replace_map(match: re.Match[str]) -> str:
        old_map = int(match.group(1))
        if old_map not in mapping:
            mapping[old_map] = len(mapping) + 1
        return f":{mapping[old_map]}]"

    return re.sub(r":(\d+)\]", replace_map, smarts)


def parse_mol(smiles: str, label: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        raise ValueError(f"could not parse {label}: {smiles!r}")
    input_markers = {
        int(atom_map): marker
        for token in re.findall(r"\[([^\]]+)\]", str(smiles))
        for marker_match in [re.search(r"(@@?)", token)]
        for map_match in [re.search(r":(\d+)$", token)]
        if marker_match is not None and map_match is not None
        for marker, atom_map in [
            (marker_match.group(1), map_match.group(1))
        ]
    }
    for atom in mol.GetAtoms():
        marker = input_markers.get(int(atom.GetAtomMapNum()))
        if marker:
            atom.SetProp(_INPUT_TETRA_MARKER_PROP, marker)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def map_to_atom(mol: Chem.Mol) -> dict[int, Chem.Atom]:
    result: dict[int, Chem.Atom] = {}
    for atom in mol.GetAtoms():
        atom_map = int(atom.GetAtomMapNum())
        if atom_map <= 0:
            continue
        if atom_map in result:
            raise ValueError(f"duplicate atom map number {atom_map}")
        result[atom_map] = atom
    return result


def special_h_count(atom: Chem.Atom) -> int | None:
    """Return H only when it is not safely implied by normal valence.

    `[nH]` and charged heteroatom hydrogen states are chemically meaningful.
    Ordinary alcohol/amine/carbon hydrogens remain absent from the visible rule.
    """

    total_h = int(atom.GetTotalNumHs())
    if total_h <= 0:
        return None
    if atom.GetIsAromatic() and int(atom.GetAtomicNum()) in {7, 15}:
        return total_h
    if int(atom.GetFormalCharge()) != 0:
        return total_h
    return None


def normalized_bond_stereo(bond: Chem.Bond) -> str:
    stereo = bond.GetStereo()
    if stereo in {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS}:
        return "E"
    if stereo in {Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS}:
        return "Z"
    if stereo == Chem.BondStereo.STEREOANY:
        return "ANY"
    return ""


def bond_type_name(bond: Chem.Bond | None) -> str:
    if bond is None:
        return ""
    return str(bond.GetBondType())


def atom_core_state(atom: Chem.Atom) -> tuple[Any, ...]:
    return (
        int(atom.GetAtomicNum()),
        bool(atom.GetIsAromatic()),
        int(atom.GetFormalCharge()),
        int(atom.GetNumRadicalElectrons()),
        int(atom.GetTotalNumHs()),
    )


def parity(values: Sequence[int]) -> int:
    inversions = 0
    for i, left in enumerate(values):
        for right in values[i + 1 :]:
            inversions += int(left > right)
    return inversions % 2


def tetra_relation_from_labels(
    left_tag: Chem.ChiralType,
    left_labels: Sequence[int],
    right_tag: Chem.ChiralType,
    right_labels: Sequence[int],
) -> str:
    """Compare tetrahedral orientation using aligned neighbor-role labels."""

    unspecified = Chem.ChiralType.CHI_UNSPECIFIED
    if left_tag == unspecified and right_tag == unspecified:
        return "none"
    if left_tag == unspecified:
        return "create"
    if right_tag == unspecified:
        return "erase"

    left_labels = list(left_labels)
    right_labels = list(right_labels)
    if len(left_labels) == 3:
        left_labels.append(-1)
    if len(right_labels) == 3:
        right_labels.append(-1)
    if len(left_labels) != 4 or len(right_labels) != 4:
        return "ambiguous"

    only_left = [value for value in left_labels if value not in right_labels]
    only_right = [value for value in right_labels if value not in left_labels]
    if len(only_left) > 1 or len(only_right) > 1:
        return "ambiguous"
    aligned_right = list(right_labels)
    if only_left and only_right:
        aligned_right = [
            only_left[0] if value == only_right[0] else value
            for value in aligned_right
        ]
    if set(left_labels) != set(aligned_right):
        return "ambiguous"

    left_order = {
        value: index for index, value in enumerate(sorted(left_labels))
    }
    left_parity = parity([left_order[value] for value in left_labels])
    right_parity = parity([left_order[value] for value in aligned_right])
    tags_same = left_tag == right_tag
    equivalent = (left_parity == right_parity) == tags_same
    return "retain" if equivalent else "invert"


def mapped_tetra_relation(left: Chem.Atom, right: Chem.Atom) -> str:
    """Compare mapped tetrahedral orientation independent of source notation."""

    left_tag = left.GetChiralTag()
    right_tag = right.GetChiralTag()
    unspecified = Chem.ChiralType.CHI_UNSPECIFIED
    if left_tag == unspecified and right_tag == unspecified:
        return "none"
    if left_tag == unspecified:
        return "create"
    if right_tag == unspecified:
        return "erase"

    left_labels = [
        int(atom.GetAtomMapNum()) or (1000000 + int(atom.GetIdx()))
        for atom in left.GetNeighbors()
    ]
    right_labels = [
        int(atom.GetAtomMapNum()) or (2000000 + int(atom.GetIdx()))
        for atom in right.GetNeighbors()
    ]
    return tetra_relation_from_labels(
        left_tag,
        left_labels,
        right_tag,
        right_labels,
    )


def tetra_relation(left: Chem.Atom, right: Chem.Atom) -> str:
    """Return the conservative stereochemical action used by ARIES rules.

    An inversion enters the compact template only when the mapped source
    reaction explicitly changes its ``@``/``@@`` marker. Same-marker mapped
    inversions are retained only as source-forward provenance.
    """

    relation = mapped_tetra_relation(left, right)
    if relation == "invert":
        left_marker = (
            left.GetProp(_INPUT_TETRA_MARKER_PROP)
            if left.HasProp(_INPUT_TETRA_MARKER_PROP)
            else ""
        )
        right_marker = (
            right.GetProp(_INPUT_TETRA_MARKER_PROP)
            if right.HasProp(_INPUT_TETRA_MARKER_PROP)
            else ""
        )
        if left_marker and left_marker == right_marker:
            return "retain"
    return relation


def atom_symbol(atom: Chem.Atom) -> str:
    if int(atom.GetAtomicNum()) == 1:
        return "#1"
    symbol = atom.GetSymbol()
    if atom.GetIsAromatic() and symbol in {"B", "C", "N", "O", "P", "S"}:
        return symbol.lower()
    return symbol


def charge_field(charge: int) -> str:
    if charge == 1:
        return "+"
    if charge == -1:
        return "-"
    if charge > 1:
        return f"+{charge}"
    return str(charge)


def stereo_marker(token: str) -> str:
    if "@@" in token:
        return "@@"
    if "@" in token:
        return "@"
    return ""


def atom_query_token(
    atom: Chem.Atom,
    local_id: int,
    marker: str = "",
    degree: int | None = None,
    allow_attachments: bool = False,
) -> str:
    fields = [f"{atom_symbol(atom)}{marker}"]
    special_h = special_h_count(atom)
    if special_h is not None:
        fields.append("H" if special_h == 1 else f"H{special_h}")
    charge = int(atom.GetFormalCharge())
    if charge != 0:
        fields.append(charge_field(charge))
    degree_value = int(atom.GetDegree()) if degree is None else int(degree)
    suffix = "+" if allow_attachments else ""
    fields.append(f"D{degree_value}{suffix}")
    return f"[{';'.join(fields)}:{local_id}]"


def parse_aries_smarts(smarts: str) -> Chem.Mol | None:
    """Parse ARIES SMARTS after lowering custom attachment-tolerant degrees.

    ``Dn+`` means that ``n`` bonds are represented in the rule and additional
    unchanged saturated attachments may be present. RDKit performs the broad
    atom/bond match; :func:`_match_has_valid_attachments` enforces the custom
    boundary-bond condition.
    """

    lowered = re.sub(r";D\d+\+", "", smarts)
    return Chem.MolFromSmarts(lowered)


def fragment_smarts(
    mol: Chem.Mol,
    keys: set[NodeKey],
    side_key_by_idx: dict[int, NodeKey],
    local_ids: dict[NodeKey, int],
    retained_keys: set[NodeKey],
    attachment_keys: set[NodeKey],
    suppress_tetra_keys: set[NodeKey],
) -> str:
    indices = sorted(
        idx for idx, key in side_key_by_idx.items() if key in keys and key in local_ids
    )
    if not indices:
        return ""

    selected = set(indices)
    working = Chem.Mol(mol)
    for bond in working.GetBonds():
        # E/Z is serialized separately as structured metadata.  Clearing the
        # direction before SMILES generation preserves the underlying bond
        # type; replacing "/" or "\" afterwards can turn an aromatic bond
        # into an explicit single bond.
        bond.SetBondDir(Chem.BondDir.NONE)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    for atom in working.GetAtoms():
        atom.SetAtomMapNum(0)
    temp_to_idx: dict[int, int] = {}
    for atom_idx in indices:
        key = side_key_by_idx[atom_idx]
        temp_map = local_ids[key]
        working.GetAtomWithIdx(atom_idx).SetAtomMapNum(temp_map)
        temp_to_idx[temp_map] = atom_idx

    # RDKit's canonical traversal can start from a different atom when an atom
    # or bond state changes.  Root both sides at the same stable rule-local ID
    # so one graph rewrite always has one serialized representation.
    atom_mappings: list[tuple[int, ...]] = []
    source_fragments = Chem.GetMolFrags(
        working,
        asMols=True,
        sanitizeFrags=False,
        fragsMolAtomMapping=atom_mappings,
    )
    serialized_fragments = []
    selected_indices = set(indices)
    for source_fragment, atom_mapping in zip(
        source_fragments,
        atom_mappings,
    ):
        original_indices = [
            idx for idx in atom_mapping if idx in selected_indices
        ]
        if not original_indices:
            continue
        retained_indices = [
            idx
            for idx in original_indices
            if side_key_by_idx[idx] in retained_keys
        ]
        root_original_idx = min(
            retained_indices or original_indices,
            key=lambda idx: local_ids[side_key_by_idx[idx]],
        )
        original_to_fragment = {
            original_idx: fragment_idx
            for fragment_idx, original_idx in enumerate(atom_mapping)
        }
        serialized_fragments.append(
            Chem.MolFragmentToSmiles(
                source_fragment,
                atomsToUse=[
                    original_to_fragment[idx]
                    for idx in original_indices
                ],
                rootedAtAtom=original_to_fragment[root_original_idx],
                allBondsExplicit=True,
                canonical=True,
                isomericSmiles=True,
            )
        )
    fragment = ".".join(serialized_fragments)

    def replace(match: re.Match[str]) -> str:
        token = match.group(0)
        atom_idx = temp_to_idx[int(match.group(1))]
        key = side_key_by_idx[atom_idx]
        local_degree = sum(
            side_key_by_idx[int(neighbor.GetIdx())] in keys
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
        )
        return atom_query_token(
            mol.GetAtomWithIdx(atom_idx),
            local_ids[key],
            "" if key in suppress_tetra_keys else stereo_marker(token),
            degree=local_degree if key in attachment_keys else None,
            allow_attachments=key in attachment_keys,
        )

    serialized = re.sub(r"\[[^\]]*:(\d+)\]", replace, fragment)
    components = serialized.split(".")
    if len(components) > 1:
        def component_key(component: str) -> tuple[int, str]:
            local_maps = [
                int(value)
                for value in re.findall(r":(\d+)\]", component)
            ]
            return (min(local_maps, default=10**9), component)

        serialized = ".".join(sorted(components, key=component_key))
    return serialized


def _flip_local_tetra_marker(smarts: str, local_id: int) -> str:
    pattern = re.compile(rf"\[[^\]]*:{int(local_id)}\]")

    def replace(match: re.Match[str]) -> str:
        token = match.group(0)
        if "@@" in token:
            return token.replace("@@", "@", 1)
        if "@" in token:
            return token.replace("@", "@@", 1)
        return token

    return pattern.sub(replace, smarts, count=1)


def correct_fragment_tetra_smarts(
    smarts: str,
    mol: Chem.Mol,
    side_key_by_idx: dict[int, NodeKey],
    local_ids: dict[NodeKey, int],
) -> str:
    """Make the serialized fragment match the exact source stereochemistry.

    RDKit occasionally emits the opposite tetrahedral marker for a partial
    fragment, especially at a ring boundary.  RDChiral handles this with a
    randomized flip loop.  Here the intended atom correspondence is known, so
    each inconsistent local marker can be corrected directly and deterministically.
    """

    local_to_target = {
        local_ids[key]: atom_idx
        for atom_idx, key in side_key_by_idx.items()
        if key in local_ids
    }
    output = smarts
    for _ in range(3):
        query = parse_aries_smarts(output)
        if query is None:
            return output
        query_to_target = tuple(
            local_to_target[int(atom.GetAtomMapNum())] for atom in query.GetAtoms()
        )
        mismatched = []
        for query_atom in query.GetAtoms():
            if query_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            target_atom = mol.GetAtomWithIdx(
                query_to_target[int(query_atom.GetIdx())]
            )
            if not _chirality_matches_query(
                query_atom,
                target_atom,
                query_to_target,
            ):
                mismatched.append(int(query_atom.GetAtomMapNum()))
        if not mismatched:
            return output
        for local_id in mismatched:
            output = _flip_local_tetra_marker(output, local_id)
    return output


def strip_directional_bond_markers(smarts: str) -> str:
    """Keep E/Z out of inline SMARTS; it is explicit rule metadata instead."""

    return smarts.replace("/", "-").replace("\\", "-")


def double_bond_stereo_specs(
    mol: Chem.Mol,
    selected: set[NodeKey],
    key_by_idx: dict[int, NodeKey],
    local_ids: dict[NodeKey, int],
    allowed_bond_keys: set[frozenset[NodeKey]] | None = None,
) -> tuple[DoubleBondStereoSpec, ...]:
    specs = []
    for bond in mol.GetBonds():
        stereo = normalized_bond_stereo(bond)
        if stereo not in {"E", "Z"}:
            continue
        begin_key = key_by_idx[int(bond.GetBeginAtomIdx())]
        end_key = key_by_idx[int(bond.GetEndAtomIdx())]
        pair = frozenset((begin_key, end_key))
        if (
            allowed_bond_keys is not None
            and pair not in allowed_bond_keys
        ):
            continue
        if begin_key not in selected or end_key not in selected:
            continue
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) != 2:
            continue
        begin_neighbor_key = key_by_idx.get(int(stereo_atoms[0]))
        end_neighbor_key = key_by_idx.get(int(stereo_atoms[1]))
        if (
            begin_neighbor_key not in selected
            or end_neighbor_key not in selected
        ):
            continue
        specs.append(
            DoubleBondStereoSpec(
                begin_id=local_ids[begin_key],
                end_id=local_ids[end_key],
                begin_neighbor_id=local_ids[begin_neighbor_key],
                end_neighbor_id=local_ids[end_neighbor_key],
                stereo=stereo,
            )
        )
    return tuple(
        sorted(
            specs,
            key=lambda value: (
                value.begin_id,
                value.end_id,
                value.begin_neighbor_id,
                value.end_neighbor_id,
            ),
        )
    )


def side_keys(
    mol: Chem.Mol,
    side: Literal["r", "p"],
    retained_maps: set[int],
    reconciled_r_to_p: dict[int, int] | None = None,
) -> tuple[dict[int, NodeKey], dict[NodeKey, int]]:
    reconciled_r_to_p = reconciled_r_to_p or {}
    reconciled_p_to_r = {
        product_map: reactant_map
        for reactant_map, product_map in reconciled_r_to_p.items()
    }
    by_idx: dict[int, NodeKey] = {}
    by_key: dict[NodeKey, int] = {}
    for atom in mol.GetAtoms():
        atom_idx = int(atom.GetIdx())
        atom_map = int(atom.GetAtomMapNum())
        if atom_map > 0 and atom_map in retained_maps:
            key = ("m", atom_map)
        elif side == "r" and atom_map in reconciled_r_to_p:
            key = ("x", atom_map)
        elif side == "p" and atom_map in reconciled_p_to_r:
            key = ("x", reconciled_p_to_r[atom_map])
        elif atom_map > 0:
            key = (f"{side}m", atom_map)
        else:
            key = (f"{side}i", atom_idx)
        by_idx[atom_idx] = key
        by_key[key] = atom_idx
    return by_idx, by_key


def reconcile_unmatched_atom_maps(
    reactant: Chem.Mol,
    product: Chem.Mol,
    r_maps: dict[int, Chem.Atom],
    p_maps: dict[int, Chem.Atom],
) -> dict[int, int]:
    """Pair conservatively remapped atoms using already paired neighbors.

    Existing equal map numbers are authoritative.  A repair requires the same
    element and at least one mapped neighboring atom whose correspondence is
    already known.  This handles local map-number replacement at reaction sites
    without attempting a fresh global atom mapping.
    """

    common = set(r_maps) & set(p_maps)
    correspondence = {atom_map: atom_map for atom_map in common}
    repaired: dict[int, int] = {}
    remaining_r = set(r_maps) - common
    remaining_p = set(p_maps) - common

    while remaining_r and remaining_p:
        candidates: list[
            tuple[
                tuple[int, int, int, int, int],
                int,
                int,
            ]
        ] = []
        for reactant_map in remaining_r:
            r_atom = r_maps[reactant_map]
            known_r_neighbors = {
                int(neighbor.GetAtomMapNum()): bond_type_name(
                    reactant.GetBondBetweenAtoms(
                        int(r_atom.GetIdx()),
                        int(neighbor.GetIdx()),
                    )
                )
                for neighbor in r_atom.GetNeighbors()
                if int(neighbor.GetAtomMapNum()) in correspondence
            }
            if not known_r_neighbors:
                continue
            for product_map in remaining_p:
                p_atom = p_maps[product_map]
                if int(r_atom.GetAtomicNum()) != int(p_atom.GetAtomicNum()):
                    continue
                p_neighbor_bonds = {
                    int(neighbor.GetAtomMapNum()): bond_type_name(
                        product.GetBondBetweenAtoms(
                            int(p_atom.GetIdx()),
                            int(neighbor.GetIdx()),
                        )
                    )
                    for neighbor in p_atom.GetNeighbors()
                }
                matched_neighbors = 0
                matched_bond_types = 0
                for r_neighbor_map, r_bond_type in known_r_neighbors.items():
                    expected_p_map = correspondence[r_neighbor_map]
                    if expected_p_map not in p_neighbor_bonds:
                        continue
                    matched_neighbors += 1
                    matched_bond_types += int(
                        p_neighbor_bonds[expected_p_map] == r_bond_type
                    )
                if matched_neighbors == 0:
                    continue
                score = (
                    matched_neighbors,
                    matched_bond_types,
                    int(r_atom.GetIsAromatic() == p_atom.GetIsAromatic()),
                    -abs(int(r_atom.GetDegree()) - int(p_atom.GetDegree())),
                    -abs(
                        int(r_atom.GetFormalCharge())
                        - int(p_atom.GetFormalCharge())
                    ),
                )
                candidates.append((score, reactant_map, product_map))

        if not candidates:
            break
        candidates.sort(reverse=True)
        used_r: set[int] = set()
        used_p: set[int] = set()
        accepted: list[tuple[int, int]] = []
        for score, reactant_map, product_map in candidates:
            if reactant_map in used_r or product_map in used_p:
                continue
            # A single matched neighbor is accepted only when the local degree
            # and aromatic class agree; otherwise it may be a true replacement.
            if score[0] == 1 and (score[2] == 0 or score[3] < -1):
                continue
            used_r.add(reactant_map)
            used_p.add(product_map)
            accepted.append((reactant_map, product_map))
        if not accepted:
            break
        for reactant_map, product_map in accepted:
            correspondence[reactant_map] = product_map
            repaired[reactant_map] = product_map
            remaining_r.remove(reactant_map)
            remaining_p.remove(product_map)
    return repaired


def bond_lookup(
    mol: Chem.Mol,
    key_by_idx: dict[int, NodeKey],
) -> dict[frozenset[NodeKey], Chem.Bond]:
    result = {}
    for bond in mol.GetBonds():
        left = key_by_idx[int(bond.GetBeginAtomIdx())]
        right = key_by_idx[int(bond.GetEndAtomIdx())]
        result[frozenset((left, right))] = bond
    return result


def collapsible_reactant_payload_keys(
    r_idx_by_key: dict[NodeKey, int],
    p_idx_by_key: dict[NodeKey, int],
    r_bonds: dict[frozenset[NodeKey], Chem.Bond],
) -> set[NodeKey]:
    """Find internal atoms of deleted branches with one retained attachment.

    The atom directly bonded to the retained graph remains in the reacting
    center and defines the leaving-group chemistry.  Its downstream branch is
    application payload, not reusable template context.
    """

    reactant_only = set(r_idx_by_key) - set(p_idx_by_key)
    adjacency: dict[NodeKey, set[NodeKey]] = {
        key: set() for key in reactant_only
    }
    for pair in r_bonds:
        values = tuple(pair)
        if len(values) != 2:
            continue
        left, right = values
        if left in reactant_only and right in reactant_only:
            adjacency[left].add(right)
            adjacency[right].add(left)

    payload: set[NodeKey] = set()
    unseen = set(reactant_only)
    while unseen:
        start = unseen.pop()
        component = {start}
        frontier = [start]
        while frontier:
            current = frontier.pop()
            for neighbor in adjacency[current]:
                if neighbor not in unseen:
                    continue
                unseen.remove(neighbor)
                component.add(neighbor)
                frontier.append(neighbor)

        boundary_edges: list[tuple[NodeKey, NodeKey]] = []
        for pair in r_bonds:
            if not pair.intersection(component):
                continue
            values = tuple(pair)
            if len(values) != 2:
                continue
            left, right = values
            if left in component and right not in component:
                boundary_edges.append((left, right))
            elif right in component and left not in component:
                boundary_edges.append((right, left))
        if len(boundary_edges) != 1:
            continue
        root, retained_neighbor = boundary_edges[0]
        if retained_neighbor not in p_idx_by_key:
            continue
        payload.update(component - {root})
    return payload


def bonds_differ(left: Chem.Bond | None, right: Chem.Bond | None) -> bool:
    if left is None or right is None:
        return left is not right
    return (
        bond_type_name(left) != bond_type_name(right)
        or bool(left.GetIsAromatic()) != bool(right.GetIsAromatic())
        or normalized_bond_stereo(left) != normalized_bond_stereo(right)
    )


def is_hetero_context(atom: Chem.Atom) -> bool:
    """Use retained one-hop heteroatoms as explicit context."""

    return int(atom.GetAtomicNum()) not in {1, 6}


def _accfg_lite():
    """Construct the default AccFG-lite detector without ring groups."""

    global _ACCFG_LITE
    if _ACCFG_LITE is None:
        from accfg import AccFG

        _ACCFG_LITE = AccFG(
            lite=True,
            exclude_fgs=["rings"],
        )
    return _ACCFG_LITE


def accfg_atom_groups(mol: Chem.Mol) -> tuple[frozenset[int], ...]:
    """Return deterministic distinct AccFG-lite atom-index groups."""

    matches = _accfg_lite().run_mol(
        mol,
        show_atoms=True,
        show_graph=False,
    )
    groups = {
        frozenset(int(atom_idx) for atom_idx in atom_indices)
        for name in sorted(matches)
        for atom_indices in matches[name]
        if atom_indices
    }
    return tuple(
        sorted(
            groups,
            key=lambda group: (len(group), tuple(sorted(group))),
        )
    )


def _add_functional_group_atoms(
    selected: set[NodeKey],
    keys: Iterable[NodeKey],
    mol: Chem.Mol,
    idx_by_key: dict[NodeKey, int],
    key_by_idx: dict[int, NodeKey],
) -> set[NodeKey]:
    """Add every AccFG-lite group containing at least one selected center."""

    center_indices = {
        int(idx_by_key[key])
        for key in keys
        if key in idx_by_key
    }
    if not center_indices:
        return set()
    added: set[NodeKey] = set()
    for group in accfg_atom_groups(mol):
        if group.isdisjoint(center_indices):
            continue
        for atom_idx in group:
            key = key_by_idx[int(atom_idx)]
            if key not in selected:
                added.add(key)
    selected.update(added)
    return added


def _external_bonds_are_saturated(
    mol: Chem.Mol,
    key: NodeKey,
    selected: set[NodeKey],
    idx_by_key: dict[NodeKey, int],
    key_by_idx: dict[int, NodeKey],
) -> bool:
    atom_idx = idx_by_key.get(key)
    if atom_idx is None:
        return False
    atom = mol.GetAtomWithIdx(atom_idx)
    for bond in atom.GetBonds():
        neighbor_idx = int(bond.GetOtherAtomIdx(atom_idx))
        if key_by_idx[neighbor_idx] in selected:
            continue
        if bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
            return False
    return True


def attachment_tolerant_keys(
    reactant: Chem.Mol,
    product: Chem.Mol,
    selected: set[NodeKey],
    center: set[NodeKey],
    r_idx_by_key: dict[NodeKey, int],
    p_idx_by_key: dict[NodeKey, int],
    r_key_by_idx: dict[int, NodeKey],
    p_key_by_idx: dict[int, NodeKey],
) -> set[NodeKey]:
    """Find retained center atoms that may carry extra saturated substituents."""

    result: set[NodeKey] = set()
    for key in center:
        if key not in r_idx_by_key or key not in p_idx_by_key:
            continue
        r_atom = reactant.GetAtomWithIdx(r_idx_by_key[key])
        p_atom = product.GetAtomWithIdx(p_idx_by_key[key])
        if int(r_atom.GetAtomicNum()) <= 1:
            continue
        if int(r_atom.GetAtomicNum()) != int(p_atom.GetAtomicNum()):
            continue
        if r_atom.GetIsAromatic() or p_atom.GetIsAromatic():
            continue
        if not _external_bonds_are_saturated(
            reactant,
            key,
            selected,
            r_idx_by_key,
            r_key_by_idx,
        ):
            continue
        if not _external_bonds_are_saturated(
            product,
            key,
            selected,
            p_idx_by_key,
            p_key_by_idx,
        ):
            continue
        result.add(key)
    return result


def _add_neighbors(
    selected: set[NodeKey],
    keys: Iterable[NodeKey],
    mol: Chem.Mol,
    idx_by_key: dict[NodeKey, int],
    key_by_idx: dict[int, NodeKey],
    predicate,
) -> set[NodeKey]:
    added: set[NodeKey] = set()
    for key in list(keys):
        atom_idx = idx_by_key.get(key)
        if atom_idx is None:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if predicate(neighbor):
                neighbor_key = key_by_idx[int(neighbor.GetIdx())]
                if neighbor_key not in selected:
                    added.add(neighbor_key)
    selected.update(added)
    return added


def _stereo_support(
    reactant: Chem.Mol,
    product: Chem.Mol,
    selected: set[NodeKey],
    center: set[NodeKey],
    r_idx_by_key: dict[NodeKey, int],
    p_idx_by_key: dict[NodeKey, int],
    r_key_by_idx: dict[int, NodeKey],
    p_key_by_idx: dict[int, NodeKey],
) -> set[NodeKey]:
    """Add atoms needed to represent stereochemistry at a reacting center.

    Context atoms describe only the immediate retained neighborhood.  They do
    not seed this expansion, because doing so can recursively pull an entire
    unchanged E/Z motif into an otherwise local context rule.
    """

    support: set[NodeKey] = set()
    for mol, idx_by_key, key_by_idx in (
        (reactant, r_idx_by_key, r_key_by_idx),
        (product, p_idx_by_key, p_key_by_idx),
    ):
        stereo_centers: set[NodeKey] = set()
        for key in list(center):
            atom_idx = idx_by_key.get(key)
            if atom_idx is None:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                if key in r_idx_by_key and key in p_idx_by_key:
                    relation = tetra_relation(
                        reactant.GetAtomWithIdx(r_idx_by_key[key]),
                        product.GetAtomWithIdx(p_idx_by_key[key]),
                    )
                else:
                    relation = "side-only"
                if relation not in {"erase", "retain"}:
                    stereo_centers.add(key)
            for bond in atom.GetBonds():
                current_stereo = normalized_bond_stereo(bond)
                if current_stereo not in {"E", "Z"}:
                    continue
                other_key = key_by_idx[
                    int(bond.GetOtherAtomIdx(atom_idx))
                ]
                counterpart = None
                if mol is reactant:
                    if key in p_idx_by_key and other_key in p_idx_by_key:
                        counterpart = product.GetBondBetweenAtoms(
                            p_idx_by_key[key],
                            p_idx_by_key[other_key],
                        )
                elif key in r_idx_by_key and other_key in r_idx_by_key:
                    counterpart = reactant.GetBondBetweenAtoms(
                        r_idx_by_key[key],
                        r_idx_by_key[other_key],
                    )
                counterpart_stereo = (
                    normalized_bond_stereo(counterpart)
                    if counterpart is not None
                    else ""
                )
                if current_stereo == counterpart_stereo:
                    continue
                # Erased reactant E/Z must not constrain application.
                if (
                    mol is reactant
                    and counterpart_stereo not in {"E", "Z"}
                ):
                    continue
                stereo_centers.add(key)
                stereo_centers.add(other_key)

        for key in stereo_centers:
            atom_idx = idx_by_key.get(key)
            if atom_idx is None:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if key not in selected:
                support.add(key)
            for neighbor in atom.GetNeighbors():
                support.add(key_by_idx[int(neighbor.GetIdx())])

    support.difference_update(selected)
    selected.update(support)
    return support


def _node_sort_key(
    key: NodeKey,
    graph: nx.Graph,
) -> tuple[str, tuple[str, ...], str]:
    label = str(graph.nodes[key]["label"])
    neighbors = tuple(
        sorted(
            f"{graph.edges[key, neighbor]['label']}|{graph.nodes[neighbor]['label']}"
            for neighbor in graph.neighbors(key)
        )
    )
    return label, neighbors, repr(key)


def _graph_hash(graph: nx.Graph) -> str:
    return nx.weisfeiler_lehman_graph_hash(
        graph,
        node_attr="label",
        edge_attr="label",
        iterations=4,
    )


def _graphs_match(left: nx.Graph, right: nx.Graph) -> bool:
    return nx.is_isomorphic(
        left,
        right,
        node_match=lambda a, b: a["label"] == b["label"],
        edge_match=lambda a, b: a["label"] == b["label"],
    )


def _source_assignment_in_reference(
    reference: _CandidateRule,
    candidate: _CandidateRule,
) -> tuple[tuple[tuple[int, int], ...], ...]:
    """Express one candidate's observed embedding in the reference local IDs."""

    if reference is candidate:
        return (candidate.source_assignment,)
    candidate_local_to_map = dict(candidate.source_assignment)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        reference.graph,
        candidate.graph,
        node_match=lambda a, b: a["label"] == b["label"],
        edge_match=lambda a, b: a["label"] == b["label"],
    )
    assignments = []
    seen = set()
    for mapping_index, mapping in enumerate(matcher.isomorphisms_iter()):
        if mapping_index >= 128:
            break
        assignment = []
        for reference_key, candidate_key in mapping.items():
            candidate_local = candidate.local_ids[candidate_key]
            atom_map = candidate_local_to_map.get(candidate_local)
            if atom_map is None:
                continue
            assignment.append((reference.local_ids[reference_key], atom_map))
        normalized = tuple(sorted(assignment))
        if normalized not in seen:
            seen.add(normalized)
            assignments.append(normalized)
    return tuple(assignments)


def _candidate_graph(
    reactant: Chem.Mol,
    product: Chem.Mol,
    selected: set[NodeKey],
    center: set[NodeKey],
    attachment: set[NodeKey],
    context: set[NodeKey],
    stereo: set[NodeKey],
    r_idx_by_key: dict[NodeKey, int],
    p_idx_by_key: dict[NodeKey, int],
    r_bonds: dict[frozenset[NodeKey], Chem.Bond],
    p_bonds: dict[frozenset[NodeKey], Chem.Bond],
) -> nx.Graph:
    def graph_atom_state(
        atom: Chem.Atom,
        key: NodeKey,
        bonds: dict[frozenset[NodeKey], Chem.Bond],
    ) -> AtomState:
        state = AtomState.from_atom(atom)
        if key not in attachment:
            return state
        local_degree = sum(
            key in pair and pair.issubset(selected)
            for pair in bonds
        )
        return replace(state, degree=int(local_degree))

    graph = nx.Graph()
    for key in selected:
        r_atom = (
            reactant.GetAtomWithIdx(r_idx_by_key[key]) if key in r_idx_by_key else None
        )
        p_atom = product.GetAtomWithIdx(p_idx_by_key[key]) if key in p_idx_by_key else None
        if key in center:
            role = "center"
        elif key in context:
            role = "context"
        else:
            role = "stereo"
        relation = (
            tetra_relation(r_atom, p_atom)
            if r_atom is not None and p_atom is not None
            else "side-only"
        )
        r_cip = r_atom.GetProp("_CIPCode") if r_atom and r_atom.HasProp("_CIPCode") else ""
        p_cip = p_atom.GetProp("_CIPCode") if p_atom and p_atom.HasProp("_CIPCode") else ""
        if relation == "retain":
            relation = "none"
            r_cip = ""
            p_cip = ""
        label = json.dumps(
            {
                "role": role,
                "r": (
                    asdict(graph_atom_state(r_atom, key, r_bonds))
                    if r_atom
                    else None
                ),
                "p": (
                    asdict(graph_atom_state(p_atom, key, p_bonds))
                    if p_atom
                    else None
                ),
                "tetra": relation,
                "r_cip": r_cip,
                "p_cip": p_cip,
            },
            sort_keys=True,
        )
        graph.add_node(key, label=label)

    all_pairs = set(r_bonds) | set(p_bonds)
    for pair in all_pairs:
        if not pair.issubset(selected):
            continue
        values = tuple(pair)
        if len(values) != 2:
            continue
        r_bond = r_bonds.get(pair)
        p_bond = p_bonds.get(pair)
        r_state = BondState.from_bond(r_bond) if r_bond else None
        p_state = BondState.from_bond(p_bond) if p_bond else None
        if (
            r_state is not None
            and p_state is not None
            and r_state.stereo == p_state.stereo
        ):
            r_state = replace(r_state, stereo="")
            p_state = replace(p_state, stereo="")
        label = json.dumps(
            {
                "r": asdict(r_state) if r_state else None,
                "p": asdict(p_state) if p_state else None,
            },
            sort_keys=True,
        )
        graph.add_edge(values[0], values[1], label=label)
    return graph


def _extract_candidates(
    reactant: Chem.Mol,
    product: Chem.Mol,
    mode: Mode,
) -> tuple[list[_CandidateRule], tuple[SourceTetraSpec, ...]]:
    r_maps = map_to_atom(reactant)
    p_maps = map_to_atom(product)
    retained_maps = set(r_maps) & set(p_maps)
    reconciled_r_to_p = reconcile_unmatched_atom_maps(
        reactant,
        product,
        r_maps,
        p_maps,
    )
    r_key_by_idx, r_idx_by_key = side_keys(
        reactant,
        "r",
        retained_maps,
        reconciled_r_to_p,
    )
    p_key_by_idx, p_idx_by_key = side_keys(
        product,
        "p",
        retained_maps,
        reconciled_r_to_p,
    )
    r_bonds = bond_lookup(reactant, r_key_by_idx)
    p_bonds = bond_lookup(product, p_key_by_idx)
    reactant_payload_keys = collapsible_reactant_payload_keys(
        r_idx_by_key,
        p_idx_by_key,
        r_bonds,
    )
    reactant_active_double_stereo: set[frozenset[NodeKey]] = set()
    product_active_double_stereo: set[frozenset[NodeKey]] = set()
    erased_double_stereo: set[frozenset[NodeKey]] = set()
    for pair in set(r_bonds) | set(p_bonds):
        r_stereo = (
            normalized_bond_stereo(r_bonds[pair])
            if pair in r_bonds
            else ""
        )
        p_stereo = (
            normalized_bond_stereo(p_bonds[pair])
            if pair in p_bonds
            else ""
        )
        if p_stereo in {"E", "Z"} and p_stereo != r_stereo:
            product_active_double_stereo.add(pair)
        if r_stereo in {"E", "Z"} and p_stereo not in {"E", "Z"}:
            erased_double_stereo.add(pair)
        if (
            r_stereo in {"E", "Z"}
            and p_stereo in {"E", "Z"}
            and r_stereo != p_stereo
        ):
            reactant_active_double_stereo.add(pair)

    changes = _DisjointSet()
    tetra_relations: dict[NodeKey, str] = {}
    mapped_tetra_relations: dict[NodeKey, str] = {}

    for atom_map in retained_maps:
        key = ("m", atom_map)
        left = r_maps[atom_map]
        right = p_maps[atom_map]
        relation = tetra_relation(left, right)
        tetra_relations[key] = relation
        mapped_tetra_relations[key] = mapped_tetra_relation(left, right)
        if atom_core_state(left) != atom_core_state(right):
            changes.add(key)
    for reactant_map, product_map in reconciled_r_to_p.items():
        key = ("x", reactant_map)
        left = r_maps[reactant_map]
        right = p_maps[product_map]
        relation = tetra_relation(left, right)
        tetra_relations[key] = relation
        mapped_tetra_relations[key] = mapped_tetra_relation(left, right)
        if atom_core_state(left) != atom_core_state(right):
            changes.add(key)

    all_bond_pairs = set(r_bonds) | set(p_bonds)
    for pair in all_bond_pairs:
        if pair.intersection(reactant_payload_keys):
            continue
        left = r_bonds.get(pair)
        right = p_bonds.get(pair)
        if not bonds_differ(left, right):
            continue
        values = tuple(pair)
        if len(values) == 1:
            changes.add(values[0])
        else:
            changes.union(values[0], values[1])

    chemical_changed_keys = set(changes.parent)
    stereo_change_keys = {
        key
        for key, relation in tetra_relations.items()
        if relation in {"create", "erase", "invert"}
    }
    promoted_stereo_keys: set[NodeKey] = set()
    if not chemical_changed_keys:
        for key in stereo_change_keys:
            changes.add(key)
            promoted_stereo_keys.add(key)
    else:
        for key in stereo_change_keys:
            if key in chemical_changed_keys:
                promoted_stereo_keys.add(key)
                continue
            adjacent_changed = set()
            for mol, idx_by_key, key_by_idx in (
                (reactant, r_idx_by_key, r_key_by_idx),
                (product, p_idx_by_key, p_key_by_idx),
            ):
                atom_idx = idx_by_key.get(key)
                if atom_idx is None:
                    continue
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    neighbor_key = key_by_idx[int(neighbor.GetIdx())]
                    if neighbor_key in chemical_changed_keys:
                        adjacent_changed.add(neighbor_key)
            for neighbor_key in adjacent_changed:
                changes.union(key, neighbor_key)
            if adjacent_changed:
                promoted_stereo_keys.add(key)

    if mode == "fg":
        # One functional group owns one rewrite. If an AccFG match spans
        # otherwise independent changed components, merge those components
        # before adding the remaining group atoms as context. This prevents
        # duplicate creation/deletion of the same side-specific FG atoms.
        for mol, key_by_idx in (
            (reactant, r_key_by_idx),
            (product, p_key_by_idx),
        ):
            for group in accfg_atom_groups(mol):
                changed_group_keys = sorted(
                    (
                        key_by_idx[int(atom_idx)]
                        for atom_idx in group
                        if key_by_idx[int(atom_idx)] in changes.parent
                    ),
                    key=str,
                )
                if len(changed_group_keys) < 2:
                    continue
                anchor = changed_group_keys[0]
                for key in changed_group_keys[1:]:
                    changes.union(anchor, key)

    source_only_tetra_keys = {
        key
        for key, relation in mapped_tetra_relations.items()
        if relation == "invert"
        and tetra_relations.get(key) == "retain"
    }
    remote_tetra_specs = []
    for key in sorted(
        (stereo_change_keys - promoted_stereo_keys)
        | source_only_tetra_keys
    ):
        if key[0] not in {"m", "x"} or int(key[1]) <= 0:
            continue
        product_idx = p_idx_by_key.get(key)
        if product_idx is None:
            continue
        product_atom = product.GetAtomWithIdx(product_idx)
        neighbor_maps = []
        for neighbor in product_atom.GetNeighbors():
            neighbor_key = p_key_by_idx[int(neighbor.GetIdx())]
            neighbor_maps.append(
                int(neighbor_key[1])
                if neighbor_key[0] in {"m", "x"}
                else 0
            )
        remote_tetra_specs.append(
            SourceTetraSpec(
                atom_map=int(key[1]),
                action=(
                    mapped_tetra_relations[key]
                    if key in source_only_tetra_keys
                    else tetra_relations[key]
                ),
                product_tag=str(product_atom.GetChiralTag()),
                product_neighbor_maps=tuple(neighbor_maps),
            )
        )

    components = sorted(
        changes.components(),
        key=lambda component: tuple(sorted(component)),
    )
    if not components:
        raise ValueError("no heavy-atom, state, or stereochemical change found")

    candidates: list[_CandidateRule] = []
    for center in components:
        selected = set(center)
        context: set[NodeKey] = set()
        if mode == "context":
            context.update(
                _add_neighbors(
                    selected,
                    center,
                    reactant,
                    r_idx_by_key,
                    r_key_by_idx,
                    is_hetero_context,
                )
            )
            # Context describes retained surroundings. A side-only atom is an
            # edit owned by its own change component, not context for this one.
            side_only_context = {
                key
                for key in context
                if key not in r_idx_by_key or key not in p_idx_by_key
            }
            context.difference_update(side_only_context)
            selected.difference_update(side_only_context)
            context.update(
                _add_neighbors(
                    selected,
                    center,
                    product,
                    p_idx_by_key,
                    p_key_by_idx,
                    is_hetero_context,
                )
            )
        elif mode == "fg":
            context.update(
                _add_functional_group_atoms(
                    selected,
                    center,
                    reactant,
                    r_idx_by_key,
                    r_key_by_idx,
                )
            )
            context.update(
                _add_functional_group_atoms(
                    selected,
                    center,
                    product,
                    p_idx_by_key,
                    p_key_by_idx,
                )
            )
            # Functional groups are side-specific. A complete group may
            # therefore contain atoms present only in the reactant or product.
            # Keep those atoms explicit without adding them to the minimal
            # reacting-center atom set.

        stereo = _stereo_support(
            reactant,
            product,
            selected,
            center,
            r_idx_by_key,
            p_idx_by_key,
            r_key_by_idx,
            p_key_by_idx,
        )
        # As with context, side-only stereo neighbors belong to the component
        # that creates or deletes them.  Retained support atoms are sufficient
        # for matching and parity alignment.
        side_only_stereo = {
            key
            for key in stereo
            if key not in r_idx_by_key or key not in p_idx_by_key
        }
        stereo.difference_update(side_only_stereo)
        selected.difference_update(side_only_stereo)

        attachment = attachment_tolerant_keys(
            reactant,
            product,
            selected,
            center,
            r_idx_by_key,
            p_idx_by_key,
            r_key_by_idx,
            p_key_by_idx,
        )
        suppress_tetra = {
            key
            for key in selected
            if key in r_idx_by_key
            and key in p_idx_by_key
            and tetra_relation(
                reactant.GetAtomWithIdx(r_idx_by_key[key]),
                product.GetAtomWithIdx(p_idx_by_key[key]),
            )
            in {"erase", "retain"}
        }
        graph = _candidate_graph(
            reactant,
            product,
            selected,
            center,
            attachment,
            context,
            stereo,
            r_idx_by_key,
            p_idx_by_key,
            r_bonds,
            p_bonds,
        )
        ordered_keys = sorted(selected, key=lambda key: _node_sort_key(key, graph))
        local_ids = {key: index for index, key in enumerate(ordered_keys, start=1)}
        retained_keys = set(r_idx_by_key).intersection(p_idx_by_key, selected)
        reactant_smarts = fragment_smarts(
            reactant,
            selected,
            r_key_by_idx,
            local_ids,
            retained_keys,
            attachment,
            suppress_tetra,
        )
        reactant_smarts = correct_fragment_tetra_smarts(
            reactant_smarts,
            reactant,
            r_key_by_idx,
            local_ids,
        )
        reactant_double_stereo = double_bond_stereo_specs(
            reactant,
            selected,
            r_key_by_idx,
            local_ids,
            reactant_active_double_stereo,
        )
        reactant_smarts = strip_directional_bond_markers(reactant_smarts)
        product_smarts = fragment_smarts(
            product,
            selected,
            p_key_by_idx,
            local_ids,
            retained_keys,
            attachment,
            suppress_tetra,
        )
        product_double_stereo = double_bond_stereo_specs(
            product,
            selected,
            p_key_by_idx,
            local_ids,
            product_active_double_stereo,
        )
        product_double_stereo += tuple(
            DoubleBondStereoSpec(
                begin_id=local_ids[values[0]],
                end_id=local_ids[values[1]],
                begin_neighbor_id=0,
                end_neighbor_id=0,
                stereo="erase",
            )
            for pair in sorted(
                erased_double_stereo,
                key=lambda value: tuple(sorted(value)),
            )
            if pair.issubset(selected)
            and pair in p_bonds
            and p_bonds[pair].GetBondType() == Chem.BondType.DOUBLE
            for values in [tuple(sorted(pair))]
        )
        product_smarts = strip_directional_bond_markers(product_smarts)
        product_smarts = correct_fragment_tetra_smarts(
            product_smarts,
            product,
            p_key_by_idx,
            local_ids,
        )
        reactant_states = {
            local_ids[key]: AtomState.from_atom(
                reactant.GetAtomWithIdx(r_idx_by_key[key])
            )
            for key in selected
            if key in r_idx_by_key
        }
        product_states = {
            local_ids[key]: AtomState.from_atom(product.GetAtomWithIdx(p_idx_by_key[key]))
            for key in selected
            if key in p_idx_by_key
        }
        tetra_actions = {}
        for key in center:
            if key not in r_idx_by_key or key not in p_idx_by_key:
                continue
            relation = tetra_relation(
                reactant.GetAtomWithIdx(r_idx_by_key[key]),
                product.GetAtomWithIdx(p_idx_by_key[key]),
            )
            if relation in {"create", "erase", "invert"}:
                tetra_actions[local_ids[key]] = relation
        occurrence_maps = tuple(
            sorted(
                key[1]
                for key in center
                if key[0] in {"m", "x"} and key[1] > 0
            )
        )
        source_assignment = tuple(
            sorted(
                (
                    local_ids[key],
                    int(
                        reactant.GetAtomWithIdx(
                            r_idx_by_key[key]
                        ).GetAtomMapNum()
                    ),
                )
                for key in selected
                if key in r_idx_by_key
                and int(
                    reactant.GetAtomWithIdx(
                        r_idx_by_key[key]
                    ).GetAtomMapNum()
                )
                > 0
            )
        )
        candidates.append(
            _CandidateRule(
                graph=graph,
                selected_keys=selected,
                center_keys=set(center),
                attachment_keys=attachment,
                changed_bond_keys={
                    pair
                    for pair in set(r_bonds) | set(p_bonds)
                    if pair.issubset(center)
                    and bonds_differ(r_bonds.get(pair), p_bonds.get(pair))
                },
                reactant_double_stereo=reactant_double_stereo,
                product_double_stereo=product_double_stereo,
                tetra_actions=tetra_actions,
                context_keys=context,
                stereo_keys=stereo,
                occurrence_maps=occurrence_maps,
                source_assignment=source_assignment,
                local_ids=local_ids,
                reactant_smarts=reactant_smarts,
                product_smarts=product_smarts,
                reactant_states=reactant_states,
                product_states=product_states,
                graph_hash=_graph_hash(graph),
            )
        )
    return candidates, tuple(remote_tetra_specs)


def extract_aries_template(
    mapped_substrate_smiles: str,
    mapped_product_smiles: str,
    mode: Mode = "context",
) -> AriesTemplate:
    if mode not in {"core", "context", "fg"}:
        raise ValueError(f"unsupported ARIES RXN mode: {mode}")
    try:
        reactant = parse_mol(mapped_substrate_smiles, "mapped substrate")
        product = parse_mol(mapped_product_smiles, "mapped product")
        candidates, remote_tetra_specs = _extract_candidates(
            reactant,
            product,
            mode,
        )

        groups: list[list[_CandidateRule]] = []
        for candidate in candidates:
            placed = False
            for group in groups:
                reference = group[0]
                if (
                    candidate.graph_hash == reference.graph_hash
                    and _graphs_match(candidate.graph, reference.graph)
                ):
                    group.append(candidate)
                    placed = True
                    break
            if not placed:
                groups.append([candidate])

        rules = []
        for group in groups:
            reference = group[0]
            rules.append(
                AriesRule(
                    reactant_smarts=reference.reactant_smarts,
                    product_smarts=reference.product_smarts,
                    center_ids=tuple(
                        sorted(reference.local_ids[key] for key in reference.center_keys)
                    ),
                    attachment_ids=tuple(
                        sorted(
                            reference.local_ids[key]
                            for key in reference.attachment_keys
                        )
                    ),
                    changed_bond_ids=tuple(
                        sorted(
                            tuple(
                                sorted(reference.local_ids[key] for key in pair)
                            )
                            for pair in reference.changed_bond_keys
                            if len(pair) == 2
                        )
                    ),
                    reactant_double_stereo=reference.reactant_double_stereo,
                    product_double_stereo=reference.product_double_stereo,
                    tetra_actions=reference.tetra_actions,
                    context_ids=tuple(
                        sorted(reference.local_ids[key] for key in reference.context_keys)
                    ),
                    stereo_support_ids=tuple(
                        sorted(reference.local_ids[key] for key in reference.stereo_keys)
                    ),
                    observed_repeats=len(group),
                    source_occurrences=tuple(
                        candidate.occurrence_maps for candidate in group
                    ),
                    source_assignments=tuple(
                        assignment
                        for candidate in group
                        for assignment in _source_assignment_in_reference(
                            reference,
                            candidate,
                        )
                    ),
                    graph_hash=reference.graph_hash,
                    reactant_states=reference.reactant_states,
                    product_states=reference.product_states,
                )
            )
        rules.sort(key=lambda rule: (rule.graph_hash, rule.smarts))
        return AriesTemplate(
            mode=mode,
            rules=tuple(rules),
            source_mapped_substrate=canonical_connectivity_smiles_from_mol(
                reactant,
                keep_maps=True,
            ),
            source_remote_tetra_specs=remote_tetra_specs,
        )
    except Exception as exc:
        return AriesTemplate(
            mode=mode,
            rules=(),
            status="error",
            error=f"{type(exc).__name__}: {exc}",
        )


def query_atom_maps(query: Chem.Mol) -> dict[int, int]:
    result = {}
    for atom in query.GetAtoms():
        local_id = int(atom.GetAtomMapNum())
        if local_id > 0:
            result[local_id] = int(atom.GetIdx())
    return result


def _chirality_matches_query(
    query_atom: Chem.Atom,
    target_atom: Chem.Atom,
    query_to_target: Sequence[int],
) -> bool:
    if query_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
        return True
    if target_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
        return not target_atom.HasProp("_ChiralityPossible")

    target_to_query = {
        int(target_idx): int(query_idx)
        for query_idx, target_idx in enumerate(query_to_target)
    }
    query_labels = [
        int(neighbor.GetAtomMapNum())
        for neighbor in query_atom.GetNeighbors()
    ]
    target_labels = []
    for neighbor in target_atom.GetNeighbors():
        query_idx = target_to_query.get(int(neighbor.GetIdx()))
        if query_idx is None:
            target_labels.append(1000000 + int(neighbor.GetIdx()))
        else:
            target_labels.append(
                int(query_atom.GetOwningMol().GetAtomWithIdx(query_idx).GetAtomMapNum())
            )
    if len(query_labels) == 3:
        query_labels.append(-1)
    if len(target_labels) == 3:
        target_labels.append(-1)
    if len(query_labels) != 4 or len(target_labels) != 4:
        return True
    if set(query_labels) != set(target_labels):
        return True
    rank = {value: index for index, value in enumerate(sorted(query_labels))}
    query_parity = parity([rank[value] for value in query_labels])
    target_parity = parity([rank[value] for value in target_labels])
    tags_same = query_atom.GetChiralTag() == target_atom.GetChiralTag()
    return (query_parity == target_parity) == tags_same


def _bond_stereo_matches_query(
    query_bond: Chem.Bond,
    target_bond: Chem.Bond,
) -> bool:
    expected = normalized_bond_stereo(query_bond)
    if not expected:
        return True
    observed = normalized_bond_stereo(target_bond)
    if not observed:
        return False
    return expected == observed


def _match_has_required_stereo(
    query: Chem.Mol,
    substrate: Chem.Mol,
    match: Sequence[int],
    ignored_tetra_local_ids: Collection[int] = (),
) -> bool:
    ignored = {int(value) for value in ignored_tetra_local_ids}
    for query_atom in query.GetAtoms():
        if int(query_atom.GetAtomMapNum()) in ignored:
            continue
        target_atom = substrate.GetAtomWithIdx(int(match[int(query_atom.GetIdx())]))
        if not _chirality_matches_query(query_atom, target_atom, match):
            return False
    for query_bond in query.GetBonds():
        if not normalized_bond_stereo(query_bond):
            continue
        begin = int(match[int(query_bond.GetBeginAtomIdx())])
        end = int(match[int(query_bond.GetEndAtomIdx())])
        target_bond = substrate.GetBondBetweenAtoms(begin, end)
        if target_bond is None or not _bond_stereo_matches_query(query_bond, target_bond):
            return False
    return True


def _match_has_valid_attachments(
    rule: AriesRule,
    query: Chem.Mol,
    substrate: Chem.Mol,
    match: Sequence[int],
) -> bool:
    """Enforce the saturated boundary condition represented by ``Dn+``."""

    if not rule.attachment_ids:
        return True
    local_to_qidx = query_atom_maps(query)
    represented_target_bonds = {
        tuple(
            sorted(
                (
                    int(match[int(bond.GetBeginAtomIdx())]),
                    int(match[int(bond.GetEndAtomIdx())]),
                )
            )
        )
        for bond in query.GetBonds()
    }
    for local_id in rule.attachment_ids:
        query_idx = local_to_qidx.get(local_id)
        if query_idx is None:
            return False
        target_idx = int(match[query_idx])
        target_atom = substrate.GetAtomWithIdx(target_idx)
        if int(target_atom.GetDegree()) < int(query.GetAtomWithIdx(query_idx).GetDegree()):
            return False
        for bond in target_atom.GetBonds():
            pair = tuple(
                sorted((target_idx, int(bond.GetOtherAtomIdx(target_idx))))
            )
            if pair in represented_target_bonds:
                continue
            if bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
                return False
    return True


def _mapping_matches_query(
    query: Chem.Mol,
    substrate: Chem.Mol,
    match: Sequence[int],
    require_stereo: bool = True,
) -> bool:
    if len(match) != query.GetNumAtoms() or len(set(match)) != len(match):
        return False
    for query_atom in query.GetAtoms():
        target_atom = substrate.GetAtomWithIdx(
            int(match[int(query_atom.GetIdx())])
        )
        if not query_atom.Match(target_atom):
            return False
    for query_bond in query.GetBonds():
        target_bond = substrate.GetBondBetweenAtoms(
            int(match[int(query_bond.GetBeginAtomIdx())]),
            int(match[int(query_bond.GetEndAtomIdx())]),
        )
        if target_bond is None or not query_bond.Match(target_bond):
            return False
    return (
        _match_has_required_stereo(query, substrate, match)
        if require_stereo
        else True
    )


def _complete_observed_mapping(
    query: Chem.Mol,
    substrate: Chem.Mol,
    partial_match: Sequence[int | None],
    max_completions: int = 64,
) -> list[tuple[int, ...]]:
    """Complete a map-anchored embedding, including explicit unmapped H/D."""

    working = list(partial_match)
    used = {int(value) for value in working if value is not None}
    completions: list[tuple[int, ...]] = []

    def search() -> None:
        if len(completions) >= int(max_completions):
            return
        missing = [
            index for index, value in enumerate(working) if value is None
        ]
        if not missing:
            match = tuple(int(value) for value in working if value is not None)
            if _mapping_matches_query(
                query,
                substrate,
                match,
                require_stereo=False,
            ):
                completions.append(match)
            return

        query_idx = max(
            missing,
            key=lambda index: sum(
                working[int(neighbor.GetIdx())] is not None
                for neighbor in query.GetAtomWithIdx(index).GetNeighbors()
            ),
        )
        query_atom = query.GetAtomWithIdx(query_idx)
        assigned_neighbors = [
            int(neighbor.GetIdx())
            for neighbor in query_atom.GetNeighbors()
            if working[int(neighbor.GetIdx())] is not None
        ]
        if assigned_neighbors:
            anchor_target = int(working[assigned_neighbors[0]])
            candidates = [
                int(atom.GetIdx())
                for atom in substrate.GetAtomWithIdx(
                    anchor_target
                ).GetNeighbors()
            ]
        else:
            candidates = list(range(substrate.GetNumAtoms()))

        for target_idx in candidates:
            if target_idx in used:
                continue
            target_atom = substrate.GetAtomWithIdx(target_idx)
            if not query_atom.Match(target_atom):
                continue
            compatible = True
            for query_neighbor_idx in assigned_neighbors:
                query_bond = query.GetBondBetweenAtoms(
                    query_idx,
                    query_neighbor_idx,
                )
                target_bond = substrate.GetBondBetweenAtoms(
                    target_idx,
                    int(working[query_neighbor_idx]),
                )
                if (
                    query_bond is None
                    or target_bond is None
                    or not query_bond.Match(target_bond)
                ):
                    compatible = False
                    break
            if not compatible:
                continue
            working[query_idx] = target_idx
            used.add(target_idx)
            search()
            used.remove(target_idx)
            working[query_idx] = None

    search()
    return completions


def _observed_rule_matches(
    rule: AriesRule,
    query: Chem.Mol,
    substrate: Chem.Mol,
) -> list[tuple[int, ...]]:
    substrate_map_to_idx = {
        int(atom.GetAtomMapNum()): int(atom.GetIdx())
        for atom in substrate.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    matches = []
    for assignment in rule.source_assignments:
        local_to_map = dict(assignment)
        target_indices: list[int | None] = []
        for query_atom in query.GetAtoms():
            local_id = int(query_atom.GetAtomMapNum())
            atom_map = local_to_map.get(local_id)
            target_idx = substrate_map_to_idx.get(atom_map) if atom_map else None
            target_indices.append(target_idx)
        complete = all(value is not None for value in target_indices)
        match = tuple(
            int(value) for value in target_indices if value is not None
        )
        if complete and _mapping_matches_query(
            query,
            substrate,
            match,
            require_stereo=False,
        ):
            matches.append(match)
            continue

        matches.extend(
            _complete_observed_mapping(
                query,
                substrate,
                target_indices,
            )
        )
    return matches


def find_rule_matches(
    rule: AriesRule,
    substrate: Chem.Mol,
    rule_index: int = 0,
    max_matches: int = 1000,
) -> list[RuleMatch]:
    query = parse_aries_smarts(rule.reactant_smarts)
    if query is None:
        raise ValueError(f"could not parse ARIES reactant SMARTS: {rule.reactant_smarts}")
    local_to_qidx = query_atom_maps(query)
    center_qidx = [
        local_to_qidx[local_id]
        for local_id in rule.center_ids
        if local_id in local_to_qidx
    ]
    observed_match_list = _observed_rule_matches(rule, query, substrate)
    observed_matches = set(observed_match_list)
    relative_tetra_ids = {
        int(local_id)
        for local_id, action in rule.tetra_actions.items()
        if action in {"retain", "invert"}
    }
    raw_matches = list(observed_match_list)
    raw_matches.extend(
        substrate.GetSubstructMatches(
            query,
            useChirality=False,
            uniquify=False,
            maxMatches=int(max_matches),
        )
    )
    result = []
    seen_assignments: set[tuple[int, ...]] = set()
    for match in raw_matches:
        is_observed = tuple(match) in observed_matches
        if not _match_has_valid_attachments(
            rule,
            query,
            substrate,
            match,
        ):
            continue
        if not is_observed and not _match_has_required_stereo(
            query,
            substrate,
            match,
            ignored_tetra_local_ids=relative_tetra_ids,
        ):
            continue
        stereo_matches = True
        for spec in () if is_observed else rule.reactant_double_stereo:
            begin_qidx = local_to_qidx.get(spec.begin_id)
            end_qidx = local_to_qidx.get(spec.end_id)
            if begin_qidx is None or end_qidx is None:
                stereo_matches = False
                break
            target_bond = substrate.GetBondBetweenAtoms(
                int(match[begin_qidx]),
                int(match[end_qidx]),
            )
            if (
                target_bond is None
                or normalized_bond_stereo(target_bond) != spec.stereo
            ):
                stereo_matches = False
                break
        if not stereo_matches:
            continue
        center_assignment = tuple(
            int(match[local_to_qidx[local_id]])
            for local_id in sorted(rule.center_ids)
            if local_id in local_to_qidx
        )
        if center_assignment in seen_assignments:
            continue
        seen_assignments.add(center_assignment)
        center_indices = tuple(sorted({int(match[idx]) for idx in center_qidx}))
        atom_maps = tuple(
            sorted(
                {
                    int(substrate.GetAtomWithIdx(idx).GetAtomMapNum())
                    for idx in center_indices
                    if int(substrate.GetAtomWithIdx(idx).GetAtomMapNum()) > 0
                }
            )
        )
        result.append(
            RuleMatch(
                rule_index=rule_index,
                query_to_substrate=tuple(int(value) for value in match),
                center_atom_indices=center_indices,
                center_atom_maps=atom_maps,
            )
        )
    return result


def _compatible_subsets(
    matches: Sequence[RuleMatch],
    maximum_size: int,
    cap: int,
) -> list[tuple[RuleMatch, ...]]:
    output: list[tuple[RuleMatch, ...]] = []
    maximum_size = min(int(maximum_size), len(matches))
    # Full observed multiplicity comes first so the training product survives a cap.
    sizes = list(range(maximum_size, 0, -1))
    for size in sizes:
        for combination in itertools.combinations(matches, size):
            occupied: set[int] = set()
            compatible = True
            for match in combination:
                current = set(match.center_atom_indices)
                if occupied & current:
                    compatible = False
                    break
                occupied.update(current)
            if compatible:
                output.append(combination)
                if len(output) >= cap:
                    return output
    return output


def _application_combinations(
    template: AriesTemplate,
    matches_by_rule: Sequence[Sequence[RuleMatch]],
    cap: int,
) -> tuple[list[tuple[RuleMatch, ...]], bool]:
    choices = []
    for rule, matches in zip(template.rules, matches_by_rule):
        subsets = _compatible_subsets(
            matches,
            maximum_size=rule.observed_repeats,
            cap=cap,
        )
        if not subsets:
            return [], False
        choices.append(subsets)

    combinations: list[tuple[RuleMatch, ...]] = []
    truncated = False
    for grouped in itertools.product(*choices):
        flat = tuple(match for subset in grouped for match in subset)
        occupied: set[int] = set()
        compatible = True
        for match in flat:
            current = set(match.center_atom_indices)
            if occupied & current:
                compatible = False
                break
            occupied.update(current)
        if not compatible:
            continue
        combinations.append(flat)
        if len(combinations) >= cap:
            truncated = True
            break
    return combinations, truncated


def rule_is_component_coupling(rule: AriesRule) -> bool:
    """Return whether a changed product bond joins two left-side components."""

    reactant_query = parse_aries_smarts(rule.reactant_smarts)
    product_query = parse_aries_smarts(rule.product_smarts)
    if reactant_query is None or product_query is None:
        return False
    reactant_fragments = Chem.GetMolFrags(reactant_query, asMols=False)
    if len(reactant_fragments) < 2:
        return False
    qidx_to_component = {
        int(atom_idx): component_index
        for component_index, atom_indices in enumerate(reactant_fragments)
        for atom_idx in atom_indices
    }
    reactant_local_to_qidx = query_atom_maps(reactant_query)
    product_local_to_qidx = query_atom_maps(product_query)
    for begin_local, end_local in rule.changed_bond_ids:
        if (
            begin_local not in reactant_local_to_qidx
            or end_local not in reactant_local_to_qidx
            or begin_local not in product_local_to_qidx
            or end_local not in product_local_to_qidx
        ):
            continue
        if (
            qidx_to_component[reactant_local_to_qidx[begin_local]]
            == qidx_to_component[reactant_local_to_qidx[end_local]]
        ):
            continue
        product_bond = product_query.GetBondBetweenAtoms(
            product_local_to_qidx[begin_local],
            product_local_to_qidx[end_local],
        )
        if product_bond is not None:
            return True
    return False


def _atom_fragment_membership(mol: Chem.Mol) -> dict[int, int]:
    return {
        int(atom_idx): fragment_index
        for fragment_index, atom_indices in enumerate(
            Chem.GetMolFrags(mol, asMols=False)
        )
        for atom_idx in atom_indices
    }


def _combination_spans_input_fragments(
    template: AriesTemplate,
    combination: Sequence[RuleMatch],
    atom_to_fragment: dict[int, int],
) -> bool:
    for match in combination:
        if not rule_is_component_coupling(template.rules[match.rule_index]):
            continue
        fragments = {
            atom_to_fragment[int(atom_idx)]
            for atom_idx in match.query_to_substrate
        }
        if len(fragments) > 1:
            return True
    return False


def duplicated_mapped_substrate(mol: Chem.Mol, copies: int = 2) -> str:
    copies = max(2, int(copies))
    max_map = max(
        (int(atom.GetAtomMapNum()) for atom in mol.GetAtoms()),
        default=0,
    )
    stride = max(1000, max_map + 1000)
    combined: Chem.Mol | None = None
    for copy_index in range(copies):
        current = Chem.Mol(mol)
        for atom in current.GetAtoms():
            original_map = int(atom.GetAtomMapNum())
            atom.SetIntProp("_aries_copy_id", copy_index)
            atom.SetIntProp("_aries_original_map", original_map)
            if original_map > 0:
                atom.SetAtomMapNum(original_map + copy_index * stride)
        combined = current if combined is None else Chem.CombineMols(combined, current)
    if combined is None:
        return ""
    return Chem.MolToSmiles(combined, canonical=True, isomericSmiles=True)


def _bond_type(value: str) -> Chem.BondType:
    lookup = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
        "AROMATIC": Chem.BondType.AROMATIC,
        "QUADRUPLE": Chem.BondType.QUADRUPLE,
    }
    return lookup.get(value, Chem.BondType.SINGLE)


def _query_bond_type(bond: Chem.Bond) -> Chem.BondType:
    bond_type = bond.GetBondType()
    if bond_type == Chem.BondType.UNSPECIFIED:
        return Chem.BondType.SINGLE
    return bond_type


def _query_bond_is_aromatic(bond: Chem.Bond) -> bool:
    """Handle SMARTS aromatic bonds whose aromatic flag is not initialized."""

    return (
        bond.GetBondType() == Chem.BondType.AROMATIC
        or bool(bond.GetIsAromatic())
    )


def _set_product_atom_state(
    atom: Chem.Atom,
    state: AtomState,
) -> None:
    atom.SetAtomicNum(int(state.atomic_num))
    atom.SetIsAromatic(bool(state.aromatic))
    atom.SetFormalCharge(int(state.formal_charge))
    atom.SetNumRadicalElectrons(int(state.radical_electrons))
    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    if state.special_h is None:
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    else:
        atom.SetNumExplicitHs(int(state.special_h))
        atom.SetNoImplicit(True)


def _new_product_atom(state: AtomState) -> Chem.Atom:
    atom = Chem.Atom(int(state.atomic_num))
    _set_product_atom_state(atom, state)
    return atom


def _copy_product_tetrahedral_stereo(
    product_query: Chem.Mol,
    local_to_output: dict[int, int],
    output: Chem.Mol,
    owned_local_ids: set[int],
) -> None:
    output_to_local = {value: key for key, value in local_to_output.items()}
    for query_atom in product_query.GetAtoms():
        tag = query_atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        local_id = int(query_atom.GetAtomMapNum())
        if local_id not in owned_local_ids:
            continue
        output_idx = local_to_output.get(local_id)
        if output_idx is None:
            continue
        output_atom = output.GetAtomWithIdx(output_idx)
        if output_atom.GetDegree() < 3:
            continue
        query_labels = [
            int(neighbor.GetAtomMapNum())
            for neighbor in query_atom.GetNeighbors()
        ]
        output_labels = [
            output_to_local.get(int(neighbor.GetIdx()), 1000000 + int(neighbor.GetIdx()))
            for neighbor in output_atom.GetNeighbors()
        ]
        if len(query_labels) == 3:
            query_labels.append(-1)
        if len(output_labels) == 3:
            output_labels.append(-1)
        output_atom.SetChiralTag(tag)
        if (
            len(query_labels) == 4
            and len(output_labels) == 4
            and set(query_labels) == set(output_labels)
        ):
            rank = {value: index for index, value in enumerate(sorted(query_labels))}
            query_parity = parity([rank[value] for value in query_labels])
            output_parity = parity([rank[value] for value in output_labels])
            if query_parity != output_parity:
                output_atom.InvertChirality()


def _apply_tetrahedral_actions(
    actions: dict[int, str],
    local_to_source: dict[int, int],
    local_to_output: dict[int, int],
    substrate: Chem.Mol,
    output: Chem.Mol,
    source_to_output: dict[int, int] | None = None,
) -> set[int]:
    """Apply stereo as a relation to the matched substrate configuration."""

    source_to_output = source_to_output or {}
    source_idx_to_local = {
        int(atom_idx): int(local_id)
        for local_id, atom_idx in local_to_source.items()
    }
    output_idx_to_local = {
        int(atom_idx): int(local_id)
        for local_id, atom_idx in local_to_output.items()
    }

    def source_neighbor_label(atom: Chem.Atom) -> int:
        atom_idx = int(atom.GetIdx())
        if atom_idx in source_idx_to_local:
            return source_idx_to_local[atom_idx]
        if atom_idx in source_to_output:
            return 1000000 + int(source_to_output[atom_idx])
        atom_map = int(atom.GetAtomMapNum())
        if atom_map > 0:
            return 2000000 + atom_map
        return 3000000 + atom_idx

    def output_neighbor_label(atom: Chem.Atom) -> int:
        atom_idx = int(atom.GetIdx())
        if atom_idx in output_idx_to_local:
            return output_idx_to_local[atom_idx]
        return 1000000 + atom_idx

    handled: set[int] = set()
    for local_id, action in actions.items():
        output_idx = local_to_output.get(local_id)
        if output_idx is None:
            continue
        output_atom = output.GetAtomWithIdx(output_idx)
        if action == "erase":
            output_atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            handled.add(local_id)
            continue
        if action not in {"retain", "invert"}:
            continue
        source_idx = local_to_source.get(local_id)
        if source_idx is None:
            continue
        source_atom = substrate.GetAtomWithIdx(source_idx)
        source_tag = source_atom.GetChiralTag()
        if source_tag == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        output_atom.SetChiralTag(source_tag)
        observed_relation = tetra_relation_from_labels(
            source_tag,
            [source_neighbor_label(atom) for atom in source_atom.GetNeighbors()],
            output_atom.GetChiralTag(),
            [output_neighbor_label(atom) for atom in output_atom.GetNeighbors()],
        )
        if observed_relation not in {"retain", "invert"}:
            continue
        if observed_relation != action:
            output_atom.InvertChirality()
        handled.add(local_id)
    return handled


def _chiral_tag_from_text(value: str) -> Chem.ChiralType:
    if value == "CHI_TETRAHEDRAL_CW":
        return Chem.ChiralType.CHI_TETRAHEDRAL_CW
    if value == "CHI_TETRAHEDRAL_CCW":
        return Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    return Chem.ChiralType.CHI_UNSPECIFIED


def _source_remote_stereo_variant(
    template: AriesTemplate,
    substrate: Chem.Mol,
    product: Chem.Mol,
) -> Chem.Mol | None:
    """Apply source-only remote stereo provenance to the exact source graph."""

    if not template.source_remote_tetra_specs:
        return None
    if (
        canonical_connectivity_smiles_from_mol(
            substrate,
            keep_maps=True,
        )
        != template.source_mapped_substrate
    ):
        return None

    corrected = Chem.Mol(product)
    source_by_map = {
        int(atom.GetAtomMapNum()): int(atom.GetIdx())
        for atom in substrate.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    output_by_map = {
        int(atom.GetAtomMapNum()): int(atom.GetIdx())
        for atom in corrected.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    relative_actions = {
        spec.atom_map: spec.action
        for spec in template.source_remote_tetra_specs
        if spec.action in {"retain", "invert", "erase"}
        and spec.atom_map in output_by_map
    }
    _apply_tetrahedral_actions(
        relative_actions,
        source_by_map,
        output_by_map,
        substrate,
        corrected,
    )

    for spec in template.source_remote_tetra_specs:
        if spec.action != "create":
            continue
        output_idx = output_by_map.get(spec.atom_map)
        if output_idx is None:
            continue
        output_atom = corrected.GetAtomWithIdx(output_idx)
        tag = _chiral_tag_from_text(spec.product_tag)
        if tag == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        expected_labels = list(spec.product_neighbor_maps)
        output_labels = [
            int(neighbor.GetAtomMapNum())
            for neighbor in output_atom.GetNeighbors()
        ]
        if len(expected_labels) == 3:
            expected_labels.append(-1)
        if len(output_labels) == 3:
            output_labels.append(-1)
        output_atom.SetChiralTag(tag)
        if (
            len(expected_labels) == 4
            and len(output_labels) == 4
            and len(set(expected_labels)) == 4
            and set(expected_labels) == set(output_labels)
        ):
            rank = {
                value: index
                for index, value in enumerate(sorted(expected_labels))
            }
            expected_parity = parity(
                [rank[value] for value in expected_labels]
            )
            output_parity = parity([rank[value] for value in output_labels])
            if expected_parity != output_parity:
                output_atom.InvertChirality()

    try:
        Chem.AssignStereochemistry(corrected, cleanIt=False, force=True)
        Chem.SanitizeMol(corrected)
    except Exception:
        return None
    return corrected


def _copy_product_bond_stereo(
    specs: Sequence[DoubleBondStereoSpec],
    local_to_output: dict[int, int],
    output: Chem.Mol,
) -> None:
    for spec in specs:
        begin_output = local_to_output.get(spec.begin_id)
        end_output = local_to_output.get(spec.end_id)
        if begin_output is None or end_output is None:
            continue
        output_bond = output.GetBondBetweenAtoms(begin_output, end_output)
        if output_bond is None:
            continue
        if spec.stereo == "erase":
            output_bond.SetBondDir(Chem.BondDir.NONE)
            output_bond.SetStereo(Chem.BondStereo.STEREONONE)
            continue
        stereo_output = [
            local_to_output.get(spec.begin_neighbor_id),
            local_to_output.get(spec.end_neighbor_id),
        ]
        if any(value is None for value in stereo_output):
            continue
        first_stereo = int(stereo_output[0])
        second_stereo = int(stereo_output[1])
        if output_bond.GetBeginAtomIdx() != begin_output:
            first_stereo, second_stereo = second_stereo, first_stereo
        try:
            output_bond.SetStereoAtoms(first_stereo, second_stereo)
        except RuntimeError:
            continue
        output_bond.SetStereo(
            Chem.BondStereo.STEREOE
            if spec.stereo == "E"
            else Chem.BondStereo.STEREOZ
        )


def _capture_double_bond_stereo(
    mol: Chem.Mol,
) -> list[tuple[int, int, int, int, str]]:
    captured = []
    for bond in mol.GetBonds():
        stereo = normalized_bond_stereo(bond)
        stereo_atoms = list(bond.GetStereoAtoms())
        if stereo not in {"E", "Z"} or len(stereo_atoms) != 2:
            continue
        captured.append(
            (
                int(bond.GetBeginAtomIdx()),
                int(bond.GetEndAtomIdx()),
                int(stereo_atoms[0]),
                int(stereo_atoms[1]),
                stereo,
            )
        )
    return captured


def _clear_all_bond_stereo(mol: Chem.RWMol) -> None:
    for bond in mol.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)
        bond.SetStereo(Chem.BondStereo.STEREONONE)


def _restore_unchanged_double_bond_stereo(
    captured: Sequence[tuple[int, int, int, int, str]],
    old_to_new: dict[int, int],
    edited_old_pairs: set[tuple[int, int]],
    output: Chem.Mol,
) -> None:
    Chem.AssignStereochemistry(
        output,
        cleanIt=False,
        force=True,
        flagPossibleStereoCenters=True,
    )
    cip_ranks = [
        int(atom.GetProp("_CIPRank"))
        if atom.HasProp("_CIPRank")
        else int(atom.GetAtomicNum())
        for atom in output.GetAtoms()
    ]

    def restored_neighbor(
        endpoint: int,
        opposite: int,
    ) -> int | None:
        candidates = sorted(
            int(neighbor.GetIdx())
            for neighbor in output.GetAtomWithIdx(endpoint).GetNeighbors()
            if int(neighbor.GetIdx()) != opposite
        )
        if not candidates:
            return None
        return max(
            candidates,
            key=lambda value: (cip_ranks[value], value),
        )

    for old_begin, old_end, old_first, old_second, stereo in captured:
        if tuple(sorted((old_begin, old_end))) in edited_old_pairs:
            continue
        if old_begin not in old_to_new or old_end not in old_to_new:
            continue
        begin = old_to_new[old_begin]
        end = old_to_new[old_end]
        bond = output.GetBondBetweenAtoms(begin, end)
        if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        first = restored_neighbor(begin, end)
        second = restored_neighbor(end, begin)
        if first is None or second is None:
            continue
        if bond.GetBeginAtomIdx() != begin:
            first, second = second, first
        try:
            bond.SetStereoAtoms(first, second)
        except RuntimeError:
            continue
        bond.SetStereo(
            Chem.BondStereo.STEREOE
            if stereo == "E"
            else Chem.BondStereo.STEREOZ
        )


def _strip_maps(mol: Chem.Mol) -> Chem.Mol:
    result = Chem.Mol(mol)
    for atom in result.GetAtoms():
        atom.SetAtomMapNum(0)
    return result


def canonical_smiles_from_mol(mol: Chem.Mol, keep_maps: bool = False) -> str:
    result = Chem.Mol(mol)
    if not keep_maps:
        result = _strip_maps(result)
    try:
        Chem.SanitizeMol(result)
        Chem.AssignStereochemistry(result, cleanIt=True, force=True)
        return Chem.MolToSmiles(result, canonical=True, isomericSmiles=True)
    except Exception:
        return ""


def canonical_connectivity_smiles_from_mol(
    mol: Chem.Mol,
    keep_maps: bool = False,
) -> str:
    result = Chem.Mol(mol)
    Chem.RemoveStereochemistry(result)
    return canonical_smiles_from_mol(
        result,
        keep_maps=keep_maps,
    )


def canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return ""
    return canonical_smiles_from_mol(mol)


def _apply_match_combination(
    template: AriesTemplate,
    substrate: Chem.Mol,
    combination: Sequence[RuleMatch],
) -> tuple[Chem.Mol | None, str]:
    parsed: dict[int, tuple[Chem.Mol, Chem.Mol, dict[int, int], dict[int, int]]] = {}
    for rule_index in {match.rule_index for match in combination}:
        rule = template.rules[rule_index]
        reactant_query = parse_aries_smarts(rule.reactant_smarts)
        product_query = parse_aries_smarts(rule.product_smarts)
        if reactant_query is None or product_query is None:
            return None, "template query parse failed"
        parsed[rule_index] = (
            reactant_query,
            product_query,
            query_atom_maps(reactant_query),
            query_atom_maps(product_query),
        )

    delete_indices: set[int] = set()
    payload_root_indices: set[int] = set()
    protected_indices: set[int] = set()
    edited_old_pairs: set[tuple[int, int]] = set()
    for match in combination:
        rule = template.rules[match.rule_index]
        reactant_query, _, r_local_to_qidx, p_local_to_qidx = parsed[
            match.rule_index
        ]
        for local_id, query_idx in r_local_to_qidx.items():
            target_idx = int(match.query_to_substrate[query_idx])
            if local_id in p_local_to_qidx:
                protected_indices.add(target_idx)
                continue
            # Every explicitly represented reactant-only atom is deleted.
            # In FG mode this includes side-specific functional-group atoms
            # kept outside the minimal reacting-center definition.
            if (
                local_id not in rule.center_ids
                and (
                    template.mode != "fg"
                    or local_id not in rule.context_ids
                )
            ):
                continue
            delete_indices.add(target_idx)
            if (
                substrate.GetAtomWithIdx(target_idx).GetDegree()
                > reactant_query.GetAtomWithIdx(query_idx).GetDegree()
            ):
                payload_root_indices.add(target_idx)
        for begin_local, end_local in rule.changed_bond_ids:
            begin_qidx = r_local_to_qidx.get(begin_local)
            end_qidx = r_local_to_qidx.get(end_local)
            if begin_qidx is None or end_qidx is None:
                continue
            edited_old_pairs.add(
                tuple(
                    sorted(
                        (
                            int(match.query_to_substrate[begin_qidx]),
                            int(match.query_to_substrate[end_qidx]),
                        )
                    )
                )
            )

    # A compact leaving-group root can stand for a larger one-sided branch.
    # Remove that branch without crossing any atom retained by the product.
    frontier = list(payload_root_indices)
    while frontier:
        current_idx = frontier.pop()
        for neighbor in substrate.GetAtomWithIdx(current_idx).GetNeighbors():
            neighbor_idx = int(neighbor.GetIdx())
            if (
                neighbor_idx in protected_indices
                or neighbor_idx in delete_indices
            ):
                continue
            delete_indices.add(neighbor_idx)
            frontier.append(neighbor_idx)

    old_to_new: dict[int, int] = {}
    next_idx = 0
    for old_idx in range(substrate.GetNumAtoms()):
        if old_idx in delete_indices:
            continue
        old_to_new[old_idx] = next_idx
        next_idx += 1

    captured_double_stereo = _capture_double_bond_stereo(substrate)
    rw_mol = Chem.RWMol(substrate)
    _clear_all_bond_stereo(rw_mol)
    for atom_idx in sorted(delete_indices, reverse=True):
        rw_mol.RemoveAtom(atom_idx)

    stereo_jobs: list[
        tuple[
            Chem.Mol,
            tuple[DoubleBondStereoSpec, ...],
            dict[int, str],
            dict[int, int],
            dict[int, int],
            set[int],
        ]
    ] = []
    for match in combination:
        rule = template.rules[match.rule_index]
        reactant_query, product_query, r_local_to_qidx, p_local_to_qidx = parsed[
            match.rule_index
        ]
        local_to_output: dict[int, int] = {}
        local_to_source = {
            local_id: int(match.query_to_substrate[query_idx])
            for local_id, query_idx in r_local_to_qidx.items()
        }
        for local_id, query_idx in r_local_to_qidx.items():
            old_idx = int(match.query_to_substrate[query_idx])
            if old_idx in old_to_new and local_id in p_local_to_qidx:
                local_to_output[local_id] = old_to_new[old_idx]

        for local_id in sorted(p_local_to_qidx):
            if local_id in local_to_output:
                continue
            if local_id not in rule.center_ids:
                continue
            state = rule.product_states[local_id]
            new_idx = int(rw_mol.AddAtom(_new_product_atom(state)))
            local_to_output[local_id] = new_idx

        # Remove or replace bonds explicitly represented on the reactant side.
        product_bonds: dict[frozenset[int], Chem.Bond] = {}
        for bond in product_query.GetBonds():
            pair = frozenset(
                (
                    int(bond.GetBeginAtom().GetAtomMapNum()),
                    int(bond.GetEndAtom().GetAtomMapNum()),
                )
            )
            product_bonds[pair] = bond

        for bond in reactant_query.GetBonds():
            begin_local = int(bond.GetBeginAtom().GetAtomMapNum())
            end_local = int(bond.GetEndAtom().GetAtomMapNum())
            if tuple(sorted((begin_local, end_local))) not in rule.changed_bond_ids:
                continue
            if begin_local not in local_to_output or end_local not in local_to_output:
                continue
            pair = frozenset((begin_local, end_local))
            product_bond = product_bonds.get(pair)
            begin_output = local_to_output[begin_local]
            end_output = local_to_output[end_local]
            existing = rw_mol.GetBondBetweenAtoms(begin_output, end_output)
            if existing is None:
                continue
            if (
                product_bond is None
                or _query_bond_type(existing) != _query_bond_type(product_bond)
                or bool(existing.GetIsAromatic())
                != _query_bond_is_aromatic(product_bond)
            ):
                rw_mol.RemoveBond(begin_output, end_output)

        for bond in product_query.GetBonds():
            begin_local = int(bond.GetBeginAtom().GetAtomMapNum())
            end_local = int(bond.GetEndAtom().GetAtomMapNum())
            if tuple(sorted((begin_local, end_local))) not in rule.changed_bond_ids:
                continue
            begin_output = local_to_output[begin_local]
            end_output = local_to_output[end_local]
            desired_type = _query_bond_type(bond)
            existing = rw_mol.GetBondBetweenAtoms(begin_output, end_output)
            if existing is None:
                rw_mol.AddBond(begin_output, end_output, desired_type)
                existing = rw_mol.GetBondBetweenAtoms(begin_output, end_output)
            elif (
                _query_bond_type(existing) != desired_type
                or bool(existing.GetIsAromatic()) != _query_bond_is_aromatic(bond)
            ):
                rw_mol.RemoveBond(begin_output, end_output)
                rw_mol.AddBond(begin_output, end_output, desired_type)
                existing = rw_mol.GetBondBetweenAtoms(begin_output, end_output)
            if existing is not None:
                existing.SetIsAromatic(_query_bond_is_aromatic(bond))
                existing.SetBondDir(Chem.BondDir.NONE)
                existing.SetStereo(Chem.BondStereo.STEREONONE)

        for local_id, output_idx in local_to_output.items():
            if local_id not in rule.center_ids:
                continue
            state = rule.product_states.get(local_id)
            if state is not None:
                _set_product_atom_state(rw_mol.GetAtomWithIdx(output_idx), state)

        stereo_jobs.append(
            (
                product_query,
                rule.product_double_stereo,
                rule.tetra_actions,
                local_to_source,
                local_to_output,
                set(rule.center_ids),
            )
        )

    try:
        product = rw_mol.GetMol()
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
    except Exception as exc:
        return None, f"pre-stereo {type(exc).__name__}: {exc}"

    _restore_unchanged_double_bond_stereo(
        captured_double_stereo,
        old_to_new,
        edited_old_pairs,
        product,
    )
    for (
        product_query,
        product_double_stereo,
        tetra_actions,
        local_to_source,
        local_to_output,
        center_ids,
    ) in stereo_jobs:
        effective_tetra_actions = dict(tetra_actions)
        for local_id in center_ids:
            if local_id in effective_tetra_actions:
                continue
            source_idx = local_to_source.get(local_id)
            output_idx = local_to_output.get(local_id)
            if source_idx is None or output_idx is None:
                continue
            source_atom = substrate.GetAtomWithIdx(source_idx)
            if (
                source_atom.GetChiralTag()
                != Chem.ChiralType.CHI_UNSPECIFIED
            ):
                # Retaining matched substrate configuration is the default
                # applicator behavior, not a template stereo constraint.
                effective_tetra_actions[local_id] = "retain"
        handled_tetra = _apply_tetrahedral_actions(
            effective_tetra_actions,
            local_to_source,
            local_to_output,
            substrate,
            product,
            source_to_output=old_to_new,
        )
        _copy_product_tetrahedral_stereo(
            product_query,
            local_to_output,
            product,
            center_ids - handled_tetra,
        )
        _copy_product_bond_stereo(
            product_double_stereo,
            local_to_output,
            product,
        )

    try:
        Chem.SetDoubleBondNeighborDirections(product)
        Chem.AssignStereochemistry(product, cleanIt=False, force=True)
        Chem.SanitizeMol(product)
    except Exception as exc:
        return None, f"post-stereo {type(exc).__name__}: {exc}"
    return product, ""


def apply_aries_template(
    template: AriesTemplate | str,
    mapped_substrate_smiles: str,
    max_products: int = 1000,
    max_matches_per_rule: int = 1000,
    allow_intramolecular: bool = True,
    allow_inter_copy: bool = True,
    _require_interfragment: bool = False,
) -> ApplicationResult:
    if isinstance(template, str):
        template = AriesTemplate.from_json(template)
    result = ApplicationResult()
    if template.status != "ok":
        result.error = template.error or "template extraction failed"
        return result
    if not template.rules:
        result.error = "template contains no rules"
        return result

    try:
        substrate = parse_mol(mapped_substrate_smiles, "mapped substrate")
        matches_by_rule = [
            find_rule_matches(
                rule,
                substrate,
                rule_index=rule_index,
                max_matches=max_matches_per_rule,
            )
            for rule_index, rule in enumerate(template.rules)
        ]
    except Exception as exc:
        result.error = f"{type(exc).__name__}: {exc}"
        return result

    result.match_count = sum(len(matches) for matches in matches_by_rule)
    for matches in matches_by_rule:
        for match in matches:
            if match.center_atom_maps:
                result.predicted_centers.add(match.center_atom_maps)
    combinations, truncated = _application_combinations(
        template,
        matches_by_rule,
        cap=max(int(max_products) * 4, int(max_products)),
    )
    result.truncated = truncated
    if not combinations:
        result.error = "one or more rules had no compatible match"
        return result

    atom_to_fragment = _atom_fragment_membership(substrate)
    for combination in combinations:
        if len(result.products) >= int(max_products):
            result.truncated = True
            break
        spans_fragments = _combination_spans_input_fragments(
            template,
            combination,
            atom_to_fragment,
        )
        if _require_interfragment and not spans_fragments:
            continue
        if not _require_interfragment and not spans_fragments and not allow_intramolecular:
            continue
        result.attempted_combinations += 1
        center_maps = tuple(
            sorted(
                {
                    atom_map
                    for match in combination
                    for atom_map in match.center_atom_maps
                }
            )
        )
        if center_maps:
            result.predicted_centers.add(center_maps)
        product, invalid_error = _apply_match_combination(
            template,
            substrate,
            combination,
        )
        if product is None:
            result.invalid_product_count += 1
            result.invalid_product_errors[invalid_error] = (
                result.invalid_product_errors.get(invalid_error, 0) + 1
            )
            continue
        product_variants = [product]
        source_stereo_variant = _source_remote_stereo_variant(
            template,
            substrate,
            product,
        )
        if source_stereo_variant is not None:
            product_variants.append(source_stereo_variant)
        for product_variant in product_variants:
            product_smiles = canonical_smiles_from_mol(product_variant)
            mapped_product = canonical_smiles_from_mol(
                product_variant,
                keep_maps=True,
            )
            if product_smiles:
                result.products.add(product_smiles)
                if center_maps:
                    result.product_centers.setdefault(product_smiles, set()).add(
                        center_maps
                    )
                if spans_fragments:
                    result.intermolecular_products.add(product_smiles)
                else:
                    result.intramolecular_products.add(product_smiles)
            if mapped_product:
                result.mapped_products.add(mapped_product)

    has_coupling_rule = any(rule_is_component_coupling(rule) for rule in template.rules)
    substrate_fragment_count = len(Chem.GetMolFrags(substrate, asMols=False))
    if (
        allow_inter_copy
        and has_coupling_rule
        and substrate_fragment_count == 1
        and len(result.products) < int(max_products)
    ):
        duplicated = duplicated_mapped_substrate(substrate, copies=2)
        inter_result = apply_aries_template(
            template,
            duplicated,
            max_products=max(1, int(max_products) - len(result.products)),
            max_matches_per_rule=max_matches_per_rule,
            allow_intramolecular=False,
            allow_inter_copy=False,
            _require_interfragment=True,
        )
        result.products.update(inter_result.products)
        result.intermolecular_products.update(inter_result.products)
        result.mapped_products.update(inter_result.mapped_products)
        for product_smiles, centers in inter_result.product_centers.items():
            result.product_centers.setdefault(product_smiles, set()).update(centers)
        result.predicted_centers.update(inter_result.predicted_centers)
        result.match_count += inter_result.match_count
        result.attempted_combinations += inter_result.attempted_combinations
        result.invalid_product_count += inter_result.invalid_product_count
        for error, count in inter_result.invalid_product_errors.items():
            result.invalid_product_errors[error] = (
                result.invalid_product_errors.get(error, 0) + count
            )
        result.truncated = result.truncated or inter_result.truncated
        if not result.error and inter_result.error and not result.products:
            result.error = inter_result.error
    return result


def template_atom_count(template: AriesTemplate) -> int:
    return sum(
        len(re.findall(r"\[[^\]]+\]", rule.reactant_smarts))
        + len(re.findall(r"\[[^\]]+\]", rule.product_smarts))
        for rule in template.rules
    )


def template_char_count(template: AriesTemplate) -> int:
    return len(template.compact_smarts)
