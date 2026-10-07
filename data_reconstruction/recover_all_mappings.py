#!/usr/bin/env python3
"""Recover protein identifiers for every published 2--5 layer PPI instance.

Run ``recover_4layer_mapping.py`` first.  The published numeric identifiers are
nested: the m+1-layer instance uses a prefix of the identifiers in the m-layer
instance.  We verify that property edge by edge, extend the recovered 4-layer
mapping upward to 3 and 2 layers by exact adjacency fingerprints, and restrict
it downward to 5 layers.  Finally, all five numeric layers are translated back
to systematic yeast ORF names and canonical numeric edge lists are rebuilt.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from recover_4layer_mapping import (
    ROOT,
    adjacency,
    mapped_edges,
    read_named_graph,
    read_numeric_graph,
)


REFERENCE = ROOT.parent.parent / "k-MLRDP-main" / "Instances" / "PPI"
SOURCES = ROOT / "sources" / "SupportnessOfPPIComplexes" / "networks"
OUTPUT = ROOT / "output"
ORDER = ("STRING", "iRefIndex", "wiphi", "BioGRID", "MINT")
LAYERS = {
    2: ORDER[:2],
    3: ORDER[:3],
    4: ORDER[:4],
    5: ORDER,
}


def load_mapping(path: Path):
    result = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[int(row["internal_id"])] = row["systematic_name"]
    return result


def read_sources():
    string_nodes, string_edges = read_named_graph(
        SOURCES / "String" / "4932.protein.links.v9.1.txt", True
    )
    wiphi_nodes, wiphi_edges = read_named_graph(SOURCES / "wiphi.txt")
    biogrid_nodes, biogrid_edges = read_named_graph(SOURCES / "BioGRID.txt")
    return {
        "STRING": (string_nodes, string_edges),
        "wiphi": (wiphi_nodes, wiphi_edges),
        "BioGRID": (biogrid_nodes, biogrid_edges),
    }


def target_graph(layer_count: int, name: str):
    return read_numeric_graph(REFERENCE / f"{layer_count}Layers" / f"{name}_mapped_names.txt")


def fp(node, anchor_index, adjs):
    # Iterate over the node's neighbours rather than over thousands of anchors.
    # This changes the extension from O(|V|*|anchors|) to O(|E|).
    return tuple(
        tuple(sorted(anchor_index[neighbour] for neighbour in adj[node] if neighbour in anchor_index))
        for adj in adjs
    )


def extend_prefix_mapping(initial, target_n, candidate_nodes, target_adjs, candidate_adjs):
    mapping = dict(initial)
    reverse = {name: node for node, name in mapping.items()}
    stats = []
    while True:
        anchors = sorted(mapping.items())
        target_anchor_index = {node: i for i, (node, _) in enumerate(anchors)}
        candidate_anchor_index = {name: i for i, (_, name) in enumerate(anchors)}
        target_groups = defaultdict(list)
        candidate_groups = defaultdict(list)
        for node in set(range(1, target_n + 1)) - mapping.keys():
            target_groups[fp(node, target_anchor_index, target_adjs)].append(node)
        for node in candidate_nodes - reverse.keys():
            candidate_groups[fp(node, candidate_anchor_index, candidate_adjs)].append(node)
        additions = []
        for signature in target_groups.keys() & candidate_groups.keys():
            left, right = target_groups[signature], candidate_groups[signature]
            if len(left) == len(right) == 1:
                additions.append((left[0], right[0]))
        if not additions:
            break
        for node, name in additions:
            mapping[node] = name
            reverse[name] = node
        stats.append(len(additions))
    return mapping, stats


def validate_named_layers(layer_count, mapping, source_data, names):
    selected = set(mapping.values())
    results = []
    exact = True
    for name in names:
        _, target_edges = target_graph(layer_count, name)
        source_edges = source_data[name][1]
        expected = {(u, v) for u, v in source_edges if u in selected and v in selected}
        observed = mapped_edges(target_edges, mapping)
        missing, extra = expected - observed, observed - expected
        layer_exact = not missing and not extra
        exact &= layer_exact
        results.append(
            {
                "layer": name,
                "expected_edges": len(expected),
                "observed_edges": len(observed),
                "missing_edges": len(missing),
                "extra_edges": len(extra),
                "exact": layer_exact,
            }
        )
    return exact, results


def nesting_report():
    reports = []
    for low, high in ((2, 3), (3, 4), (4, 5)):
        common = LAYERS[low]
        for name in common:
            _, low_edges = target_graph(low, name)
            high_n, high_edges = target_graph(high, name)
            induced = {(u, v) for u, v in low_edges if u <= high_n and v <= high_n}
            reports.append(
                {
                    "larger_vertex_set": f"{low}Layers",
                    "prefix_vertex_set": f"{high}Layers",
                    "layer": name,
                    "edge_subset_exact": high_edges <= low_edges,
                    "prefix_induced_exact": high_edges == induced,
                }
            )
    return reports


def write_mapping(layer_count, mapping, method):
    folder = OUTPUT / f"{layer_count}Layers"
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "node_mapping.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("internal_id", "systematic_name", "method"))
        for node, name in sorted(mapping.items()):
            writer.writerow((node, name, method))


def write_reconstructed_layers(layer_count, mapping):
    complete = len(mapping) == target_graph(layer_count, "STRING")[0]
    named_dir = OUTPUT / f"{layer_count}Layers" / ("named_edges" if complete else "named_edges_partial")
    rebuilt_dir = OUTPUT / f"{layer_count}Layers" / "rebuilt_numeric"
    named_dir.mkdir(parents=True, exist_ok=True)
    rebuilt_dir.mkdir(parents=True, exist_ok=True)
    for layer in LAYERS[layer_count]:
        n, edges = target_graph(layer_count, layer)
        mappable_edges = {(u, v) for u, v in edges if u in mapping and v in mapping}
        named = sorted(mapped_edges(mappable_edges, mapping))
        with (named_dir / f"{layer}.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("protein_a", "protein_b"))
            writer.writerows(named)
        with (rebuilt_dir / f"{layer}_mapped_names.txt").open("w", encoding="utf-8", newline="") as handle:
            handle.write(f"{n}\t{len(edges)}\n")
            for u, v in sorted(edges):
                handle.write(f"{u}\t{v}\n")


def unresolved_groups(mapping, target_n, candidate_nodes, target_adjs, candidate_adjs):
    anchors = sorted(mapping.items())
    target_anchor_index = {node: i for i, (node, _) in enumerate(anchors)}
    candidate_anchor_index = {name: i for i, (_, name) in enumerate(anchors)}
    target_groups = defaultdict(list)
    candidate_groups = defaultdict(list)
    for node in set(range(1, target_n + 1)) - mapping.keys():
        target_groups[fp(node, target_anchor_index, target_adjs)].append(node)
    for name in candidate_nodes - set(mapping.values()):
        candidate_groups[fp(name, candidate_anchor_index, candidate_adjs)].append(name)
    result = []
    for signature, ids in target_groups.items():
        result.append({"internal_ids": sorted(ids), "candidate_names": sorted(candidate_groups.get(signature, []))})
    return sorted(result, key=lambda item: item["internal_ids"])


def main():
    source_data = read_sources()
    mapping4 = load_mapping(OUTPUT / "4Layers" / "node_mapping.tsv")
    if len(mapping4) != 5237:
        raise SystemExit("complete 4-layer mapping is required; run recover_4layer_mapping.py first")

    # The exact edge-subset test proves that identifiers are nested prefixes.
    nested = nesting_report()
    if not all(r["edge_subset_exact"] and r["prefix_induced_exact"] for r in nested):
        raise SystemExit("published identifiers are not nested as expected")

    # 3 layers: candidate nodes are the STRING/WI-PHI intersection.
    candidate3 = source_data["STRING"][0] & source_data["wiphi"][0]
    source_adjs3 = [
        adjacency(candidate3, {e for e in source_data[name][1] if e[0] in candidate3 and e[1] in candidate3})
        for name in ("STRING", "wiphi")
    ]
    target_adjs3 = []
    for name in ("STRING", "wiphi"):
        n3, edges = target_graph(3, name)
        target_adjs3.append(adjacency(set(range(1, n3 + 1)), edges))
    mapping3, rounds3 = extend_prefix_mapping(mapping4, n3, candidate3, target_adjs3, source_adjs3)

    # 2 layers: candidate nodes are all STRING v9.1 proteins.
    candidate2 = source_data["STRING"][0]
    source_adj2 = [adjacency(candidate2, source_data["STRING"][1])]
    n2, edges2 = target_graph(2, "STRING")
    target_adj2 = [adjacency(set(range(1, n2 + 1)), edges2)]
    mapping2, rounds2 = extend_prefix_mapping(mapping3, n2, candidate2, target_adj2, source_adj2)
    unresolved2 = unresolved_groups(mapping2, n2, candidate2, target_adj2, source_adj2)

    # 5 layers use exactly the first 4,223 identifiers from the 4-layer map.
    mapping5 = {node: mapping4[node] for node in range(1, 4224)}
    mappings = {2: mapping2, 3: mapping3, 4: mapping4, 5: mapping5}
    expected_sizes = {3: 5708, 4: 5237, 5: 4223}
    for count, expected in expected_sizes.items():
        if len(mappings[count]) != expected:
            raise SystemExit(f"{count}-layer mapping incomplete: {len(mappings[count])}/{expected}")

    validations = {}
    known_layers = {3: ("STRING", "wiphi"), 4: ("STRING", "wiphi", "BioGRID"), 5: ("STRING", "wiphi", "BioGRID")}
    for count, names in known_layers.items():
        exact, layers = validate_named_layers(count, mappings[count], source_data, names)
        validations[count] = {"exact": exact, "layers": layers}
        if not exact:
            raise SystemExit(f"{count}-layer mapping failed exact edge validation")

    for count, mapping in mappings.items():
        method = "exact_multilayer_graph_match" if count == 4 else "verified_nested_prefix_and_exact_adjacency_fingerprint"
        write_mapping(count, mapping, method)
        write_reconstructed_layers(count, mapping)

    unresolved_path = OUTPUT / "2Layers" / "unresolved_candidates.tsv"
    with unresolved_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("internal_ids", "candidate_systematic_names", "reason"))
        for group in unresolved2:
            writer.writerow(
                (
                    ",".join(map(str, group["internal_ids"])),
                    ",".join(group["candidate_names"]),
                    "indistinguishable_in_STRING; historical_named_iRefIndex_required",
                )
            )

    report = {
        "source_versions": {
            "STRING": "v9.1 (exact size and edge fingerprint)",
            "WI-PHI": "author repository file wiphi.txt (exact size and edge fingerprint)",
            "BioGRID": "author repository file BioGRID.txt (exact size and edge fingerprint)",
            "iRefIndex": "release not stated; restricted named edges recovered from published numeric graph",
            "MINT": "release not stated; restricted named edges recovered from published numeric graph",
        },
        "mapping_sizes": {str(k): len(v) for k, v in mappings.items()},
        "mapping_status": {
            "2Layers": f"partial: {len(mapping2)}/{n2}; {n2-len(mapping2)} identifiers unresolved",
            "3Layers": "complete and exact",
            "4Layers": "complete and exact",
            "5Layers": "complete and exact",
        },
        "2layer_unresolved_groups": unresolved2,
        "extension_rounds": {"3Layers": rounds3, "2Layers": rounds2},
        "nesting": nested,
        "known_source_edge_validation": validations,
    }
    (OUTPUT / "reconstruction_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
