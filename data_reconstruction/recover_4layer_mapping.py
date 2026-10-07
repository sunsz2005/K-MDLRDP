#!/usr/bin/env python3
"""Recover the 4-layer numeric-to-protein mapping from three named source graphs.

The published instance contains numeric vertex identifiers only.  STRING v9.1,
WI-PHI and BioGRID are available with systematic yeast ORF names, and the
4-layer instance is an induced subgraph of their common vertex set.  This
script recovers the lost permutation by exact multilayer graph matching and
then verifies every edge, not merely graph sizes or degree sequences.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_INSTANCES = ROOT.parent.parent / "k-MLRDP-main" / "Instances" / "PPI" / "4Layers"
DEFAULT_SOURCES = ROOT / "sources" / "SupportnessOfPPIComplexes" / "networks"
LAYER_NAMES = ("STRING", "wiphi", "BioGRID")


def normalise_name(value: str) -> str:
    value = value.strip()
    return value[5:] if value.startswith("4932.") else value


def read_named_graph(path: Path, string_format: bool = False):
    nodes: set[str] = set()
    edges: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8") as handle:
        if string_format:
            next(handle)
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            u, v = normalise_name(fields[0]), normalise_name(fields[1])
            nodes.update((u, v))
            if u != v:
                edges.add((u, v) if u < v else (v, u))
    return nodes, edges


def read_numeric_graph(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        first = handle.readline().split()
        n, declared_m = map(int, first[:2])
        edges: set[tuple[int, int]] = set()
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            u, v = map(int, fields[:2])
            if u != v:
                edges.add((u, v) if u < v else (v, u))
    if len(edges) != declared_m:
        raise ValueError(f"{path}: declared {declared_m} edges but read {len(edges)}")
    return n, edges


def adjacency(nodes, edges):
    result = {u: set() for u in nodes}
    for u, v in edges:
        result[u].add(v)
        result[v].add(u)
    return result


def degree_signatures(nodes, adjs):
    return {u: tuple(len(adj[u]) for adj in adjs) for u in nodes}


def unique_exact_degree_proposals(target_nodes, candidate_nodes, target_adjs, candidate_adjs):
    target_by_sig = defaultdict(list)
    candidate_by_sig = defaultdict(list)
    for node, sig in degree_signatures(target_nodes, target_adjs).items():
        target_by_sig[sig].append(node)
    for node, sig in degree_signatures(candidate_nodes, candidate_adjs).items():
        candidate_by_sig[sig].append(node)
    proposals = []
    for sig in target_by_sig.keys() & candidate_by_sig.keys():
        if len(target_by_sig[sig]) == len(candidate_by_sig[sig]) == 1:
            proposals.append((target_by_sig[sig][0], candidate_by_sig[sig][0]))
    return proposals


def compatible(left, right, target_adjs, candidate_adjs):
    t1, c1 = left
    t2, c2 = right
    return all((t2 in ta[t1]) == (c2 in ca[c1]) for ta, ca in zip(target_adjs, candidate_adjs))


def prune_to_consistent_seed(proposals, target_adjs, candidate_adjs):
    """Find the large mutually consistent core among unique degree matches.

    Correct proposals form a clique because deleting vertices never changes
    adjacency between retained vertices.  False degree coincidences normally
    conflict with many true proposals.  Repeatedly remove the proposal with
    most conflicts; deterministic tie-breaking keeps runs reproducible.
    """
    active = set(range(len(proposals)))
    conflicts = [set() for _ in proposals]
    for i in range(len(proposals)):
        for j in range(i + 1, len(proposals)):
            if not compatible(proposals[i], proposals[j], target_adjs, candidate_adjs):
                conflicts[i].add(j)
                conflicts[j].add(i)

    removed = 0
    while True:
        worst = max(active, key=lambda i: (len(conflicts[i] & active), i))
        count = len(conflicts[worst] & active)
        if count == 0:
            break
        active.remove(worst)
        removed += 1
    seed = [proposals[i] for i in sorted(active)]
    return seed, removed


def fingerprint(node, anchors, adjs, target_side: bool):
    # Tuples of anchor positions are compact for the sparse WI-PHI/BioGRID
    # layers and remain collision-free (unlike a hash).
    if target_side:
        return tuple(tuple(i for i, (t, _) in enumerate(anchors) if t in adj[node]) for adj in adjs)
    return tuple(tuple(i for i, (_, c) in enumerate(anchors) if c in adj[node]) for adj in adjs)


def expand_mapping(seed, target_nodes, candidate_nodes, target_adjs, candidate_adjs):
    mapping = dict(seed)
    reverse = {candidate: target for target, candidate in seed}
    round_stats = []
    while True:
        anchors = sorted(mapping.items())
        unmatched_target = target_nodes - mapping.keys()
        unmatched_candidate = candidate_nodes - reverse.keys()
        target_groups = defaultdict(list)
        candidate_groups = defaultdict(list)
        for node in unmatched_target:
            target_groups[fingerprint(node, anchors, target_adjs, True)].append(node)
        for node in unmatched_candidate:
            candidate_groups[fingerprint(node, anchors, candidate_adjs, False)].append(node)

        additions = []
        for sig in target_groups.keys() & candidate_groups.keys():
            ts, cs = target_groups[sig], candidate_groups[sig]
            if len(ts) == len(cs) == 1:
                additions.append((ts[0], cs[0]))
        if not additions:
            break
        for target, candidate in additions:
            mapping[target] = candidate
            reverse[candidate] = target
        round_stats.append(len(additions))
    return mapping, round_stats


def mapped_edges(edges, mapping):
    result = set()
    for u, v in edges:
        a, b = mapping[u], mapping[v]
        result.add((a, b) if a < b else (b, a))
    return result


def verify(mapping, target_edges, source_edges):
    selected = set(mapping.values())
    reports = []
    all_exact = True
    for layer, t_edges, s_edges in zip(LAYER_NAMES, target_edges, source_edges):
        expected = {(u, v) for u, v in s_edges if u in selected and v in selected}
        observed = mapped_edges(t_edges, mapping)
        missing = expected - observed
        extra = observed - expected
        exact = not missing and not extra
        all_exact &= exact
        reports.append(
            {
                "layer": layer,
                "expected_edges": len(expected),
                "observed_edges": len(observed),
                "missing_edges": len(missing),
                "extra_edges": len(extra),
                "exact": exact,
                "missing_examples": sorted(missing)[:5],
                "extra_examples": sorted(extra)[:5],
            }
        )
    return all_exact, reports


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", type=Path, default=DEFAULT_INSTANCES)
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "4Layers")
    args = parser.parse_args()

    source_paths = (
        args.sources / "String" / "4932.protein.links.v9.1.txt",
        args.sources / "wiphi.txt",
        args.sources / "BioGRID.txt",
    )
    target_paths = tuple(args.instances / f"{name}_mapped_names.txt" for name in LAYER_NAMES)

    named = [read_named_graph(path, i == 0) for i, path in enumerate(source_paths)]
    candidate_nodes = set.intersection(*(nodes for nodes, _ in named))
    source_edges = [edges for _, edges in named]
    candidate_edges = [{edge for edge in edges if edge[0] in candidate_nodes and edge[1] in candidate_nodes} for edges in source_edges]
    candidate_adjs = [adjacency(candidate_nodes, edges) for edges in candidate_edges]

    numeric = [read_numeric_graph(path) for path in target_paths]
    sizes = {n for n, _ in numeric}
    if len(sizes) != 1:
        raise ValueError(f"target layers disagree on vertex count: {sizes}")
    n = sizes.pop()
    target_nodes = set(range(1, n + 1))
    target_edges = [edges for _, edges in numeric]
    target_adjs = [adjacency(target_nodes, edges) for edges in target_edges]

    print(f"candidate common nodes: {len(candidate_nodes)}; target nodes: {n}")
    proposals = unique_exact_degree_proposals(target_nodes, candidate_nodes, target_adjs, candidate_adjs)
    print(f"unique exact degree proposals: {len(proposals)}")
    seed, removed = prune_to_consistent_seed(proposals, target_adjs, candidate_adjs)
    print(f"mutually consistent seed: {len(seed)} (removed {removed})")
    mapping, rounds = expand_mapping(seed, target_nodes, candidate_nodes, target_adjs, candidate_adjs)
    print(f"fingerprint expansion rounds: {rounds}; mapped: {len(mapping)}/{n}")

    exact = False
    reports = []
    if len(mapping) == n:
        exact, reports = verify(mapping, target_edges, source_edges)
    print(f"exact edge-set verification: {exact}")
    for report in reports:
        print(report)

    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "node_mapping.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("internal_id", "systematic_name", "method"))
        for internal_id, name in sorted(mapping.items()):
            writer.writerow((internal_id, name, "exact_multilayer_graph_match"))
    report = {
        "candidate_nodes": len(candidate_nodes),
        "target_nodes": n,
        "initial_proposals": len(proposals),
        "consistent_seed": len(seed),
        "expansion_rounds": rounds,
        "mapped_nodes": len(mapping),
        "exact_all_layers": exact,
        "layers": reports,
    }
    (args.output / "validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not exact:
        raise SystemExit("mapping is incomplete or failed exact verification")


if __name__ == "__main__":
    main()
