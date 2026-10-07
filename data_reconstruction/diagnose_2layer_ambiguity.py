#!/usr/bin/env python3
"""Inspect the few 2-layer vertices that STRING alone cannot identify."""

from __future__ import annotations

from collections import defaultdict

from recover_4layer_mapping import adjacency
from recover_all_mappings import (
    OUTPUT,
    extend_prefix_mapping,
    fp,
    load_mapping,
    read_sources,
    target_graph,
)


def main():
    sources = read_sources()
    mapping4 = load_mapping(OUTPUT / "4Layers" / "node_mapping.tsv")

    candidate3 = sources["STRING"][0] & sources["wiphi"][0]
    source_adjs3 = [
        adjacency(candidate3, {e for e in sources[name][1] if e[0] in candidate3 and e[1] in candidate3})
        for name in ("STRING", "wiphi")
    ]
    target_adjs3 = []
    for name in ("STRING", "wiphi"):
        n3, edges = target_graph(3, name)
        target_adjs3.append(adjacency(set(range(1, n3 + 1)), edges))
    mapping3, _ = extend_prefix_mapping(mapping4, n3, candidate3, target_adjs3, source_adjs3)

    candidate2 = sources["STRING"][0]
    source_adj2 = [adjacency(candidate2, sources["STRING"][1])]
    n2, edges2 = target_graph(2, "STRING")
    target_adj2 = [adjacency(set(range(1, n2 + 1)), edges2)]
    mapping2, _ = extend_prefix_mapping(mapping3, n2, candidate2, target_adj2, source_adj2)

    anchors = sorted(mapping2.items())
    tidx = {node: i for i, (node, _) in enumerate(anchors)}
    cidx = {name: i for i, (_, name) in enumerate(anchors)}
    target_groups = defaultdict(list)
    candidate_groups = defaultdict(list)
    for node in set(range(1, n2 + 1)) - mapping2.keys():
        target_groups[fp(node, tidx, target_adj2)].append(node)
    for name in candidate2 - set(mapping2.values()):
        candidate_groups[fp(name, cidx, source_adj2)].append(name)

    _, iref_edges = target_graph(2, "iRefIndex")
    iref_adj = adjacency(set(range(1, n2 + 1)), iref_edges)
    bg_adj = adjacency(sources["BioGRID"][0], sources["BioGRID"][1])
    wi_adj = adjacency(sources["wiphi"][0], sources["wiphi"][1])
    mapped_names = set(mapping2.values())

    print(f"mapped={len(mapping2)} unresolved={n2-len(mapping2)}")
    for signature, target_ids in sorted(target_groups.items(), key=lambda item: item[1]):
        candidates = sorted(candidate_groups.get(signature, []))
        print(f"target_ids={target_ids} string-compatible candidates={len(candidates)}")
        for target_id in target_ids:
            iref_named_neighbours = {mapping2[x] for x in iref_adj[target_id] if x in mapping2}
            scored = []
            for name in candidates:
                bg = bg_adj.get(name, set()) & mapped_names
                wi = wi_adj.get(name, set()) & mapped_names
                bg_overlap = len(bg & iref_named_neighbours)
                wi_overlap = len(wi & iref_named_neighbours)
                scored.append(
                    (
                        bg_overlap + wi_overlap,
                        bg_overlap,
                        len(bg),
                        wi_overlap,
                        len(wi),
                        name,
                    )
                )
            for row in sorted(scored, reverse=True)[:15]:
                print("  score=%d BG=%d/%d WI=%d/%d %s" % row)


if __name__ == "__main__":
    main()
