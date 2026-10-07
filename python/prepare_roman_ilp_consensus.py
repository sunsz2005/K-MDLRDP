#!/usr/bin/env python3
"""Create the ordinary-Roman dominant set as the intersection of two ILPs."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    formal = root / "Results" / "formal"
    biology = root / "Results" / "biology"
    solutions = [
        json.loads((formal / "ppi5_roman_ilp_direct.json").read_text(encoding="utf-8")),
        json.loads((formal / "ppi5_roman_ilp_cumulative.json").read_text(encoding="utf-8")),
    ]
    labels = [[int(value) for value in item["labels"].split()] for item in solutions]
    if any(len(vector) != 4223 for vector in labels):
        raise ValueError("ordinary-Roman ILP label vector has the wrong length")
    mapping_path = root / "data_reconstruction" / "output" / "5Layers" / "node_mapping.tsv"
    with mapping_path.open(newline="", encoding="utf-8-sig") as stream:
        mapping = list(csv.DictReader(stream, delimiter="\t"))
    records, strict_names = [], []
    for index, mapped in enumerate(mapping):
        count = sum(vector[index] == 2 for vector in labels)
        selected = count == 2
        if selected:
            strict_names.append(mapped["systematic_name"])
        records.append({
            "internal_id": mapped["internal_id"],
            "systematic_name": mapped["systematic_name"],
            "R2_count": count,
            "R2_strict": int(selected),
            "R2_8of10": int(selected)
        })
    with (biology / "ppi5_roman_k5_protein_consensus.tsv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    for name in ("R2_strict.txt", "R2_8of10.txt"):
        (biology / name).write_text("\n".join(strict_names) + "\n", encoding="utf-8")
    manifest = {
        "definition": "label 2 in both ordinary-Roman ILP formulations",
        "set_size": len(strict_names),
        "solutions": [{key: item.get(key) for key in
                       ("formulation", "objective", "status", "success", "mip_gap", "mip_dual_bound")}
                      for item in solutions],
    }
    (biology / "ppi5_roman_k5_consensus_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
