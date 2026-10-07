#!/usr/bin/env python3
"""Create seed-consensus protein sets for five-layer k=M biological analyses."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    result_file = root / "Results" / "formal" / "ppi_5layer_k5_vns_10seeds.csv"
    mapping_file = (root / "data_reconstruction" / "output" / "5Layers" /
                    "node_mapping.tsv")
    output = root / "Results" / "biology"
    output.mkdir(parents=True, exist_ok=True)

    with result_file.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 10:
        raise ValueError(f"expected 10 solution rows, found {len(rows)}")
    vectors = [[int(value) for value in row["labels"].split()] for row in rows]
    if {len(vector) for vector in vectors} != {4223}:
        raise ValueError("unexpected label-vector length")

    with mapping_file.open(newline="", encoding="utf-8-sig") as stream:
        mapping = list(csv.DictReader(stream, delimiter="\t"))
    if len(mapping) != 4223:
        raise ValueError(f"expected 4223 mapped proteins, found {len(mapping)}")
    if [int(row["internal_id"]) for row in mapping] != list(range(1, 4224)):
        raise ValueError("mapping IDs are not the expected consecutive 1-based IDs")

    records = []
    set_names = {
        "C3_strict": [], "C2_strict": [], "C23_strict": [],
        "C3_8of10": [], "C2_8of10": [], "C23_8of10": [],
    }
    for index, mapped in enumerate(mapping):
        counts = Counter(vector[index] for vector in vectors)
        count_23 = counts[2] + counts[3]
        flags = {
            "C3_strict": counts[3] == 10,
            "C2_strict": counts[2] == 10,
            "C23_strict": count_23 == 10,
            "C3_8of10": counts[3] >= 8,
            "C2_8of10": counts[2] >= 8,
            "C23_8of10": count_23 >= 8,
        }
        for name, selected in flags.items():
            if selected:
                set_names[name].append(mapped["systematic_name"])
        records.append({
            "internal_id": mapped["internal_id"],
            "systematic_name": mapped["systematic_name"],
            "label0_count": counts[0],
            "label1_count": counts[1],
            "label2_count": counts[2],
            "label3_count": counts[3],
            "label2_or_3_count": count_23,
            **{name: int(selected) for name, selected in flags.items()},
        })

    table = output / "ppi5_k5_protein_consensus.tsv"
    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    for name, proteins in set_names.items():
        (output / f"{name}.txt").write_text("\n".join(proteins) + "\n", encoding="utf-8")

    manifest = {
        "background": "five-layer common PPI vertex set",
        "background_size": len(mapping),
        "solution_count": len(vectors),
        "sets": {name: len(proteins) for name, proteins in set_names.items()},
        "mapping_coverage": len(mapping),
    }
    (output / "ppi5_k5_consensus_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
