#!/usr/bin/env python3
"""Convert ten ordinary-Roman runs into protein-level label-2 consensus sets."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_nodes(value: str) -> set[int]:
    return {int(token) for token in (value or "").split()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--source", default="reference ordinary-Roman VNS-G-FI")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    mapping_path = root / "data_reconstruction" / "output" / "5Layers" / "node_mapping.tsv"
    with args.input_csv.open(newline="", encoding="utf-8-sig") as stream:
        runs = list(csv.DictReader(stream))
    if len(runs) != 10:
        raise ValueError(f"expected exactly 10 runs, found {len(runs)}")
    required = {"label2_nodes", "label1_nodes"}
    if not required.issubset(runs[0]):
        raise ValueError(f"missing columns: {sorted(required - set(runs[0]))}")

    label2_sets: list[set[int]] = []
    for index, run in enumerate(runs):
        label2 = parse_nodes(run["label2_nodes"])
        label1 = parse_nodes(run["label1_nodes"])
        if label2 & label1:
            raise ValueError(f"run {index} has nodes labeled both 1 and 2")
        if not label2 | label1 <= set(range(1, 4224)):
            raise ValueError(f"run {index} contains an invalid node ID")
        objective_key = next((key for key in ("objective", "vns_obj", "greedy_obj") if run.get(key)), None)
        if objective_key and int(float(run[objective_key])) != 2 * len(label2) + len(label1):
            raise ValueError(f"run {index} objective does not match its labels")
        feasible_key = next((key for key in ("feasible", "vns_feasible", "greedy_feasible") if key in run), None)
        if feasible_key and str(run[feasible_key]).strip().lower() not in {"1", "true", "yes"}:
            raise ValueError(f"run {index} is not recorded as feasible")
        label2_sets.append(label2)

    with mapping_path.open(newline="", encoding="utf-8-sig") as stream:
        mapping = list(csv.DictReader(stream, delimiter="\t"))
    if len(mapping) != 4223:
        raise ValueError(f"expected 4223 mapped proteins, found {len(mapping)}")

    records = []
    strict_names: list[str] = []
    robust_names: list[str] = []
    for mapped in mapping:
        node = int(mapped["internal_id"])
        count = sum(node in nodes for nodes in label2_sets)
        strict = count == 10
        robust = count >= 8
        if strict:
            strict_names.append(mapped["systematic_name"])
        if robust:
            robust_names.append(mapped["systematic_name"])
        records.append({
            "internal_id": node,
            "systematic_name": mapped["systematic_name"],
            "R2_count": count,
            "R2_strict": int(strict),
            "R2_8of10": int(robust),
        })

    output = biology / "ppi5_roman_k5_protein_consensus.tsv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    (biology / "R2_strict.txt").write_text("\n".join(strict_names) + "\n", encoding="utf-8")
    (biology / "R2_8of10.txt").write_text("\n".join(robust_names) + "\n", encoding="utf-8")
    manifest = {
        "source": args.source,
        "input_csv": str(args.input_csv.resolve()),
        "run_count": 10,
        "background_size": 4223,
        "sets": {"R2_strict": len(strict_names), "R2_8of10": len(robust_names)},
        "objectives": [
            float(run[next(key for key in ("objective", "vns_obj", "greedy_obj") if run.get(key))])
            for run in runs
        ],
    }
    (biology / "ppi5_roman_k5_consensus_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
