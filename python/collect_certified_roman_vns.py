#!/usr/bin/env python3
"""Collect, independently validate, and tabulate certified Roman VNS runs."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "Results" / "formal" / "ordinary_roman_vns_certified"
INSTANCE = ROOT / "Instances" / "PPI" / "PPI-5Layers"
OUTPUT = ROOT / "Results" / "formal" / "ppi_5layer_k5_roman_vns_certified_10seeds.csv"
SUMMARY = ROOT / "Results" / "formal" / "ppi_5layer_k5_roman_vns_certified_summary.json"
SEEDS = list(range(20261006, 20261016))


def field(text: str, name: str) -> str:
    match = re.search(rf"^{re.escape(name)}:\s*(.*?)\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing {name}")
    return match.group(1)


def load_layers() -> tuple[int, list[list[set[int]]]]:
    layers: list[list[set[int]]] = []
    n_expected = None
    for path in sorted(INSTANCE.glob("layer_*.txt"), key=lambda p: int(p.stem.split("_")[1])):
        with path.open(encoding="utf-8") as stream:
            n, _ = map(int, stream.readline().split())
            adjacency = [set() for _ in range(n)]
            for line in stream:
                if not line.strip():
                    continue
                u, v = map(int, line.split())
                u -= 1
                v -= 1
                adjacency[u].add(v)
                adjacency[v].add(u)
        if n_expected is not None and n != n_expected:
            raise ValueError("layer vertex counts differ")
        n_expected = n
        layers.append(adjacency)
    if len(layers) != 5 or n_expected is None:
        raise ValueError(f"expected five layers, found {len(layers)}")
    return n_expected, layers


def valid_layer(labels: list[int], adjacency: list[set[int]]) -> bool:
    return all(label != 0 or any(labels[j] == 2 for j in adjacency[i])
               for i, label in enumerate(labels))


def main() -> int:
    n, layers = load_layers()
    records: list[dict[str, object]] = []
    errors: list[str] = []
    for seed in SEEDS:
        path = RUNS / f"seed_{seed}" / "stdout.log"
        if not path.exists():
            errors.append(f"seed {seed}: missing stdout.log")
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        try:
            solution = re.search(r"^SOLUTION:\s*\n\s*(.*?)\s*$", text, re.MULTILINE)
            if not solution:
                raise ValueError("missing SOLUTION")
            labels = [int(value) for value in solution.group(1).split()]
            objective = int(field(text, "OBJECTIVE"))
            if len(labels) != n:
                raise ValueError(f"label count {len(labels)} != {n}")
            if any(value not in (0, 1, 2) for value in labels):
                raise ValueError("label outside {0,1,2}")
            if sum(labels) != objective:
                raise ValueError(f"label sum {sum(labels)} != objective {objective}")
            valid = sum(valid_layer(labels, adjacency) for adjacency in layers)
            if valid != 5:
                raise ValueError(f"only {valid}/5 layers independently validate")
            if field(text, "FEASIBLE") != "YES":
                raise ValueError("solver did not report feasible")
            records.append({
                "seed": seed,
                "objective": objective,
                "valid_layers": valid,
                "feasible": True,
                "iterations": int(field(text, "ITERATIONS")),
                "elapsed_seconds": float(field(text, "ELAPSED_SECONDS")),
                "time_to_best_seconds": float(field(text, "TIME_TO_BEST_SECONDS")),
                "label2_nodes": " ".join(str(i + 1) for i, value in enumerate(labels) if value == 2),
                "label1_nodes": " ".join(str(i + 1) for i, value in enumerate(labels) if value == 1),
                "labels": " ".join(map(str, labels)),
            })
        except ValueError as exc:
            errors.append(f"seed {seed}: {exc}")

    if records:
        with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    summary = {
        "complete": len(records) == len(SEEDS) and not errors,
        "expected_runs": len(SEEDS),
        "validated_runs": len(records),
        "validation_errors": errors,
        "seeds": SEEDS,
        "objectives": [record["objective"] for record in records],
        "output": str(OUTPUT.relative_to(ROOT)).replace("\\", "/"),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if not summary["complete"]:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
