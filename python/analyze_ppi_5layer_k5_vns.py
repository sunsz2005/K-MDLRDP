#!/usr/bin/env python3
"""Validate and summarize the ten-seed five-layer PPI k=5 VNS experiment."""

from __future__ import annotations

import csv
import itertools
import json
import statistics
from collections import Counter
from pathlib import Path

from analyze_ppi_gate import FIELDS, validate, write_csv


SEEDS = (13579, 24680, 97531, 86420, 112358,
         271828, 314159, 161803, 424242, 987654)


def read(path: Path, seed_filter: int | None = None) -> tuple[list[dict], list[str]]:
    rows, errors = [], []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for line_number, row in enumerate(csv.DictReader(stream), 2):
            try:
                for field in ("n", "layers", "k", "seed", "objective", "valid_layers"):
                    row[field] = int(row[field])
                row["elapsed_seconds"] = float(row["elapsed_seconds"])
                row["time_to_best_seconds"] = float(row["time_to_best_seconds"])
            except (KeyError, TypeError, ValueError) as error:
                errors.append(f"{path.name}:{line_number}: {error}")
                continue
            if seed_filter is not None and row["seed"] != seed_filter:
                continue
            row["instance"] = str(Path(row["instance"]).resolve())
            rows.append(row)
    return rows, errors


def jaccard(left: set[int], right: set[int]) -> float:
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def stability(label_sets: list[set[int]], n: int) -> dict:
    pairwise = [jaccard(a, b) for a, b in itertools.combinations(label_sets, 2)]
    frequencies = Counter(vertex for values in label_sets for vertex in values)
    return {
        "mean_set_size": statistics.mean(map(len, label_sets)),
        "pairwise_jaccard_mean": statistics.mean(pairwise),
        "pairwise_jaccard_min": min(pairwise),
        "pairwise_jaccard_max": max(pairwise),
        "vertices_in_all_10": sum(count == 10 for count in frequencies.values()),
        "vertices_in_at_least_8": sum(count >= 8 for count in frequencies.values()),
        "vertices_in_at_least_5": sum(count >= 5 for count in frequencies.values()),
        "union_size": len(frequencies),
        "background_n": n,
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    results = root / "Results" / "formal"
    gate_all, errors = read(results / "ppi_gate_vns_g_fi.csv", seed_filter=13579)
    gate = [row for row in gate_all if row["layers"] == 5 and row["k"] == 5]
    extra, extra_errors = read(results / "ppi_5layer_k5_vns_extra.csv")
    errors.extend(extra_errors)
    rows = gate + extra

    indexed = {}
    for row in rows:
        key = (Path(row["instance"]).name, row["k"], row["method"], row["seed"])
        if key in indexed:
            errors.append(f"duplicate key {key}")
        indexed[key] = row
    expected = {("PPI-5Layers", 5, "vns-g-fi", seed) for seed in SEEDS}
    missing = sorted(expected - set(indexed))
    unexpected = sorted(set(indexed) - expected)
    validation_errors = validate(indexed)
    complete = not (errors or missing or unexpected or validation_errors)

    report = {
        "complete": complete,
        "rows": len(indexed),
        "missing_count": len(missing),
        "unexpected_count": len(unexpected),
        "structural_error_count": len(errors),
        "solution_validation_error_count": len(validation_errors),
        "errors": (errors + validation_errors)[:25],
    }
    if complete:
        ordered = [indexed[("PPI-5Layers", 5, "vns-g-fi", seed)] for seed in SEEDS]
        objectives = [row["objective"] for row in ordered]
        runtimes = [row["elapsed_seconds"] for row in ordered]
        best_times = [row["time_to_best_seconds"] for row in ordered]
        labels = [[int(value) for value in row["labels"].split()] for row in ordered]
        label_counts = [Counter(vector) for vector in labels]
        greedy, hybrid = 2759, 2721
        report["objective_summary"] = {
            "mean": statistics.mean(objectives),
            "sample_std": statistics.stdev(objectives),
            "min": min(objectives),
            "max": max(objectives),
            "median": statistics.median(objectives),
            "wins_ties_losses_vs_greedy": [
                sum(value < greedy for value in objectives),
                sum(value == greedy for value in objectives),
                sum(value > greedy for value in objectives),
            ],
            "mean_improvement_vs_greedy_percent":
                100.0 * (greedy - statistics.mean(objectives)) / greedy,
            "wins_ties_losses_vs_gate_hybrid": [
                sum(value < hybrid for value in objectives),
                sum(value == hybrid for value in objectives),
                sum(value > hybrid for value in objectives),
            ],
            "mean_improvement_vs_gate_hybrid_percent":
                100.0 * (hybrid - statistics.mean(objectives)) / hybrid,
            "mean_runtime_seconds": statistics.mean(runtimes),
            "mean_time_to_best_seconds": statistics.mean(best_times),
        }
        report["mean_label_counts"] = {
            str(value): statistics.mean(counts[value] for counts in label_counts)
            for value in range(4)
        }
        report["stability"] = {
            "label_3": stability([
                {i for i, value in enumerate(vector) if value == 3}
                for vector in labels
            ], len(labels[0])),
            "label_2": stability([
                {i for i, value in enumerate(vector) if value == 2}
                for vector in labels
            ], len(labels[0])),
            "support_union_2_or_3": stability([
                {i for i, value in enumerate(vector) if value >= 2}
                for vector in labels
            ], len(labels[0])),
        }
        merged = [{field: row.get(field, "") for field in FIELDS} for row in ordered]
        write_csv(results / "ppi_5layer_k5_vns_10seeds.csv", merged, FIELDS)

    (results / "ppi_5layer_k5_vns_summary.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
