#!/usr/bin/env python3
"""Validate and summarize the 10-seed GA-G-Hybrid PPI k=1 experiment."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

from analyze_ppi_gate import FIELDS, validate, write_csv


SEEDS = (13579, 24680, 97531, 86420, 112358,
         271828, 314159, 161803, 424242, 987654)
INSTANCES = ("PPI-3Layers", "PPI-4Layers", "PPI-5Layers")


def load(path: Path, *, only_k1: bool = False) -> tuple[list[dict], list[str]]:
    rows, errors = [], []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for line_number, row in enumerate(csv.DictReader(stream), 2):
            source = f"{path.name}:{line_number}"
            try:
                for field in ("n", "layers", "k", "seed", "objective", "valid_layers"):
                    row[field] = int(row[field])
                row["elapsed_seconds"] = float(row["elapsed_seconds"])
            except (KeyError, TypeError, ValueError) as error:
                errors.append(f"{source}: invalid numeric field: {error}")
                continue
            if only_k1 and row["k"] != 1:
                continue
            row["instance"] = str(Path(row["instance"]).resolve())
            rows.append(row)
    return rows, errors


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    results = root / "Results" / "formal"
    gate, errors = load(results / "ppi_gate_ga_g_hybrid.csv", only_k1=True)
    extra, extra_errors = load(results / "ppi_formal_ga_g_hybrid_k1_extra.csv")
    errors.extend(extra_errors)
    rows = gate + extra

    indexed = {}
    for row in rows:
        key = (Path(row["instance"]).name, row["k"], row["method"], row["seed"])
        if key in indexed:
            errors.append(f"duplicate key {key}")
        indexed[key] = row

    expected = {
        (instance, 1, "ga-g-hybrid", seed)
        for instance in INSTANCES for seed in SEEDS
    }
    missing = sorted(expected - set(indexed))
    unexpected = sorted(set(indexed) - expected)
    validation_errors = validate(indexed)

    greedy_rows, greedy_errors = load(results / "ppi_gate_greedy.csv", only_k1=True)
    errors.extend(greedy_errors)
    greedy = {Path(row["instance"]).name: row["objective"] for row in greedy_rows}

    summaries = []
    if not (errors or missing or unexpected or validation_errors):
        for instance in INSTANCES:
            group = [indexed[(instance, 1, "ga-g-hybrid", seed)] for seed in SEEDS]
            objectives = [row["objective"] for row in group]
            runtimes = [row["elapsed_seconds"] for row in group]
            baseline = greedy[instance]
            summaries.append({
                "instance": instance,
                "layers": int(instance[4]),
                "runs": len(group),
                "greedy_objective": baseline,
                "hybrid_mean": statistics.mean(objectives),
                "hybrid_sample_std": statistics.stdev(objectives),
                "hybrid_min": min(objectives),
                "hybrid_max": max(objectives),
                "mean_improvement_percent": 100.0 * (baseline - statistics.mean(objectives)) / baseline,
                "worst_run_improvement_percent": 100.0 * (baseline - max(objectives)) / baseline,
                "mean_runtime_seconds": statistics.mean(runtimes),
                "min_runtime_seconds": min(runtimes),
                "max_runtime_seconds": max(runtimes),
                "wins_vs_greedy": sum(value < baseline for value in objectives),
                "ties_vs_greedy": sum(value == baseline for value in objectives),
                "losses_vs_greedy": sum(value > baseline for value in objectives),
            })

    complete = not (errors or missing or unexpected or validation_errors)
    report = {
        "complete": complete,
        "rows": len(indexed),
        "missing_count": len(missing),
        "unexpected_count": len(unexpected),
        "structural_error_count": len(errors),
        "solution_validation_error_count": len(validation_errors),
        "errors": (errors + validation_errors)[:25],
        "summary": summaries,
    }
    if complete:
        merged = [{field: indexed[key].get(field, "") for field in FIELDS}
                  for key in sorted(indexed)]
        write_csv(results / "ppi_formal_ga_g_hybrid_k1_10seeds.csv", merged, FIELDS)
        write_csv(results / "ppi_formal_ga_g_hybrid_k1_summary.csv",
                  summaries, list(summaries[0]))
    (results / "ppi_formal_ga_g_hybrid_k1_summary.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
