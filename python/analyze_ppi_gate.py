#!/usr/bin/env python3
"""Validate and summarize the one-seed real-PPI k-MLDRDP gate."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from exact_ilp import layer_violations, load_instance


METHOD_FILES = {
    "greedy": "ppi_gate_greedy.csv",
    "vns-g-fi": "ppi_gate_vns_g_fi.csv",
    "ga-g-hybrid": "ppi_gate_ga_g_hybrid.csv",
}
EXPECTED = {
    ("PPI-3Layers", k, method, 13579)
    for k in range(1, 4) for method in METHOD_FILES
} | {
    ("PPI-4Layers", k, method, 13579)
    for k in range(1, 5) for method in METHOD_FILES
} | {
    ("PPI-5Layers", k, method, 13579)
    for k in range(1, 6) for method in METHOD_FILES
}
FIELDS = [
    "instance", "n", "layers", "k", "method", "seed", "status",
    "feasible", "objective", "valid_layers", "elapsed_seconds",
    "time_to_best_seconds", "optimal", "relative_gap", "labels", "error",
]


def read_rows(results_dir: Path) -> tuple[dict[tuple, dict], list[str]]:
    indexed = {}
    errors = []
    for method, name in METHOD_FILES.items():
        path = results_dir / name
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8-sig") as stream:
            for line_number, row in enumerate(csv.DictReader(stream), 2):
                source = f"{name}:{line_number}"
                try:
                    row["n"] = int(row["n"])
                    row["layers"] = int(row["layers"])
                    row["k"] = int(row["k"])
                    row["seed"] = int(row["seed"])
                    row["objective"] = int(row["objective"])
                    row["valid_layers"] = int(row["valid_layers"])
                    row["elapsed_seconds"] = float(row["elapsed_seconds"])
                except (KeyError, TypeError, ValueError) as error:
                    errors.append(f"{source}: invalid numeric field: {error}")
                    continue
                if row.get("method") != method:
                    errors.append(f"{source}: expected {method}, got {row.get('method')}")
                row["instance"] = str(Path(row["instance"]).resolve())
                key = (Path(row["instance"]).name, row["k"], method, row["seed"])
                if key in indexed:
                    errors.append(f"{source}: duplicate key {key}")
                indexed[key] = row
    return indexed, errors


def validate(indexed: dict[tuple, dict]) -> list[str]:
    errors = []
    graphs = {}
    for key, row in indexed.items():
        source = f"{key}"
        if row.get("status") != "ok":
            errors.append(f"{source}: status {row.get('status')}")
            continue
        if str(row.get("feasible", "")).lower() not in {"true", "1", "yes"}:
            errors.append(f"{source}: feasible flag false")
        try:
            labels = [int(value) for value in row.get("labels", "").split()]
        except ValueError as error:
            errors.append(f"{source}: invalid labels: {error}")
            continue
        instance = row["instance"]
        if instance not in graphs:
            graphs[instance] = load_instance(instance)
        graph = graphs[instance]
        if row["n"] != graph.n or row["layers"] != graph.m:
            errors.append(f"{source}: graph dimensions mismatch")
            continue
        if len(labels) != graph.n:
            errors.append(f"{source}: {len(labels)} labels, expected {graph.n}")
            continue
        if any(value not in (0, 1, 2, 3) for value in labels):
            errors.append(f"{source}: label outside 0..3")
            continue
        if sum(labels) != row["objective"]:
            errors.append(f"{source}: label sum differs from objective")
        valid = sum(value == 0 for value in layer_violations(graph, labels))
        if valid != row["valid_layers"]:
            errors.append(f"{source}: valid-layer count {valid} != {row['valid_layers']}")
        if valid < row["k"]:
            errors.append(f"{source}: only {valid} valid layers for k={row['k']}")
    return errors


def paired(indexed: dict[tuple, dict], a: str, b: str) -> dict:
    matches = []
    for instance, k, method, seed in sorted(indexed):
        if method != a:
            continue
        other = (instance, k, b, seed)
        if other in indexed:
            matches.append((indexed[(instance, k, a, seed)], indexed[other]))
    if not matches:
        return {"a": a, "b": b, "count": 0}
    differences = [right["objective"] - left["objective"] for left, right in matches]
    return {
        "a": a,
        "b": b,
        "count": len(matches),
        "wins_ties_losses": [
            sum(value > 0 for value in differences),
            sum(value == 0 for value in differences),
            sum(value < 0 for value in differences),
        ],
        "mean_relative_improvement_percent": sum(
            100.0 * (right["objective"] - left["objective"]) / right["objective"]
            for left, right in matches
        ) / len(matches),
        "mean_runtime_a_seconds": sum(row[0]["elapsed_seconds"] for row in matches) / len(matches),
        "mean_runtime_b_seconds": sum(row[1]["elapsed_seconds"] for row in matches) / len(matches),
        "mean_runtime_ratio_a_over_b": sum(
            left["elapsed_seconds"] / right["elapsed_seconds"]
            for left, right in matches
        ) / len(matches),
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "Results" / "formal")
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()

    indexed, structural_errors = read_rows(args.results_dir)
    solution_errors = validate(indexed)
    missing = sorted(EXPECTED - set(indexed))
    unexpected = sorted(set(indexed) - EXPECTED)
    complete = not structural_errors and not solution_errors and not missing and not unexpected
    report = {
        "complete": complete,
        "rows": len(indexed),
        "rows_by_method": dict(Counter(key[2] for key in indexed)),
        "missing_count": len(missing),
        "unexpected_count": len(unexpected),
        "structural_error_count": len(structural_errors),
        "solution_validation_error_count": len(solution_errors),
        "errors": (structural_errors + solution_errors)[:25],
        "pairs": {
            "vns_vs_greedy": paired(indexed, "vns-g-fi", "greedy"),
            "hybrid_vs_greedy": paired(indexed, "ga-g-hybrid", "greedy"),
            "vns_vs_hybrid": paired(indexed, "vns-g-fi", "ga-g-hybrid"),
        },
    }
    if not complete:
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0 if args.allow_incomplete else 2

    merged = [{field: indexed[key].get(field, "") for field in FIELDS}
              for key in sorted(indexed)]
    by_config = []
    for instance in ("PPI-3Layers", "PPI-4Layers", "PPI-5Layers"):
        m = int(instance[4])
        for k in range(1, m + 1):
            row = {"instance": instance, "layers": m, "k": k}
            objectives = {}
            for method in METHOD_FILES:
                result = indexed[(instance, k, method, 13579)]
                objectives[method] = result["objective"]
                row[f"{method}_objective"] = result["objective"]
                row[f"{method}_elapsed_seconds"] = result["elapsed_seconds"]
            best = min(objectives.values())
            row["best_methods"] = ",".join(
                method for method, value in objectives.items() if value == best
            )
            by_config.append(row)
    report["by_config"] = by_config
    report["validation"] = {"checked_rows": len(indexed), "error_count": 0}
    write_csv(args.results_dir / "ppi_gate_merged.csv", merged, FIELDS)
    write_csv(args.results_dir / "ppi_gate_by_config.csv", by_config, list(by_config[0]))
    (args.results_dir / "ppi_gate_summary.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
