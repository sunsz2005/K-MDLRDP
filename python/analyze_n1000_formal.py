#!/usr/bin/env python3
"""Merge, validate, and analyze the n=1000 ten-seed formal experiment."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import friedmanchisquare, wilcoxon

from exact_ilp import layer_violations, load_instance


SEEDS = (13579, 24680, 97531, 86420, 112358,
         271828, 314159, 161803, 424242, 987654)
METHOD_FILES = {
    "greedy": ("n1000_gate_greedy.csv",),
    "vns-g-fi": (
        "n1000_gate_vns_g_fi.csv",
        "n1000_formal_vns_g_fi_extra.csv",
    ),
    "ga-g-hybrid": (
        "n1000_gate_ga_g_hybrid.csv",
        "n1000_formal_ga_g_hybrid_a.csv",
        "n1000_formal_ga_g_hybrid_b.csv",
    ),
}
FIELDS = [
    "instance", "n", "layers", "k", "method", "seed", "status",
    "feasible", "objective", "valid_layers", "elapsed_seconds",
    "time_to_best_seconds", "optimal", "relative_gap", "labels", "error",
]


def read_rows(results_dir: Path) -> tuple[list[dict], list[str]]:
    rows: list[dict] = []
    errors: list[str] = []
    for method, names in METHOD_FILES.items():
        for name in names:
            path = results_dir / name
            if not path.is_file():
                errors.append(f"missing file: {path}")
                continue
            with path.open(newline="", encoding="utf-8-sig") as stream:
                for line_number, row in enumerate(csv.DictReader(stream), 2):
                    if row.get("method") != method:
                        errors.append(
                            f"{name}:{line_number}: expected method {method}, "
                            f"got {row.get('method')!r}"
                        )
                    row["_source"] = name
                    row["_line"] = line_number
                    rows.append(row)
    return rows, errors


def normalize(rows: list[dict], errors: list[str]) -> dict[tuple, dict]:
    indexed: dict[tuple, dict] = {}
    for row in rows:
        source = f"{row['_source']}:{row['_line']}"
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
        row["instance"] = str(Path(row["instance"]).resolve())
        key = (row["instance"], row["k"], row["method"], row["seed"])
        if key in indexed:
            errors.append(f"{source}: duplicate key {key}")
        else:
            indexed[key] = row
    return indexed


def audit_design(indexed: dict[tuple, dict], errors: list[str]) -> list[tuple[str, int]]:
    greedy_configs = sorted({(key[0], key[1]) for key in indexed if key[2] == "greedy"})
    expected = {
        (instance, k, method, seed)
        for instance, k in greedy_configs
        for method in METHOD_FILES
        for seed in SEEDS
    }
    actual = set(indexed)
    for key in sorted(expected - actual):
        errors.append(f"missing run: {key}")
    for key in sorted(actual - expected):
        errors.append(f"unexpected run: {key}")
    if len(greedy_configs) != 24:
        errors.append(f"expected 24 configurations, found {len(greedy_configs)}")
    return greedy_configs


def validate_solutions(indexed: dict[tuple, dict], errors: list[str]) -> None:
    graphs = {}
    for key, row in indexed.items():
        source = f"{row['_source']}:{row['_line']}"
        if row.get("status") != "ok":
            errors.append(f"{source}: status is {row.get('status')!r}")
            continue
        if str(row.get("feasible", "")).strip().lower() not in {"true", "1", "yes"}:
            errors.append(f"{source}: feasible flag is not true")
        labels_text = row.get("labels", "")
        try:
            labels = [int(value) for value in labels_text.split()]
        except ValueError as error:
            errors.append(f"{source}: invalid label: {error}")
            continue
        instance = row["instance"]
        if instance not in graphs:
            graphs[instance] = load_instance(Path(instance))
        graph = graphs[instance]
        if row["n"] != graph.n or row["layers"] != graph.m:
            errors.append(f"{source}: graph dimensions disagree with instance")
            continue
        if len(labels) != graph.n:
            errors.append(f"{source}: expected {graph.n} labels, found {len(labels)}")
            continue
        if any(label not in (0, 1, 2, 3) for label in labels):
            errors.append(f"{source}: labels outside {{0,1,2,3}}")
            continue
        if sum(labels) != row["objective"]:
            errors.append(
                f"{source}: label sum {sum(labels)} != objective {row['objective']}"
            )
        valid_layers = sum(value == 0 for value in layer_violations(graph, labels))
        if valid_layers != row["valid_layers"]:
            errors.append(
                f"{source}: recomputed valid layers {valid_layers} != "
                f"recorded {row['valid_layers']}"
            )
        if valid_layers < row["k"]:
            errors.append(f"{source}: only {valid_layers} valid layers for k={row['k']}")


def holm_adjust(p_values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(p_values, key=p_values.get)
    adjusted = {}
    running = 0.0
    count = len(ordered)
    for rank, name in enumerate(ordered):
        value = min(1.0, (count - rank) * p_values[name])
        running = max(running, value)
        adjusted[name] = running
    return adjusted


def paired_stats(a: list[dict], b: list[dict]) -> dict:
    av = np.asarray([row["objective"] for row in a], dtype=float)
    bv = np.asarray([row["objective"] for row in b], dtype=float)
    diff = bv - av
    nonzero = diff[diff != 0]
    p_value = 1.0 if len(nonzero) == 0 else float(
        wilcoxon(
            av,
            bv,
            alternative="two-sided",
            zero_method="wilcox",
            method="approx" if np.any(diff == 0) else "auto",
        ).pvalue
    )
    return {
        "a": a[0]["method"],
        "b": b[0]["method"],
        "count": len(a),
        "wins_ties_losses": [int(np.sum(diff > 0)), int(np.sum(diff == 0)),
                              int(np.sum(diff < 0))],
        "mean_objective_a": float(np.mean(av)),
        "mean_objective_b": float(np.mean(bv)),
        "mean_objective_delta_b_minus_a": float(np.mean(diff)),
        "mean_relative_improvement_percent": float(np.mean(100.0 * diff / bv)),
        "mean_runtime_a_seconds": float(np.mean([row["elapsed_seconds"] for row in a])),
        "mean_runtime_b_seconds": float(np.mean([row["elapsed_seconds"] for row in b])),
        "mean_paired_runtime_ratio_a_over_b": float(np.mean([
            left["elapsed_seconds"] / right["elapsed_seconds"]
            for left, right in zip(a, b)
        ])),
        "wilcoxon_two_sided_p": p_value,
    }


def available_pair_stats(indexed: dict[tuple, dict], method_a: str,
                         method_b: str) -> dict:
    """Summarize only currently available matched configuration/seed runs."""
    a_rows = {
        (instance, k, seed): row
        for (instance, k, method, seed), row in indexed.items()
        if method == method_a
    }
    b_rows = {
        (instance, k, seed): row
        for (instance, k, method, seed), row in indexed.items()
        if method == method_b
    }
    matched = sorted(set(a_rows) & set(b_rows))
    if not matched:
        return {"a": method_a, "b": method_b, "count": 0}
    return paired_stats(
        [a_rows[key] for key in matched],
        [b_rows[key] for key in matched],
    )


def analyze(indexed: dict[tuple, dict], configs: list[tuple[str, int]]) -> tuple[dict, list[dict]]:
    methods = list(METHOD_FILES)
    by_run = {}
    for instance, k in configs:
        for seed in SEEDS:
            by_run[(instance, k, seed)] = {
                method: indexed[(instance, k, method, seed)] for method in methods
            }

    run_pairs = {}
    config_pairs = {}
    config_rows = []
    config_method_rows: dict[str, list[dict]] = defaultdict(list)
    for instance, k in configs:
        rows_by_method = {
            method: [indexed[(instance, k, method, seed)] for seed in SEEDS]
            for method in methods
        }
        output_row = {
            "instance": instance,
            "layers": rows_by_method["greedy"][0]["layers"],
            "k": k,
        }
        for method, rows in rows_by_method.items():
            objectives = [row["objective"] for row in rows]
            runtimes = [row["elapsed_seconds"] for row in rows]
            summary_row = {
                "method": method,
                "objective": float(np.mean(objectives)),
                "elapsed_seconds": float(np.mean(runtimes)),
            }
            config_method_rows[method].append(summary_row)
            output_row[f"{method}_objective_mean"] = summary_row["objective"]
            output_row[f"{method}_objective_sd"] = float(np.std(objectives, ddof=1))
            output_row[f"{method}_runtime_mean"] = summary_row["elapsed_seconds"]
        config_rows.append(output_row)

    pair_names = [
        ("vns-g-fi", "greedy"),
        ("ga-g-hybrid", "greedy"),
        ("vns-g-fi", "ga-g-hybrid"),
    ]
    raw_run_p = {}
    raw_config_p = {}
    for method_a, method_b in pair_names:
        name = f"{method_a}_vs_{method_b}"
        run_a, run_b = [], []
        for key in sorted(by_run):
            run_a.append(by_run[key][method_a])
            run_b.append(by_run[key][method_b])
        run_pairs[name] = paired_stats(run_a, run_b)
        raw_run_p[name] = run_pairs[name]["wilcoxon_two_sided_p"]
        config_pairs[name] = paired_stats(
            config_method_rows[method_a], config_method_rows[method_b]
        )
        raw_config_p[name] = config_pairs[name]["wilcoxon_two_sided_p"]

    for name, adjusted in holm_adjust(raw_run_p).items():
        run_pairs[name]["holm_adjusted_p"] = adjusted
    for name, adjusted in holm_adjust(raw_config_p).items():
        config_pairs[name]["holm_adjusted_p"] = adjusted

    arrays = [
        [row["objective"] for row in config_method_rows[method]] for method in methods
    ]
    friedman = friedmanchisquare(*arrays)
    best_counts = Counter()
    unique_best_counts = Counter()
    for config_index in range(len(configs)):
        values = {
            method: config_method_rows[method][config_index]["objective"]
            for method in methods
        }
        best = min(values.values())
        winners = [method for method, value in values.items() if value == best]
        best_counts.update(winners)
        if len(winners) == 1:
            unique_best_counts.update(winners)

    report = {
        "design": {
            "n": 1000,
            "configurations": len(configs),
            "seeds": list(SEEDS),
            "methods": methods,
            "runs": len(indexed),
            "runs_by_method": Counter(key[2] for key in indexed),
        },
        "run_level_pairs": run_pairs,
        "configuration_mean_pairs": config_pairs,
        "friedman_configuration_means": {
            "statistic": float(friedman.statistic),
            "p_value": float(friedman.pvalue),
        },
        "configuration_mean_best_counts": {
            method: {
                "best_or_tied": best_counts[method],
                "unique_best": unique_best_counts[method],
            }
            for method in methods
        },
    }
    return report, config_rows


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0])
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

    rows, errors = read_rows(args.results_dir)
    indexed = normalize(rows, errors)
    configs = audit_design(indexed, errors)
    solution_errors: list[str] = []
    validate_solutions(indexed, solution_errors)
    complete = not errors and not solution_errors and len(indexed) == 720
    if not complete:
        counts = Counter(key[2] for key in indexed)
        available_pairs = {
            f"{method_a}_vs_{method_b}": available_pair_stats(
                indexed, method_a, method_b
            )
            for method_a, method_b in (
                ("vns-g-fi", "greedy"),
                ("ga-g-hybrid", "greedy"),
                ("vns-g-fi", "ga-g-hybrid"),
            )
        }
        print(json.dumps({
            "complete": False,
            "rows": len(indexed),
            "rows_by_method": counts,
            "error_count": len(errors),
            "errors": errors[:25],
            "checked_available_solutions": len(indexed),
            "solution_validation_error_count": len(solution_errors),
            "solution_validation_errors": solution_errors[:25],
            "available_matched_pairs": available_pairs,
        }, indent=2, ensure_ascii=False, default=dict))
        return 0 if args.allow_incomplete else 2

    merged_rows = []
    for key in sorted(indexed):
        merged_rows.append({field: indexed[key].get(field, "") for field in FIELDS})
    report, config_rows = analyze(indexed, configs)
    report["validation"] = {
        "checked_rows": len(indexed),
        "error_count": 0,
        "labels_and_constraints_recomputed": True,
    }
    write_csv(args.results_dir / "n1000_formal_10seed.csv", merged_rows, FIELDS)
    write_csv(args.results_dir / "n1000_formal_by_config.csv", config_rows)
    summary_path = args.results_dir / "n1000_formal_summary.json"
    summary_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=dict) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, ensure_ascii=False, default=dict))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
