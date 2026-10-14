#!/usr/bin/env python3
"""Reproducible batch runner for the k-MLDRDP manuscript experiments."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from pathlib import Path

from exact_ilp import load_instance


CPP_METHODS = {
    "greedy", "ga-r-fi", "ga-r-bi", "ga-g-fi", "ga-g-bi",
    "vns-r-fi", "vns-r-bi", "vns-g-fi", "vns-g-bi",
}
HYBRID_METHODS = {"ga-r-hybrid", "ga-g-hybrid"}
EXACT_METHODS = {"dr-new1", "dr-new2"}
ALL_METHODS = CPP_METHODS | EXACT_METHODS | HYBRID_METHODS
FIELD_RE = re.compile(r"^([A-Z_]+):\s*(.*?)\s*$", re.MULTILINE)
SOLUTION_RE = re.compile(r"^SOLUTION:\r?\n[ \t]*(.+)$", re.MULTILINE)
FIELDS = [
    "instance", "n", "layers", "k", "method", "seed", "status",
    "feasible", "objective", "valid_layers", "elapsed_seconds",
    "time_to_best_seconds", "optimal", "relative_gap", "labels", "error",
]


def discover_instances(root: Path) -> list[Path]:
    found = {path.parent for path in root.rglob("layer_1.txt")}
    if (root / "layer_1.txt").is_file():
        found.add(root)
    return sorted(found, key=lambda path: str(path).lower())


def parse_csv_ints(value: str) -> list[int]:
    values = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise ValueError("empty integer list")
    return values


def cpp_command(root: Path, method: str, instance: Path, k: int,
                seed: int, limit: float) -> list[str]:
    if method == "greedy":
        return [str(root / "bin" / "kmdrdp_greedy.exe"), str(instance),
                str(k), str(seed), "10"]
    if method.startswith("ga-"):
        initialization = "greedy" if "-g-" in method else "random"
        local_search = "2" if method.endswith("-bi") else "1"
        return [
            str(root / "bin" / "kmdrdp_ga.exe"), str(instance), str(k),
            "50", "100", "2", "10", "TOURNAMENT", "1", str(seed),
            local_search, str(limit), "--init", initialization,
        ]
    initialization = "greedy" if "-g-" in method else "random"
    local_search = "2" if method.endswith("-bi") else "1"
    return [
        str(root / "bin" / "kmdrdp_vns.exe"), "vns", str(instance), str(k),
        str(seed), "2", "12", "0.4", local_search, "10000",
        "5000", str(limit), "--init", initialization,
    ]


def run_one(root: Path, python: str, method: str, instance: Path,
            k: int, seed: int, limit: float) -> dict:
    if method in CPP_METHODS:
        command = cpp_command(root, method, instance, k, seed, limit)
    elif method in EXACT_METHODS:
        model = "direct" if method == "dr-new1" else "cumulative"
        command = [python, str(root / "python" / "exact_ilp.py"),
                   str(instance), str(k), "--model", model,
                   "--time-limit", str(limit)]
    else:
        initialization = "greedy" if method == "ga-g-hybrid" else "random"
        command = [
            python, str(root / "python" / "ga_ilp_hybrid.py"),
            str(instance), str(k), "--seed", str(seed), "--time-limit",
            str(limit), "--ga-time", str(max(0.1, limit * 0.30)),
            "--subproblem-time", str(max(0.1, min(10.0, limit * 0.10))),
            "--generations", "100", "--init", initialization,
        ]

    started = time.monotonic()
    completed = subprocess.run(
        command, capture_output=True, text=True, encoding="utf-8",
        errors="replace",
    )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"exit {completed.returncode}: "
            f"{(completed.stderr or completed.stdout).strip()}"
        )

    if method in CPP_METHODS:
        fields = dict(FIELD_RE.findall(completed.stdout))
        solution = SOLUTION_RE.search(completed.stdout)
        return {
            "status": "ok",
            "feasible": fields.get("FEASIBLE", "NO") == "YES",
            "objective": int(fields["OBJECTIVE"]),
            "valid_layers": int(fields["VALID_LAYERS"]),
            "elapsed_seconds": float(fields.get("ELAPSED_SECONDS", elapsed)),
            "time_to_best_seconds": float(fields.get("TIME_TO_BEST_SECONDS", "nan")),
            "optimal": "",
            "relative_gap": "",
            "labels": "" if solution is None else solution.group(1).strip(),
        }

    data = json.loads(completed.stdout)
    return {
        "status": "ok",
        "feasible": data.get("feasible", True),
        "objective": data["objective"],
        "valid_layers": len(data["valid_layers"]),
        "elapsed_seconds": data.get("elapsed_seconds", elapsed),
        "time_to_best_seconds": "",
        "optimal": data.get("optimal", ""),
        "relative_gap": data.get("relative_gap", ""),
        "labels": " ".join(map(str, data["labels"])),
    }


def load_completed(path: Path) -> set[tuple[str, str, str, str]]:
    if not path.exists():
        return set()
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return {
            (row["instance"], row["k"], row["method"], row["seed"])
            for row in csv.DictReader(stream) if row.get("status") == "ok"
        }


def append_row(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        if not exists:
            writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in FIELDS})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("instances_root", type=Path)
    parser.add_argument("--methods", default="ga-g-fi,vns-g-fi,ga-g-hybrid")
    parser.add_argument("--k", default="1")
    parser.add_argument("--seeds", default="13579")
    parser.add_argument("--time-limit", type=float, default=60.0)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--output", type=Path, default=Path("Results/batch.csv"))
    parser.add_argument("--max-instances", type=int)
    parser.add_argument("--include-regex")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    methods = [item.strip().lower() for item in args.methods.split(",") if item.strip()]
    unknown = set(methods) - ALL_METHODS
    if unknown:
        parser.error("unknown methods: " + ", ".join(sorted(unknown)))
    protocol_k = args.k.strip().lower() == "protocol"
    ks = [] if protocol_k else parse_csv_ints(args.k)
    seeds = parse_csv_ints(args.seeds)
    root = Path(__file__).resolve().parents[1]
    instances = discover_instances(args.instances_root)
    if args.include_regex:
        pattern = re.compile(args.include_regex, re.IGNORECASE)
        instances = [path for path in instances if pattern.search(str(path))]
    if args.max_instances is not None:
        instances = instances[:args.max_instances]
    completed = load_completed(args.output)

    for instance in instances:
        graph = load_instance(instance)
        instance_ks = ({2: [1, 2], 5: [1, 3, 5], 10: [1, 5, 10]}.get(
            graph.m, list(range(1, graph.m + 1))
        ) if protocol_k else ks)
        for k in instance_ks:
            if not 1 <= k <= graph.m:
                continue
            for method in methods:
                method_seeds = seeds[:1] if method in EXACT_METHODS else seeds
                for seed in method_seeds:
                    key = (str(instance.resolve()), str(k), method, str(seed))
                    if key in completed:
                        continue
                    print(f"RUN {instance.name} k={k} {method} seed={seed}", flush=True)
                    if args.dry_run:
                        continue
                    base = {
                        "instance": str(instance.resolve()), "n": graph.n,
                        "layers": graph.m, "k": k, "method": method,
                        "seed": seed,
                    }
                    try:
                        result = run_one(
                            root, args.python, method, instance, k, seed,
                            args.time_limit,
                        )
                        row = {**base, **result, "error": ""}
                    except Exception as error:  # preserve failures in the ledger
                        row = {**base, "status": "error", "error": str(error)}
                    append_row(args.output, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
