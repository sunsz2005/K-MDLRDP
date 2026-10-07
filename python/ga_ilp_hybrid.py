#!/usr/bin/env python3
"""Support-aware GA--ILP hybrid for k-MLDRDP.

V2 preserves a complete compiled-GA run, then applies restricted DR-New1
calls. A strict ILP improvement is injected into a short follow-up GA epoch.
Neighborhoods are rooted at current label-2/3 support vertices and grow after
repeated failures. The strict cutoff prevents the MIP returning the incumbent.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Sequence

from exact_ilp import is_feasible, layer_violations, load_instance, solve_exact


SOLUTION_RE = re.compile(r"^SOLUTION:\r?\n[ \t]*(.+)$", re.MULTILINE)


def run_ga(
    executable: Path,
    instance: Path,
    k: int,
    seed: int,
    time_limit: float,
    population: int,
    generations: int,
    initialization: str,
    initial_labels: Sequence[int] | None = None,
) -> tuple[list[int], str]:
    command = [
        str(executable), str(instance), str(k), str(population),
        str(generations), "2", "10", "TOURNAMENT", "1", str(seed), "1",
        str(max(0.1, time_limit)), "--init", initialization,
    ]
    label_path: Path | None = None
    try:
        if initial_labels is not None:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".labels", encoding="ascii", delete=False
            ) as stream:
                stream.write(" ".join(map(str, initial_labels)) + "\n")
                label_path = Path(stream.name)
            command.extend(("--initial-solution", str(label_path)))
        completed = subprocess.run(
            command, check=True, capture_output=True, text=True,
            encoding="utf-8", errors="replace",
        )
    finally:
        if label_path is not None:
            label_path.unlink(missing_ok=True)
    match = SOLUTION_RE.search(completed.stdout)
    if match is None:
        raise RuntimeError(
            "could not parse the GA solution; output was:\n" + completed.stdout
        )
    labels = [int(token) for token in match.group(1).split()]
    return labels, completed.stdout


def support_aware_free_set(
    graph, labels: Sequence[int], k: int, target_size: int, rng: random.Random
) -> set[int]:
    """Build a multi-layer neighborhood around expensive support vertices."""
    target_size = min(graph.n, max(1, target_size))
    violations = layer_violations(graph, labels)
    valid_layers = [index for index, value in enumerate(violations) if value == 0]
    if len(valid_layers) > k:
        rng.shuffle(valid_layers)
        valid_layers = valid_layers[:k]
    elif not valid_layers:
        valid_layers = sorted(range(graph.m), key=lambda index: violations[index])[:k]

    active = [vertex for vertex, label in enumerate(labels) if label >= 2]
    if active:
        seed = rng.choices(active, weights=[labels[vertex] for vertex in active], k=1)[0]
    else:
        seed = rng.randrange(graph.n)
    chosen = {seed}

    # Free several incumbent support vertices so that 3->(2,2), support
    # relocation, and other multi-vertex exchanges are representable.
    other_active = [vertex for vertex in active if vertex != seed]
    rng.shuffle(other_active)
    minimum_active = min(len(active), max(1, math.ceil(target_size * 0.30)))
    chosen.update(other_active[:max(0, minimum_active - 1)])

    frontier = list(chosen)
    rng.shuffle(frontier)
    while frontier and len(chosen) < target_size:
        vertex = frontier.pop()
        candidates: list[int] = []
        for layer in valid_layers:
            candidates.extend(graph.layers[layer][vertex])
        candidates = list(dict.fromkeys(candidates))
        rng.shuffle(candidates)
        candidates.sort(key=lambda item: labels[item] < 2)
        for neighbor in candidates:
            if neighbor not in chosen:
                chosen.add(neighbor)
                frontier.append(neighbor)
                if len(chosen) == target_size:
                    break

    if len(chosen) < target_size:
        remaining = [vertex for vertex in range(graph.n) if vertex not in chosen]
        chosen.update(rng.sample(remaining, target_size - len(chosen)))
    return chosen


def restricted_improvement(
    graph,
    k: int,
    incumbent: list[int],
    target_size: int,
    rng: random.Random,
    time_limit: float,
) -> tuple[list[int], dict]:
    free_vertices = support_aware_free_set(
        graph, incumbent, k, target_size, rng
    )
    fixed = [None if vertex in free_vertices else incumbent[vertex]
             for vertex in range(graph.n)]
    before = sum(incumbent)
    try:
        result = solve_exact(
            graph,
            k,
            "direct",
            time_limit,
            fixed_labels=fixed,
            objective_upper_bound=before - 1,
        )
    except RuntimeError as error:
        return incumbent, {
            "free_vertices": len(free_vertices),
            "before": before,
            "after": None,
            "accepted": False,
            "message": str(error),
        }
    candidate = result["labels"]
    if not is_feasible(graph, candidate, k) or sum(candidate) >= before:
        raise AssertionError("restricted MIP did not return a strict feasible improvement")
    return candidate, {
        "free_vertices": len(free_vertices),
        "before": before,
        "after": sum(candidate),
        "accepted": True,
        "optimal_subproblem": result["optimal"],
        "relative_gap": result["relative_gap"],
    }


def hybrid_refine(
    graph,
    k: int,
    incumbent: list[int],
    rng: random.Random,
    deadline: float,
    subproblem_time: float,
    size_fractions: Sequence[float],
    max_subproblems: int,
    failure_window: int = 1,
) -> tuple[list[int], list[dict]]:
    """Standalone V2 intensification, also used by the unit test."""
    if not is_feasible(graph, incumbent, k):
        raise ValueError("GA incumbent is infeasible")
    sizes = sorted({min(graph.n, max(1, math.ceil(value * graph.n)))
                    for value in size_fractions})
    level = 0
    failures = 0
    history: list[dict] = []
    for iteration in range(1, max_subproblems + 1):
        remaining = deadline - time.monotonic()
        if remaining <= 0.05:
            break
        incumbent, record = restricted_improvement(
            graph, k, incumbent, sizes[level], rng,
            max(0.05, min(subproblem_time, remaining)),
        )
        record.update({"iteration": iteration, "neighborhood_level": level})
        history.append(record)
        if record["accepted"]:
            level = 0
            failures = 0
        else:
            failures += 1
            if failures >= failure_window:
                level = level + 1 if level + 1 < len(sizes) else 0
                failures = 0
    return incumbent, history


def parse_fractions(value: str) -> list[float]:
    fractions = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not fractions or any(item <= 0 or item > 1 for item in fractions):
        raise argparse.ArgumentTypeError("fractions must lie in (0,1]")
    return fractions


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("instance_dir", type=Path)
    parser.add_argument("k", type=int)
    parser.add_argument("--ga-exe", type=Path)
    parser.add_argument("--seed", type=int, default=13579)
    parser.add_argument("--time-limit", type=float, default=60.0)
    parser.add_argument("--ga-time", type=float, default=30.0)
    parser.add_argument("--subproblem-time", type=float, default=5.0)
    parser.add_argument("--population", type=int, default=50)
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--epoch-generations", type=int, default=10)
    parser.add_argument("--init", choices=("random", "greedy"), default="greedy")
    parser.add_argument(
        "--neighborhood-fractions", type=parse_fractions,
        default=parse_fractions("0.20,0.30,0.40,0.60,1.00"),
    )
    parser.add_argument("--failure-window", type=int, default=1)
    parser.add_argument("--max-subproblems", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.time_limit <= 0 or args.ga_time <= 0 or args.subproblem_time <= 0:
        parser.error("time limits must be positive")
    if args.max_subproblems <= 0 or args.epoch_generations <= 0:
        parser.error("iteration limits must be positive")
    if args.failure_window <= 0:
        parser.error("failure-window must be positive")

    project_root = Path(__file__).resolve().parents[1]
    ga_exe = args.ga_exe or project_root / "bin" / "kmdrdp_ga.exe"
    graph = load_instance(args.instance_dir)
    start = time.monotonic()
    deadline = start + args.time_limit
    rng = random.Random(args.seed)
    incumbent, _ = run_ga(
        ga_exe,
        args.instance_dir,
        args.k,
        args.seed,
        min(args.ga_time, args.time_limit),
        args.population,
        args.generations,
        args.init,
    )
    initial_ga_objective = sum(incumbent)
    history: list[dict] = []
    sizes = sorted({min(graph.n, max(1, math.ceil(value * graph.n)))
                    for value in args.neighborhood_fractions})
    level = 0
    failures = 0
    subproblem_count = 0

    while subproblem_count < args.max_subproblems and time.monotonic() < deadline:
        remaining = deadline - time.monotonic()
        if remaining <= 0.05:
            break
        incumbent, record = restricted_improvement(
            graph, args.k, incumbent, sizes[level], rng,
            max(0.05, min(args.subproblem_time, remaining)),
        )
        subproblem_count += 1
        record.update({
            "iteration": subproblem_count,
            "epoch": None,
            "neighborhood_level": level,
        })
        history.append(record)
        if record["accepted"]:
            remaining = deadline - time.monotonic()
            if remaining > 0.1:
                reinjected, _ = run_ga(
                    ga_exe,
                    args.instance_dir,
                    args.k,
                    args.seed + subproblem_count * 104729,
                    min(max(0.1, args.ga_time * 0.10), remaining),
                    args.population,
                    args.epoch_generations,
                    args.init,
                    incumbent,
                )
                record["after_reinjection"] = sum(reinjected)
                if sum(reinjected) < sum(incumbent):
                    incumbent = reinjected
            level = 0
            failures = 0
        else:
            failures += 1
            if failures >= args.failure_window:
                level = level + 1 if level + 1 < len(sizes) else 0
                failures = 0

    violations = layer_violations(graph, incumbent)
    report = {
        "algorithm": f"support-aware-GA-{args.init[0].upper()}-ILP-v2",
        "instance": str(args.instance_dir.resolve()),
        "n": graph.n,
        "layers": graph.m,
        "k": args.k,
        "seed": args.seed,
        "initial_ga_objective": initial_ga_objective,
        "objective": sum(incumbent),
        "improvement": initial_ga_objective - sum(incumbent),
        "valid_layers": [index + 1 for index, value in enumerate(violations)
                         if value == 0],
        "feasible": is_feasible(graph, incumbent, args.k),
        "elapsed_seconds": time.monotonic() - start,
        "failure_window": args.failure_window,
        "max_subproblems": args.max_subproblems,
        "neighborhood_sizes": sizes,
        "subproblems": history,
        "labels": incumbent,
    }
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
