#!/usr/bin/env python3
"""Exact DR-New1 and DR-New2 models for k-multiplex double Roman domination."""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix


@dataclass(frozen=True)
class MultiplexGraph:
    n: int
    layers: tuple[tuple[frozenset[int], ...], ...]
    name: str

    @property
    def m(self) -> int:
        return len(self.layers)


def load_instance(directory: str | Path) -> MultiplexGraph:
    directory = Path(directory)
    pattern = re.compile(r"layer_(\d+)\.txt$", re.IGNORECASE)
    files = sorted(
        (path for path in directory.iterdir() if path.is_file() and pattern.match(path.name)),
        key=lambda path: int(pattern.match(path.name).group(1)),
    )
    if not files:
        raise ValueError(f"no layer_*.txt files in {directory}")
    layers: list[tuple[frozenset[int], ...]] = []
    common_n: int | None = None
    for path in files:
        tokens = path.read_text(encoding="utf-8").split()
        if len(tokens) < 2:
            raise ValueError(f"invalid header in {path}")
        n, declared_edges = map(int, tokens[:2])
        edge_tokens = list(map(int, tokens[2:]))
        if len(edge_tokens) % 2:
            raise ValueError(f"incomplete edge in {path}")
        if len(edge_tokens) // 2 != declared_edges:
            raise ValueError(f"edge count does not match header in {path}")
        if common_n is None:
            common_n = n
        elif n != common_n:
            raise ValueError("all layers must have the same vertex set")
        adjacency = [set() for _ in range(n)]
        for u, v in zip(edge_tokens[0::2], edge_tokens[1::2]):
            u -= 1
            v -= 1
            if not (0 <= u < n and 0 <= v < n):
                raise ValueError(f"edge endpoint outside 1..n in {path}")
            if u != v:
                adjacency[u].add(v)
                adjacency[v].add(u)
        layers.append(tuple(frozenset(neighbors) for neighbors in adjacency))
    return MultiplexGraph(common_n, tuple(layers), directory.name)


def layer_violations(graph: MultiplexGraph, labels: Sequence[int]) -> list[int]:
    if len(labels) != graph.n or any(label not in (0, 1, 2, 3) for label in labels):
        raise ValueError("invalid label vector")
    counts: list[int] = []
    for layer in graph.layers:
        violations = 0
        for vertex, neighbors in enumerate(layer):
            count2 = sum(labels[neighbor] == 2 for neighbor in neighbors)
            count3 = sum(labels[neighbor] == 3 for neighbor in neighbors)
            if labels[vertex] == 0 and count3 == 0 and count2 < 2:
                violations += 1
            elif labels[vertex] == 1 and count2 + count3 == 0:
                violations += 1
        counts.append(violations)
    return counts


def is_feasible(graph: MultiplexGraph, labels: Sequence[int], k: int) -> bool:
    return sum(value == 0 for value in layer_violations(graph, labels)) >= k


def _append_row(rows, lower, upper, coefficients, lb=-math.inf, ub=math.inf):
    rows.append(coefficients)
    lower.append(lb)
    upper.append(ub)


def _validate_fixed_labels(
    graph: MultiplexGraph, fixed_labels: Sequence[int | None] | None
) -> None:
    if fixed_labels is None:
        return
    if len(fixed_labels) != graph.n:
        raise ValueError("fixed_labels must contain one entry per vertex")
    if any(value is not None and value not in (0, 1, 2, 3) for value in fixed_labels):
        raise ValueError("fixed labels must be 0, 1, 2, 3, or None")


def build_model(
    graph: MultiplexGraph,
    k: int,
    model: str,
    fixed_labels: Sequence[int | None] | None = None,
    objective_upper_bound: int | None = None,
):
    if not 1 <= k <= graph.m:
        raise ValueError("k must be in 1..M")
    if model not in {"direct", "cumulative"}:
        raise ValueError("model must be direct or cumulative")
    _validate_fixed_labels(graph, fixed_labels)
    if objective_upper_bound is not None and objective_upper_bound < 0:
        raise ValueError("objective_upper_bound must be nonnegative")
    variable_count = 3 * graph.n + graph.m
    z = lambda layer: 3 * graph.n + layer
    rows: list[dict[int, float]] = []
    lower: list[float] = []
    upper: list[float] = []
    objective = np.zeros(variable_count)

    for vertex in range(graph.n):
        objective[3 * vertex : 3 * vertex + 3] = (1, 2, 3) if model == "direct" else (1, 1, 1)

    if model == "direct":
        x = lambda vertex, label: 3 * vertex + label - 1
        for vertex in range(graph.n):
            _append_row(rows, lower, upper,
                        {x(vertex, 1): 1, x(vertex, 2): 1, x(vertex, 3): 1}, ub=1)
        for layer_index, layer in enumerate(graph.layers):
            for vertex, neighbors in enumerate(layer):
                zero = {x(vertex, 1): 2, x(vertex, 2): 2, x(vertex, 3): 2,
                        z(layer_index): -2}
                one = {x(vertex, 1): -1, z(layer_index): -1}
                for neighbor in neighbors:
                    zero[x(neighbor, 2)] = zero.get(x(neighbor, 2), 0) + 1
                    zero[x(neighbor, 3)] = zero.get(x(neighbor, 3), 0) + 2
                    one[x(neighbor, 2)] = one.get(x(neighbor, 2), 0) + 1
                    one[x(neighbor, 3)] = one.get(x(neighbor, 3), 0) + 1
                _append_row(rows, lower, upper, zero, lb=0)
                _append_row(rows, lower, upper, one, lb=-1)

        if fixed_labels is not None:
            for vertex, label in enumerate(fixed_labels):
                if label is None:
                    continue
                for candidate in (1, 2, 3):
                    expected = int(label == candidate)
                    _append_row(
                        rows, lower, upper,
                        {x(vertex, candidate): 1}, lb=expected, ub=expected,
                    )
    else:
        a = lambda vertex: 3 * vertex
        b = lambda vertex: 3 * vertex + 1
        c = lambda vertex: 3 * vertex + 2
        for vertex in range(graph.n):
            _append_row(rows, lower, upper, {c(vertex): 1, b(vertex): -1}, ub=0)
            _append_row(rows, lower, upper, {b(vertex): 1, a(vertex): -1}, ub=0)
        for layer_index, layer in enumerate(graph.layers):
            for vertex, neighbors in enumerate(layer):
                zero = {a(vertex): 2, z(layer_index): -2}
                one = {a(vertex): -1, b(vertex): 1, z(layer_index): -1}
                for neighbor in neighbors:
                    zero[b(neighbor)] = zero.get(b(neighbor), 0) + 1
                    zero[c(neighbor)] = zero.get(c(neighbor), 0) + 1
                    one[b(neighbor)] = one.get(b(neighbor), 0) + 1
                _append_row(rows, lower, upper, zero, lb=0)
                _append_row(rows, lower, upper, one, lb=-1)

        if fixed_labels is not None:
            for vertex, label in enumerate(fixed_labels):
                if label is None:
                    continue
                expected = (int(label >= 1), int(label >= 2), int(label == 3))
                for variable, value in zip(
                    (a(vertex), b(vertex), c(vertex)), expected
                ):
                    _append_row(
                        rows, lower, upper,
                        {variable: 1}, lb=value, ub=value,
                    )

    _append_row(rows, lower, upper, {z(layer): 1 for layer in range(graph.m)}, lb=k)
    if objective_upper_bound is not None:
        _append_row(
            rows, lower, upper,
            {column: float(objective[column])
             for column in range(3 * graph.n) if objective[column] != 0},
            ub=objective_upper_bound,
        )
    matrix = lil_matrix((len(rows), variable_count), dtype=float)
    for row_index, coefficients in enumerate(rows):
        for column, value in coefficients.items():
            matrix[row_index, column] = value
    return (
        objective,
        np.ones(variable_count, dtype=int),
        Bounds(np.zeros(variable_count), np.ones(variable_count)),
        LinearConstraint(matrix.tocsr(), np.asarray(lower), np.asarray(upper)),
    )


def solve_exact(
    graph: MultiplexGraph,
    k: int,
    model: str,
    time_limit: float | None = None,
    fixed_labels: Sequence[int | None] | None = None,
    objective_upper_bound: int | None = None,
):
    objective, integrality, bounds, constraints = build_model(
        graph, k, model, fixed_labels, objective_upper_bound
    )
    options = {"presolve": True}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    result = milp(
        c=objective,
        integrality=integrality,
        bounds=bounds,
        constraints=constraints,
        options=options,
    )
    if result.x is None:
        raise RuntimeError(f"solver returned no incumbent: {result.message}")
    rounded = np.rint(result.x).astype(int)
    if model == "direct":
        labels = []
        for vertex in range(graph.n):
            active = [label for label in (1, 2, 3) if rounded[3 * vertex + label - 1]]
            labels.append(active[0] if active else 0)
    else:
        labels = [int(sum(rounded[3 * vertex : 3 * vertex + 3])) for vertex in range(graph.n)]
    if not is_feasible(graph, labels, k):
        raise AssertionError("MILP incumbent failed independent DRD replay")
    return {
        "model": "DR-New1" if model == "direct" else "DR-New2",
        "status": int(result.status),
        "message": str(result.message),
        "optimal": int(result.status) == 0,
        "objective": int(sum(labels)),
        "dual_bound": None if getattr(result, "mip_node_count", None) is None else float(result.mip_dual_bound),
        "relative_gap": None if getattr(result, "mip_gap", None) is None else float(result.mip_gap),
        "valid_layers": [index + 1 for index, value in enumerate(layer_violations(graph, labels)) if value == 0],
        "labels": labels,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("instance_dir")
    parser.add_argument("k", type=int)
    parser.add_argument("--model", choices=("direct", "cumulative"), default="direct")
    parser.add_argument("--time-limit", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    graph = load_instance(args.instance_dir)
    report = {
        "instance": str(Path(args.instance_dir).resolve()),
        "n": graph.n,
        "layers": graph.m,
        "k": args.k,
        **solve_exact(graph, args.k, args.model, args.time_limit),
    }
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

