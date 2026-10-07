#!/usr/bin/env python3
"""Solve five-layer all-layer ordinary Roman domination with SciPy/HiGHS."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


def read_edges(path: Path, expected_n: int) -> list[tuple[int, int]]:
    with path.open(encoding="utf-8") as stream:
        n, m = map(int, next(stream).split())
        if n != expected_n:
            raise ValueError(f"{path}: expected {expected_n} nodes, found {n}")
        edges = []
        for line in stream:
            if line.strip():
                u, v = map(int, line.split())
                edges.append((u - 1, v - 1))
    if len(edges) != m:
        raise ValueError(f"{path}: header has {m} edges but read {len(edges)}")
    return edges


def build_model(root: Path, formulation: str):
    n, layers = 4223, 5
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    layer_dir = root / "Instances" / "PPI" / "PPI-5Layers"

    # Variables are (label2, label1) for direct and (at-least-1, label2)
    # for cumulative.  In either encoding, the first block contributes the
    # self coefficient and the label-2 block covers neighbors.
    if formulation == "direct":
        self_block, neighbor_block = 0, 0
        second_self_block = n
        objective = np.r_[np.full(n, 2.0), np.ones(n)]
    else:
        self_block, neighbor_block = 0, n
        second_self_block = None
        objective = np.ones(2 * n)

    for layer in range(layers):
        offset = layer * n
        for vertex in range(n):
            rows.append(offset + vertex)
            cols.append(self_block + vertex)
            data.append(1.0)
            if second_self_block is not None:
                rows.append(offset + vertex)
                cols.append(second_self_block + vertex)
                data.append(1.0)
        for u, v in read_edges(layer_dir / f"layer_{layer + 1}.txt", n):
            rows.extend((offset + u, offset + v))
            cols.extend((neighbor_block + v, neighbor_block + u))
            data.extend((1.0, 1.0))

    relation_offset = layers * n
    if formulation == "direct":
        # label2_i + label1_i <= 1.
        for vertex in range(n):
            rows.extend((relation_offset + vertex, relation_offset + vertex))
            cols.extend((vertex, n + vertex))
            data.extend((1.0, 1.0))
        relation_upper = np.ones(n)
    else:
        # label2_i <= at_least_1_i.
        for vertex in range(n):
            rows.extend((relation_offset + vertex, relation_offset + vertex))
            cols.extend((n + vertex, vertex))
            data.extend((1.0, -1.0))
        relation_upper = np.zeros(n)

    matrix = coo_matrix((data, (rows, cols)), shape=(layers * n + n, 2 * n)).tocsr()
    lower = np.r_[np.ones(layers * n), np.full(n, -np.inf)]
    upper = np.r_[np.full(layers * n, np.inf), relation_upper]
    return objective, matrix, lower, upper


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formulation", choices=("direct", "cumulative"), required=True)
    parser.add_argument("--time-limit", type=float, default=3600.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    started = time.time()
    objective, matrix, lower, upper = build_model(root, args.formulation)
    built = time.time()
    result = milp(
        c=objective,
        integrality=np.ones_like(objective),
        bounds=Bounds(np.zeros_like(objective), np.ones_like(objective)),
        constraints=LinearConstraint(matrix, lower, upper),
        options={"time_limit": args.time_limit, "mip_rel_gap": 1e-5, "disp": True},
    )
    finished = time.time()
    if result.x is None:
        raise RuntimeError(f"solver returned no incumbent: {result.message}")
    n = 4223
    if args.formulation == "direct":
        label2 = result.x[:n] > 0.5
        label1 = result.x[n:] > 0.5
    else:
        at_least_1 = result.x[:n] > 0.5
        label2 = result.x[n:] > 0.5
        label1 = at_least_1 & ~label2
    labels = np.where(label2, 2, np.where(label1, 1, 0)).astype(int)
    computed_objective = int(labels.sum())
    layer_dir = root / "Instances" / "PPI" / "PPI-5Layers"
    violations = []
    for layer in range(1, 6):
        covered = np.zeros(n, dtype=bool)
        for u, v in read_edges(layer_dir / f"layer_{layer}.txt", n):
            if label2[v]:
                covered[u] = True
            if label2[u]:
                covered[v] = True
        violations.append(int(np.count_nonzero((labels == 0) & ~covered)))
    if any(violations):
        raise ValueError(f"infeasible reconstructed labeling; layer violations={violations}")
    payload = {
        "formulation": args.formulation,
        "solver": "scipy.optimize.milp / HiGHS",
        "status": int(result.status),
        "success": bool(result.success),
        "message": result.message,
        "objective": float(result.fun),
        "computed_integer_objective": computed_objective,
        "mip_node_count": getattr(result, "mip_node_count", None),
        "mip_dual_bound": getattr(result, "mip_dual_bound", None),
        "mip_gap": getattr(result, "mip_gap", None),
        "model_build_seconds": built - started,
        "solve_seconds": finished - built,
        "total_seconds": finished - started,
        "time_limit_seconds": args.time_limit,
        "label_counts": {str(value): int(np.count_nonzero(labels == value)) for value in range(3)},
        "zero_vertex_violations_by_layer": violations,
        "labels": " ".join(map(str, labels.tolist())),
        "label2_nodes": " ".join(map(str, (np.flatnonzero(label2) + 1).tolist())),
        "label1_nodes": " ".join(map(str, (np.flatnonzero(label1) + 1).tolist())),
    }
    if abs(float(result.fun) - computed_objective) > 1e-6:
        raise ValueError("reported and reconstructed objectives differ")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items()
                      if key not in {"labels", "label2_nodes", "label1_nodes"}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
