#!/usr/bin/env python3
"""Compare consensus proteins with 10,000 five-layer degree-matched null sets."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


SETS = ("C3_strict", "C2_strict", "C23_strict",
        "C3_8of10", "C2_8of10", "C23_8of10")
METRICS = ("degree", "harmonic_closeness", "betweenness")
TRIALS = 10_000
SEED = 20_261_006
NEIGHBORS = 20
DISTANCE_TEMPERATURE = 0.01


def read_table(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    consensus = read_table(biology / "ppi5_k5_protein_consensus.tsv", "\t")
    n = len(consensus)
    if n != 4223:
        raise ValueError(f"expected 4223 proteins, found {n}")

    centrality = np.empty((n, 5, 3), dtype=np.float64)
    for layer in range(1, 6):
        rows = read_table(biology / f"ppi5_layer_{layer}_centrality.csv")
        if len(rows) != n or [int(row["internal_id"]) for row in rows] != list(range(1, n + 1)):
            raise ValueError(f"centrality ID mismatch in layer {layer}")
        for metric_index, metric in enumerate(METRICS):
            centrality[:, layer - 1, metric_index] = [float(row[metric]) for row in rows]

    degree_features = np.log1p(centrality[:, :, 0])
    feature_mean = degree_features.mean(axis=0)
    feature_std = degree_features.std(axis=0, ddof=1)
    standardized = (degree_features - feature_mean) / feature_std
    metric_matrix = centrality.reshape(n, -1)
    rng = np.random.default_rng(SEED)
    output_rows = []
    diagnostics = {}

    for set_name in SETS:
        selected = np.array([row[set_name] == "1" for row in consensus])
        case_ids = np.flatnonzero(selected)
        control_ids = np.flatnonzero(~selected)
        tree = cKDTree(standardized[control_ids])
        k_neighbors = min(NEIGHBORS, len(control_ids))
        distances, neighbor_positions = tree.query(standardized[case_ids], k=k_neighbors)
        if distances.ndim == 1:
            distances = distances[:, None]
        if neighbor_positions.ndim == 1:
            neighbor_positions = neighbor_positions[:, None]
        neighbor_ids = control_ids[neighbor_positions]
        weights = np.exp(-(distances - distances[:, [0]]) / DISTANCE_TEMPERATURE)
        cumulative_weights = np.cumsum(weights / weights.sum(axis=1, keepdims=True), axis=1)
        null_means = np.empty((TRIALS, metric_matrix.shape[1]), dtype=np.float64)
        matched_feature_means = np.empty((TRIALS, 5), dtype=np.float64)

        for trial in range(TRIALS):
            draws = rng.random(len(case_ids))
            choices = (cumulative_weights < draws[:, None]).sum(axis=1)
            sampled = neighbor_ids[np.arange(len(case_ids)), choices]
            null_means[trial] = metric_matrix[sampled].mean(axis=0)
            matched_feature_means[trial] = degree_features[sampled].mean(axis=0)

        observed = metric_matrix[case_ids].mean(axis=0)
        background = metric_matrix[control_ids].mean(axis=0)
        null_mean = null_means.mean(axis=0)
        null_std = null_means.std(axis=0, ddof=1)
        empirical_p = np.empty_like(observed)
        for column in range(len(observed)):
            distance = abs(observed[column] - null_mean[column])
            empirical_p[column] = (
                1 + np.count_nonzero(abs(null_means[:, column] - null_mean[column]) >= distance)
            ) / (TRIALS + 1)

        for layer in range(5):
            for metric_index, metric in enumerate(METRICS):
                column = layer * len(METRICS) + metric_index
                output_rows.append({
                    "set": set_name,
                    "set_size": len(case_ids),
                    "layer": layer + 1,
                    "metric": metric,
                    "candidate_mean": observed[column],
                    "unmatched_background_mean": background[column],
                    "degree_matched_null_mean": null_mean[column],
                    "degree_matched_null_sd": null_std[column],
                    "candidate_to_null_ratio": (
                        observed[column] / null_mean[column] if null_mean[column] else ""
                    ),
                    "z_score": (
                        (observed[column] - null_mean[column]) / null_std[column]
                        if null_std[column] else ""
                    ),
                    "empirical_p_two_sided": empirical_p[column],
                    "trials": TRIALS,
                })

        # Standardized mean differences on the five matching variables.
        case_feature_mean = degree_features[case_ids].mean(axis=0)
        matched_feature_mean = matched_feature_means.mean(axis=0)
        smd = (case_feature_mean - matched_feature_mean) / feature_std
        diagnostics[set_name] = {
            "set_size": int(len(case_ids)),
            "max_abs_standardized_mean_difference_log_degree": float(np.max(np.abs(smd))),
            "layerwise_standardized_mean_difference_log_degree": smd.tolist(),
        }
        print(set_name, diagnostics[set_name], flush=True)

    output_file = biology / "ppi5_degree_matched_centrality.csv"
    with output_file.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    report = {
        "trials": TRIALS,
        "seed": SEED,
        "matching_features": "standardized log1p degree in each of five layers",
        "sampling": "with replacement from 20 nearest non-set proteins; exponentially distance weighted",
        "distance_temperature": DISTANCE_TEMPERATURE,
        "diagnostics": diagnostics,
    }
    (biology / "ppi5_degree_matched_diagnostics.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
