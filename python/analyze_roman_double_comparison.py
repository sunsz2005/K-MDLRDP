#!/usr/bin/env python3
"""Compare ordinary-Roman label 2 with double-Roman labels 3 and 2."""

from __future__ import annotations

import csv
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact


SETS = ("R2_strict", "C3_strict", "C2_strict", "R2_8of10", "C3_8of10", "C2_8of10")
METRICS = ("degree", "harmonic_closeness", "betweenness")
TRIALS = 10_000
NEIGHBORS = 20
TEMPERATURE = 0.01


def read_table(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    double = read_table(biology / "ppi5_k5_protein_consensus.tsv", "\t")
    roman = read_table(biology / "ppi5_roman_k5_protein_consensus.tsv", "\t")
    if len(double) != 4223 or len(roman) != 4223:
        raise ValueError("expected 4223 proteins in both consensus tables")
    if [r["systematic_name"] for r in double] != [r["systematic_name"] for r in roman]:
        raise ValueError("ordinary- and double-Roman protein mappings differ")
    names = np.array([r["systematic_name"] for r in double])
    masks = {
        "R2_strict": np.array([r["R2_strict"] == "1" for r in roman]),
        "R2_8of10": np.array([r["R2_8of10"] == "1" for r in roman]),
    }
    for name in SETS:
        if name.startswith("C"):
            masks[name] = np.array([r[name] == "1" for r in double])

    centrality = np.empty((4223, 5, 3), dtype=float)
    for layer in range(1, 6):
        rows = read_table(biology / f"ppi5_layer_{layer}_centrality.csv")
        for j, metric in enumerate(METRICS):
            centrality[:, layer - 1, j] = [float(r[metric]) for r in rows]
    features_raw = np.log1p(centrality[:, :, 0])
    feature_sd = features_raw.std(axis=0, ddof=1)
    features = (features_raw - features_raw.mean(axis=0)) / feature_sd
    metric_matrix = centrality.reshape(4223, -1)

    essential = set((biology / "essential_primary_systematic_S288C_null_inviable.txt")
                    .read_text(encoding="utf-8").split())
    essential_vector = np.isin(names, list(essential))
    centrality_output = []
    essential_output = []
    diagnostics = {}

    for set_index, set_name in enumerate(SETS):
        selected = masks[set_name]
        cases = np.flatnonzero(selected)
        controls = np.flatnonzero(~selected)
        distances, positions = cKDTree(features[controls]).query(
            features[cases], k=min(NEIGHBORS, len(controls)))
        if distances.ndim == 1:
            distances, positions = distances[:, None], positions[:, None]
        neighbors = controls[positions]
        weights = np.exp(-(distances - distances[:, [0]]) / TEMPERATURE)
        cumulative = np.cumsum(weights / weights.sum(axis=1, keepdims=True), axis=1)
        rng = np.random.default_rng(20_261_100 + set_index)
        null_metrics = np.empty((TRIALS, 15), dtype=float)
        null_essential = np.empty(TRIALS, dtype=int)
        matched_features = np.empty((TRIALS, 5), dtype=float)
        for trial in range(TRIALS):
            draws = rng.random(len(cases))
            choices = (cumulative < draws[:, None]).sum(axis=1)
            sampled = neighbors[np.arange(len(cases)), choices]
            null_metrics[trial] = metric_matrix[sampled].mean(axis=0)
            null_essential[trial] = essential_vector[sampled].sum()
            matched_features[trial] = features_raw[sampled].mean(axis=0)

        observed = metric_matrix[cases].mean(axis=0)
        null_mean = null_metrics.mean(axis=0)
        for layer in range(5):
            for metric_index, metric in enumerate(METRICS):
                column = layer * 3 + metric_index
                distance = abs(observed[column] - null_mean[column])
                p_value = (1 + np.count_nonzero(
                    abs(null_metrics[:, column] - null_mean[column]) >= distance
                )) / (TRIALS + 1)
                centrality_output.append({
                    "set": set_name,
                    "set_size": len(cases),
                    "layer": layer + 1,
                    "metric": metric,
                    "candidate_mean": observed[column],
                    "degree_matched_null_mean": null_mean[column],
                    "candidate_to_null_ratio": observed[column] / null_mean[column],
                    "empirical_p_two_sided": p_value,
                    "trials": TRIALS,
                })

        observed_essential = int(essential_vector[cases].sum())
        table = [[observed_essential, len(cases) - observed_essential],
                 [int(essential_vector[controls].sum()), len(controls) - int(essential_vector[controls].sum())]]
        odds, fisher_p = fisher_exact(table)
        p_enrich = (1 + np.count_nonzero(null_essential >= observed_essential)) / (TRIALS + 1)
        p_deplete = (1 + np.count_nonzero(null_essential <= observed_essential)) / (TRIALS + 1)
        essential_output.append({
            "set": set_name,
            "set_size": len(cases),
            "essential_count": observed_essential,
            "essential_fraction": observed_essential / len(cases),
            "fisher_odds_ratio": odds,
            "fisher_p_two_sided": fisher_p,
            "degree_matched_null_mean_count": null_essential.mean(),
            "degree_matched_empirical_p_two_sided": min(1.0, 2 * min(p_enrich, p_deplete)),
            "direction_vs_matched_null": "enriched" if observed_essential > null_essential.mean() else "depleted",
            "trials": TRIALS,
        })
        smd = (features_raw[cases].mean(axis=0) - matched_features.mean(axis=0)) / feature_sd
        diagnostics[set_name] = {
            "size": int(len(cases)),
            "max_abs_smd_log_degree": float(np.max(np.abs(smd))),
            "layerwise_smd_log_degree": smd.tolist(),
        }
        print(set_name, diagnostics[set_name], flush=True)

    overlap_output = []
    for left, right in combinations(SETS, 2):
        if ("strict" in left) != ("strict" in right):
            continue
        intersection = int(np.count_nonzero(masks[left] & masks[right]))
        union = int(np.count_nonzero(masks[left] | masks[right]))
        overlap_output.append({
            "left": left,
            "right": right,
            "left_size": int(masks[left].sum()),
            "right_size": int(masks[right].sum()),
            "intersection": intersection,
            "jaccard": intersection / union if union else 0,
            "left_contained_fraction": intersection / masks[left].sum(),
            "right_contained_fraction": intersection / masks[right].sum(),
        })

    for filename, rows in (
        ("ppi5_roman_double_centrality.csv", centrality_output),
        ("ppi5_roman_double_essential.csv", essential_output),
        ("ppi5_roman_double_overlap.csv", overlap_output),
    ):
        with (biology / filename).open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (biology / "ppi5_roman_double_diagnostics.json").write_text(
        json.dumps({"trials": TRIALS, "diagnostics": diagnostics}, indent=2) + "\n",
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
