#!/usr/bin/env python3
"""Test essential-protein enrichment with five-layer degree-matched controls."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2


SETS = ("C3_strict", "C2_strict", "C23_strict",
        "C3_8of10", "C2_8of10", "C23_8of10")
TRIALS = 10_000
SEED = 20_261_007
NEIGHBORS = 20
DISTANCE_TEMPERATURE = 0.01


def read_table(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    phenotype_path = biology / "source" / "SGD_phenotype_data_2026-10-06.tab"
    consensus = read_table(biology / "ppi5_k5_protein_consensus.tsv", "\t")
    names = np.array([row["systematic_name"] for row in consensus])
    n = len(names)

    definitions = {
        "primary_systematic_S288C_null_inviable": set(),
        "sensitivity_systematic_all_strains_null_inviable": set(),
        "sensitivity_all_null_inviable": set(),
    }
    with phenotype_path.open(newline="", encoding="utf-8", errors="replace") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if len(row) != 14 or row[6] != "null" or row[9] != "inviable":
                continue
            definitions["sensitivity_all_null_inviable"].add(row[0])
            if row[5] == "systematic mutation set":
                definitions["sensitivity_systematic_all_strains_null_inviable"].add(row[0])
                if row[8].upper() == "S288C":
                    definitions["primary_systematic_S288C_null_inviable"].add(row[0])

    background = set(names)
    for definition, proteins in definitions.items():
        (biology / f"essential_{definition}.txt").write_text(
            "\n".join(sorted(proteins & background)) + "\n", encoding="utf-8"
        )

    degrees = np.empty((n, 5), dtype=np.float64)
    for layer in range(1, 6):
        rows = read_table(biology / f"ppi5_layer_{layer}_centrality.csv")
        degrees[:, layer - 1] = [float(row["degree"]) for row in rows]
    features = np.log1p(degrees)
    features = (features - features.mean(axis=0)) / features.std(axis=0, ddof=1)
    essential_matrix = np.column_stack([
        np.isin(names, list(proteins)) for proteins in definitions.values()
    ])
    rng = np.random.default_rng(SEED)
    output = []

    for set_name in SETS:
        selected = np.array([row[set_name] == "1" for row in consensus])
        case_ids = np.flatnonzero(selected)
        control_ids = np.flatnonzero(~selected)
        distances, neighbor_positions = cKDTree(features[control_ids]).query(
            features[case_ids], k=min(NEIGHBORS, len(control_ids))
        )
        if distances.ndim == 1:
            distances = distances[:, None]
        if neighbor_positions.ndim == 1:
            neighbor_positions = neighbor_positions[:, None]
        neighbor_ids = control_ids[neighbor_positions]
        weights = np.exp(-(distances - distances[:, [0]]) / DISTANCE_TEMPERATURE)
        cumulative_weights = np.cumsum(weights / weights.sum(axis=1, keepdims=True), axis=1)
        null_counts = np.empty((TRIALS, len(definitions)), dtype=np.int32)
        for trial in range(TRIALS):
            draws = rng.random(len(case_ids))
            choices = (cumulative_weights < draws[:, None]).sum(axis=1)
            sampled = neighbor_ids[np.arange(len(case_ids)), choices]
            null_counts[trial] = essential_matrix[sampled].sum(axis=0)

        for definition_index, (definition, proteins) in enumerate(definitions.items()):
            essential = essential_matrix[:, definition_index]
            case_essential = int(essential[case_ids].sum())
            control_essential = int(essential[control_ids].sum())
            table = np.array([
                [case_essential, len(case_ids) - case_essential],
                [control_essential, len(control_ids) - control_essential],
            ])
            odds_ratio, fisher_p = fisher_exact(table, alternative="two-sided")
            ci_low, ci_high = Table2x2(table).oddsratio_confint()
            null = null_counts[:, definition_index]
            empirical_p_enrichment = (
                1 + np.count_nonzero(null >= case_essential)
            ) / (TRIALS + 1)
            empirical_p_depletion = (
                1 + np.count_nonzero(null <= case_essential)
            ) / (TRIALS + 1)
            empirical_p_two_sided = min(
                1.0, 2.0 * min(empirical_p_enrichment, empirical_p_depletion)
            )
            output.append({
                "set": set_name,
                "set_size": len(case_ids),
                "essential_definition": definition,
                "essential_in_full_SGD_file": len(proteins),
                "essential_in_4223_background": int(essential.sum()),
                "essential_in_set": case_essential,
                "set_fraction": case_essential / len(case_ids),
                "background_fraction": essential.mean(),
                "fisher_odds_ratio": odds_ratio,
                "odds_ratio_ci95_low": ci_low,
                "odds_ratio_ci95_high": ci_high,
                "fisher_p_two_sided": fisher_p,
                "degree_matched_null_mean_count": null.mean(),
                "degree_matched_empirical_p_enrichment": empirical_p_enrichment,
                "degree_matched_empirical_p_depletion": empirical_p_depletion,
                "degree_matched_empirical_p_two_sided": empirical_p_two_sided,
                "trials": TRIALS,
            })
        print(set_name, "done", flush=True)

    with (biology / "ppi5_essential_protein_tests.csv").open(
            "w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    report = {
        "source_url": "https://downloads.yeastgenome.org/curation/literature/phenotype_data.tab",
        "download_date": "2026-10-06",
        "source_sha256": hashlib.sha256(phenotype_path.read_bytes()).hexdigest(),
        "primary_definition": "systematic mutation set; null; S288C; inviable",
        "background_size": n,
        "trials": TRIALS,
        "seed": SEED,
        "degree_matching": "20 nearest controls with replacement; exponential distance weights",
        "distance_temperature": DISTANCE_TEMPERATURE,
        "essential_counts_full_file": {key: len(value) for key, value in definitions.items()},
        "essential_counts_background": {
            key: len(value & background) for key, value in definitions.items()
        },
    }
    (biology / "ppi5_essential_protocol.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
