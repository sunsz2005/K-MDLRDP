#!/usr/bin/env python3
"""Empirical GO enrichment after matching five-layer degree profiles."""

from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix


SETS = ("R2_strict", "C3_strict", "C2_strict")
TRIALS = 10_000
NEIGHBORS = 20
TEMPERATURE = 0.01
MIN_BACKGROUND_COUNT = 5


def read_table(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def parse_obo(path: Path):
    terms: dict[str, dict] = {}
    current: dict | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line == "[Term]":
            if current and current.get("id"):
                terms[current["id"]] = current
            current = {"parents": []}
        elif line.startswith("["):
            if current and current.get("id"):
                terms[current["id"]] = current
            current = None
        elif current is not None:
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:]
            elif line.startswith("is_a: "):
                current["parents"].append(line[6:].split()[0])
            elif line.startswith("relationship: part_of "):
                current["parents"].append(line.split()[2])
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    if current and current.get("id"):
        terms[current["id"]] = current
    return {key: value for key, value in terms.items() if not value.get("obsolete")}


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values)
    running = 1.0
    total = len(p_values)
    for reverse_rank in range(total - 1, -1, -1):
        index = order[reverse_rank]
        rank = reverse_rank + 1
        running = min(running, p_values[index] * total / rank)
        adjusted[index] = running
    return adjusted


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    source = biology / "source"
    double = read_table(biology / "ppi5_k5_protein_consensus.tsv", "\t")
    roman = read_table(biology / "ppi5_roman_k5_protein_consensus.tsv", "\t")
    names = np.array([row["systematic_name"] for row in double])
    if len(names) != 4223 or list(names) != [row["systematic_name"] for row in roman]:
        raise ValueError("candidate tables do not share the expected protein mapping")
    masks = {
        "R2_strict": np.array([row["R2_strict"] == "1" for row in roman]),
        "C3_strict": np.array([row["C3_strict"] == "1" for row in double]),
        "C2_strict": np.array([row["C2_strict"] == "1" for row in double]),
    }

    terms = parse_obo(source / "go-basic_2026-08-08.obo")

    @lru_cache(maxsize=None)
    def ancestors(term_id: str) -> frozenset[str]:
        if term_id not in terms:
            return frozenset()
        result = {term_id}
        for parent in terms[term_id]["parents"]:
            result.update(ancestors(parent))
        return frozenset(result)

    name_to_index = {name: i for i, name in enumerate(names)}
    direct: dict[int, set[str]] = defaultdict(set)
    with gzip.open(source / "GO_sgd_2026-05-28.gaf.gz", "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15 or "NOT" in fields[3].split("|"):
                continue
            aliases = [fields[2], *fields[10].split("|")]
            protein = next((alias for alias in aliases if alias in name_to_index), None)
            if protein and fields[4] in terms:
                direct[name_to_index[protein]].add(fields[4])

    propagated: dict[int, set[str]] = {}
    background_counts: dict[str, int] = defaultdict(int)
    for protein_index, annotations in direct.items():
        expanded = set().union(*(ancestors(term) for term in annotations))
        propagated[protein_index] = expanded
        for term in expanded:
            background_counts[term] += 1
    eligible = sorted(term for term, count in background_counts.items()
                      if MIN_BACKGROUND_COUNT <= count < len(names))
    term_index = {term: i for i, term in enumerate(eligible)}
    matrix_rows, matrix_cols = [], []
    for protein_index, annotations in propagated.items():
        for term in annotations:
            if term in term_index:
                matrix_rows.append(protein_index)
                matrix_cols.append(term_index[term])
    annotation = csr_matrix(
        (np.ones(len(matrix_rows), dtype=np.int8), (matrix_rows, matrix_cols)),
        shape=(len(names), len(eligible)), dtype=np.int8)

    degrees = np.empty((len(names), 5), dtype=float)
    for layer in range(1, 6):
        rows = read_table(biology / f"ppi5_layer_{layer}_centrality.csv")
        degrees[:, layer - 1] = [float(row["degree"]) for row in rows]
    raw_features = np.log1p(degrees)
    feature_sd = raw_features.std(axis=0, ddof=1)
    features = (raw_features - raw_features.mean(axis=0)) / feature_sd

    output = []
    diagnostics = {}
    for set_index, set_name in enumerate(SETS):
        cases = np.flatnonzero(masks[set_name])
        controls = np.flatnonzero(~masks[set_name])
        observed = np.asarray(annotation[cases].sum(axis=0)).ravel().astype(int)
        distances, positions = cKDTree(features[controls]).query(
            features[cases], k=min(NEIGHBORS, len(controls)))
        if distances.ndim == 1:
            distances, positions = distances[:, None], positions[:, None]
        neighbors = controls[positions]
        weights = np.exp(-(distances - distances[:, [0]]) / TEMPERATURE)
        cumulative = np.cumsum(weights / weights.sum(axis=1, keepdims=True), axis=1)
        rng = np.random.default_rng(20_261_200 + set_index)
        null_sum = np.zeros(len(eligible), dtype=np.float64)
        null_ge = np.zeros(len(eligible), dtype=np.int32)
        matched_feature_sum = np.zeros(5, dtype=np.float64)
        for trial in range(TRIALS):
            draws = rng.random(len(cases))
            choices = (cumulative < draws[:, None]).sum(axis=1)
            sampled = neighbors[np.arange(len(cases)), choices]
            counts = np.asarray(annotation[sampled].sum(axis=0)).ravel()
            null_sum += counts
            null_ge += counts >= observed
            matched_feature_sum += raw_features[sampled].mean(axis=0)
        null_mean = null_sum / TRIALS
        p_values = (1 + null_ge) / (TRIALS + 1)
        q_values = bh_adjust(p_values)
        for j, term_id in enumerate(eligible):
            if observed[j] == 0:
                continue
            output.append({
                "set": set_name,
                "set_size": len(cases),
                "go_id": term_id,
                "term_name": terms[term_id].get("name", ""),
                "namespace": terms[term_id].get("namespace", ""),
                "candidate_count": observed[j],
                "background_count": background_counts[term_id],
                "degree_matched_null_mean_count": null_mean[j],
                "candidate_to_null_ratio": observed[j] / null_mean[j] if null_mean[j] else "",
                "empirical_p_enrichment": p_values[j],
                "bh_fdr_over_eligible_terms": q_values[j],
                "significant_fdr_0_05": q_values[j] < 0.05,
                "trials": TRIALS,
            })
        matched_mean = matched_feature_sum / TRIALS
        smd = (raw_features[cases].mean(axis=0) - matched_mean) / feature_sd
        diagnostics[set_name] = {
            "set_size": int(len(cases)),
            "significant_terms": int(np.count_nonzero((q_values < 0.05) & (observed > 0))),
            "max_abs_smd_log_degree": float(np.max(np.abs(smd))),
        }
        print(set_name, diagnostics[set_name], flush=True)

    output.sort(key=lambda row: (row["set"], float(row["bh_fdr_over_eligible_terms"]), row["go_id"]))
    with (biology / "ppi5_degree_matched_go.csv").open(
            "w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    protocol = {
        "sets": SETS,
        "trials": TRIALS,
        "random_seeds": {set_name: 20_261_200 + set_index
                         for set_index, set_name in enumerate(SETS)},
        "matching_temperature": TEMPERATURE,
        "matching_neighbors": NEIGHBORS,
        "matching": "standardized log1p degree in all five PPI layers",
        "sampling": "20 nearest non-set proteins with replacement and exponential distance weights",
        "go_propagation": "is_a and part_of ancestors",
        "negative_annotations": "NOT annotations excluded",
        "eligible_terms": len(eligible),
        "minimum_background_count": MIN_BACKGROUND_COUNT,
        "multiple_testing": "Benjamini-Hochberg across all eligible GO terms within each set",
        "diagnostics": diagnostics,
    }
    (biology / "ppi5_degree_matched_go_protocol.json").write_text(
        json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
