#!/usr/bin/env python3
"""Summarize a sparse-term sensitivity check from archived matched-GO results."""

from __future__ import annotations

import csv
import json
from pathlib import Path


SET_NAMES = ("R2_strict", "C3_strict", "C2_strict")
MINIMUM_CANDIDATE_COUNT = 5


def jaccard(first: set[str], second: set[str]) -> dict[str, float | int]:
    intersection = len(first & second)
    union = len(first | second)
    return {
        "intersection": intersection,
        "union": union,
        "jaccard": round(intersection / union, 3),
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    protocol = json.loads((biology / "ppi5_degree_matched_go_protocol.json").read_text(encoding="utf-8"))
    significant = {name: set() for name in SET_NAMES}
    filtered = {name: set() for name in SET_NAMES}
    singleton_counts = {name: 0 for name in SET_NAMES}
    with (biology / "ppi5_degree_matched_go.csv").open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            name = row["set"]
            if name not in significant or row["significant_fdr_0_05"] != "True":
                continue
            term = row["go_id"]
            if term in significant[name]:
                raise ValueError(f"duplicate GO term for {name}: {term}")
            significant[name].add(term)
            count = int(row["candidate_count"])
            if count == 1:
                singleton_counts[name] += 1
            if count >= MINIMUM_CANDIDATE_COUNT:
                filtered[name].add(term)

    summary = {
        "source": "Results/biology/ppi5_degree_matched_go.csv",
        "minimum_annotated_candidate_proteins": MINIMUM_CANDIDATE_COUNT,
        "original_bh_adjustment_retained": True,
        "bh_eligible_background_terms": protocol["eligible_terms"],
        "sets": {
            name: {
                "significant_all": len(significant[name]),
                "significant_with_minimum_count": len(filtered[name]),
                "single_candidate_significant_all": singleton_counts[name],
            }
            for name in SET_NAMES
        },
        "comparisons": {
            f"R2_strict_vs_{name}": {
                "all": jaccard(significant["R2_strict"], significant[name]),
                "minimum_count": jaccard(filtered["R2_strict"], filtered[name]),
            }
            for name in ("C3_strict", "C2_strict")
        },
    }
    output = biology / "ppi5_degree_matched_go_sparse_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
