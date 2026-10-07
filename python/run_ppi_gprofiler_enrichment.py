#!/usr/bin/env python3
"""Run version-recorded g:Profiler enrichment with the 4,223-protein background."""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import requests


API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SETS = ("C3_strict", "C2_strict", "C23_strict",
        "C3_8of10", "C2_8of10", "C23_8of10")
SOURCES = ("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    with (biology / "ppi5_k5_protein_consensus.tsv").open(
            newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    background = [row["systematic_name"] for row in rows]
    flat = []
    run_metadata = {}

    for set_name in SETS:
        query = [row["systematic_name"] for row in rows if row[set_name] == "1"]
        request_data = {
            "organism": "scerevisiae",
            "query": query,
            "sources": list(SOURCES),
            "user_threshold": 0.05,
            "significance_threshold_method": "fdr",
            "domain_scope": "custom",
            "background": background,
            "all_results": False,
            "ordered": False,
            "no_evidences": False,
        }
        (biology / f"gprofiler_{set_name}_request.json").write_text(
            json.dumps(request_data, indent=2) + "\n", encoding="utf-8"
        )
        last_error = None
        for attempt in range(1, 4):
            try:
                response = requests.post(API, json=request_data, timeout=180)
                response.raise_for_status()
                data = response.json()
                break
            except Exception as error:
                last_error = error
                if attempt == 3:
                    raise
                time.sleep(5 * attempt)
        else:
            raise RuntimeError(last_error)
        (biology / f"gprofiler_{set_name}_response.json").write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        run_metadata[set_name] = data.get("meta", {})
        mapped_genes = (data.get("meta", {}).get("genes_metadata", {})
                        .get("query", {}).get("query_1", {}).get("ensgs", []))
        for result in data.get("result", []):
            evidences = result.get("intersections", [])
            intersecting_genes = [
                gene for gene, evidence in zip(mapped_genes, evidences) if evidence
            ]
            flat.append({
                "set": set_name,
                "set_size": len(query),
                "source": result.get("source", ""),
                "term_id": result.get("native", ""),
                "term_name": result.get("name", ""),
                "adjusted_p_fdr": result.get("p_value", ""),
                "significant": result.get("significant", ""),
                "query_size": result.get("query_size", ""),
                "term_size": result.get("term_size", ""),
                "intersection_size": result.get("intersection_size", ""),
                "effective_domain_size": result.get("effective_domain_size", ""),
                "precision": result.get("precision", ""),
                "recall": result.get("recall", ""),
                "intersection": ";".join(intersecting_genes),
            })
        print(set_name, "query", len(query), "significant terms", len(data.get("result", [])), flush=True)

    fields = list(flat[0]) if flat else ["set", "set_size", "source", "term_id", "term_name"]
    with (biology / "ppi5_functional_enrichment.csv").open(
            "w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)
    (biology / "gprofiler_run_metadata.json").write_text(
        json.dumps({"api": API, "sources": SOURCES, "sets": run_metadata},
                   indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
