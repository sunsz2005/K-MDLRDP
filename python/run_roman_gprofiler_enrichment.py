#!/usr/bin/env python3
"""Run g:Profiler for ordinary-Roman label-2 consensus sets."""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import requests


API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SETS = ("R2_strict", "R2_8of10")
SOURCES = ("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    with (biology / "ppi5_roman_k5_protein_consensus.tsv").open(
            newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    background = [row["systematic_name"] for row in rows]
    flat, metadata = [], {}
    for set_name in SETS:
        query = [row["systematic_name"] for row in rows if row[set_name] == "1"]
        payload = {
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
            json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        for attempt in range(1, 4):
            try:
                response = requests.post(API, json=payload, timeout=180)
                response.raise_for_status()
                data = response.json()
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(5 * attempt)
        (biology / f"gprofiler_{set_name}_response.json").write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        metadata[set_name] = data.get("meta", {})
        mapped = (data.get("meta", {}).get("genes_metadata", {})
                  .get("query", {}).get("query_1", {}).get("ensgs", []))
        for result in data.get("result", []):
            intersection = [gene for gene, evidence in
                            zip(mapped, result.get("intersections", [])) if evidence]
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
                "intersection": ";".join(intersection),
            })
        print(set_name, len(query), len(data.get("result", [])), flush=True)
    fields = list(flat[0]) if flat else ["set", "set_size", "source", "term_id", "term_name"]
    with (biology / "ppi5_roman_functional_enrichment.csv").open(
            "w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)
    (biology / "gprofiler_roman_run_metadata.json").write_text(
        json.dumps({"api": API, "sources": SOURCES, "sets": metadata},
                   indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
