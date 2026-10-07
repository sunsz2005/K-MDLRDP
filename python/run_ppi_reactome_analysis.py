#!/usr/bin/env python3
"""Run Reactome's official over-representation service for consensus sets."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import requests


API = "https://reactome.org/AnalysisService/identifiers/"
VERSION_API = "https://reactome.org/ContentService/data/database/version"
SETS = ("C3_strict", "C2_strict", "C23_strict",
        "C3_8of10", "C2_8of10", "C23_8of10")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    biology = root / "Results" / "biology"
    with (biology / "ppi5_k5_protein_consensus.tsv").open(
            newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    version = requests.get(VERSION_API, timeout=60).text.strip()
    flat = []
    metadata = {}

    for set_name in SETS:
        query = [row["systematic_name"] for row in rows if row[set_name] == "1"]
        payload = "#Genes\n" + "\n".join(query) + "\n"
        params = {
            "pageSize": 5000,
            "page": 1,
            "sortBy": "ENTITIES_FDR",
            "order": "ASC",
            "resource": "TOTAL",
            "pValue": 1,
            "includeDisease": "false",
        }
        (biology / f"reactome_{set_name}_request.txt").write_text(payload, encoding="utf-8")
        response = requests.post(
            API, params=params, data=payload,
            headers={"Content-Type": "text/plain"}, timeout=300,
        )
        response.raise_for_status()
        data = response.json()
        (biology / f"reactome_{set_name}_response.json").write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        metadata[set_name] = {
            "query_size": len(query),
            "identifiers_not_found": data.get("identifiersNotFound"),
            "pathways_found": data.get("pathwaysFound"),
            "token": data.get("summary", {}).get("token"),
        }
        for pathway in data.get("pathways", []):
            entities = pathway.get("entities", {})
            species = pathway.get("species", {})
            if species.get("taxId") != "4932":
                continue
            flat.append({
                "set": set_name,
                "set_size": len(query),
                "pathway_id": pathway.get("stId", ""),
                "pathway_name": pathway.get("name", ""),
                "species": species.get("name", ""),
                "entities_found": entities.get("found", ""),
                "entities_total": entities.get("total", ""),
                "p_value": entities.get("pValue", ""),
                "fdr": entities.get("fdr", ""),
                "significant_fdr_0_05": float(entities.get("fdr", 1)) < 0.05,
                "lowest_level_pathway": pathway.get("llp", ""),
            })
        print(set_name, metadata[set_name], flush=True)

    with (biology / "ppi5_reactome_enrichment.csv").open(
            "w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    protocol = {
        "api": API,
        "reactome_version": version,
        "date": "2026-10-06",
        "background_note": (
            "Reactome AnalysisService uses its own mapped-entity universe; "
            "unlike the g:Profiler GO/KEGG analysis, this result does not use "
            "the custom 4,223-protein PPI background."
        ),
        "sets": metadata,
    }
    (biology / "ppi5_reactome_protocol.json").write_text(
        json.dumps(protocol, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(protocol, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
