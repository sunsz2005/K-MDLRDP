#!/usr/bin/env python3
"""Stream the Reactome all-level Ensembl mapping and retain yeast pathways."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import requests


MAPPING_URL = "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"
VERSION_URL = "https://reactome.org/ContentService/data/database/version"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    source = root / "Results" / "biology" / "source"
    source.mkdir(parents=True, exist_ok=True)
    output = source / "Reactome_Ensembl2Reactome_All_Levels_SCE_2026-10-06.tsv"

    response = requests.get(MAPPING_URL, stream=True, timeout=(30, 300))
    response.raise_for_status()
    kept = 0
    scanned = 0
    digest = hashlib.sha256()
    with output.open("wb") as stream:
        for line in response.iter_lines():
            scanned += 1
            if b"\tR-SCE-" not in line and not line.endswith(b"\tSaccharomyces cerevisiae"):
                continue
            payload = line + b"\n"
            stream.write(payload)
            digest.update(payload)
            kept += 1
            if kept % 10000 == 0:
                print("kept", kept, "scanned", scanned, flush=True)

    version_response = requests.get(VERSION_URL, timeout=60)
    version_response.raise_for_status()
    metadata = {
        "source_url": MAPPING_URL,
        "download_date": "2026-10-06",
        "reactome_version": version_response.text.strip(),
        "last_modified": response.headers.get("Last-Modified"),
        "etag": response.headers.get("ETag"),
        "scanned_lines": scanned,
        "retained_lines": kept,
        "filtered_file_sha256": digest.hexdigest(),
        "filtered_file": str(output),
    }
    (source / "Reactome_SCE_download_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
