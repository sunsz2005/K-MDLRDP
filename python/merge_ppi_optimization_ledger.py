#!/usr/bin/env python3
"""Merge the PPI gate and formal k=1 runs into one deduplicated ledger."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORMAL = ROOT / "Results" / "formal"
GATE = FORMAL / "ppi_gate_merged.csv"
CONFIRMATION = FORMAL / "ppi_formal_ga_g_hybrid_k1_10seeds.csv"
OUTPUT = FORMAL / "ppi_optimization_63run_ledger.csv"
MANIFEST = FORMAL / "ppi_optimization_63run_ledger_manifest.json"


def instance_name(value: str) -> str:
    return value.replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]


def key(row: dict[str, str]) -> tuple[str, int, str, int]:
    return (
        instance_name(row["instance"]),
        int(row["k"]),
        row["method"],
        int(row["seed"]),
    )


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def validate_row(row: dict[str, str]) -> list[str]:
    errors: list[str] = []
    labels = [int(x) for x in row["labels"].split()] if row["labels"].strip() else []
    if row["status"].lower() != "ok":
        errors.append("status is not ok")
    if row["feasible"].lower() != "true":
        errors.append("feasible is not true")
    if len(labels) != int(row["n"]):
        errors.append("label-vector length differs from n")
    if any(label not in (0, 1, 2, 3) for label in labels):
        errors.append("label outside {0,1,2,3}")
    if sum(labels) != int(float(row["objective"])):
        errors.append("label sum differs from objective")
    if int(row["valid_layers"]) < int(row["k"]):
        errors.append("valid_layers is below k")
    return errors


def main() -> None:
    gate_header, gate_rows = read_csv(GATE)
    confirmation_header, confirmation_rows = read_csv(CONFIRMATION)
    if gate_header != confirmation_header:
        raise RuntimeError("PPI source ledgers have different schemas")

    merged: dict[tuple[str, int, str, int], dict[str, str]] = {}
    sources: dict[tuple[str, int, str, int], list[str]] = {}
    for source, rows in (("gate", gate_rows), ("formal_k1", confirmation_rows)):
        for row in rows:
            row_key = key(row)
            if row_key in merged and merged[row_key] != row:
                raise RuntimeError(f"Conflicting duplicate row: {row_key}")
            merged[row_key] = row
            sources.setdefault(row_key, []).append(source)

    duplicate_keys = [row_key for row_key, values in sources.items() if len(values) > 1]
    validation_errors: list[dict[str, object]] = []
    output_rows: list[dict[str, str]] = []
    for row_key in sorted(merged, key=lambda x: (int(x[0].split("-")[1][:-6]), x[1], x[2], x[3])):
        row = dict(merged[row_key])
        row["instance"] = f"Instances/PPI/{row_key[0]}"
        errors = validate_row(row)
        if errors:
            validation_errors.append({"key": list(row_key), "errors": errors})
        output_rows.append(row)

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=gate_header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    manifest = {
        "description": "Deduplicated PPI optimization ledger used by the manuscript",
        "merge_key": ["instance", "k", "method", "seed"],
        "source_files": {
            str(GATE.relative_to(ROOT)).replace("\\", "/"): len(gate_rows),
            str(CONFIRMATION.relative_to(ROOT)).replace("\\", "/"): len(confirmation_rows),
        },
        "source_rows_total": len(gate_rows) + len(confirmation_rows),
        "duplicate_rows_removed": len(duplicate_keys),
        "duplicate_keys": [list(value) for value in duplicate_keys],
        "output_file": str(OUTPUT.relative_to(ROOT)).replace("\\", "/"),
        "output_rows": len(output_rows),
        "normalization": "Absolute instance paths were replaced by Instances/PPI/<instance>.",
        "structural_validation": {
            "checks": [
                "status=ok",
                "feasible=true",
                "label count=n",
                "labels in {0,1,2,3}",
                "sum(labels)=objective",
                "valid_layers>=k",
            ],
            "errors": validation_errors,
        },
        "source_feasibility_replay": {
            "gate_summary": "Results/formal/ppi_gate_summary.json",
            "formal_k1_summary": "Results/formal/ppi_formal_ga_g_hybrid_k1_summary.json",
            "reported_errors": 0,
        },
    }
    with MANIFEST.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    if len(output_rows) != 63 or validation_errors:
        raise RuntimeError(
            f"Expected 63 valid rows; got {len(output_rows)} with "
            f"{len(validation_errors)} validation errors"
        )
    print(f"Wrote {len(output_rows)} rows to {OUTPUT}")
    print(f"Removed {len(duplicate_keys)} identical duplicate rows")


if __name__ == "__main__":
    main()
