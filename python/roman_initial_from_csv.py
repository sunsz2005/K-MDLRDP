#!/usr/bin/env python3
"""Convert one reference ordinary-Roman greedy CSV row to a label vector."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_nodes(value: str) -> set[int]:
    return {int(token) for token in (value or "").split()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_labels", type=Path)
    parser.add_argument("--row", type=int, default=None,
                        help="zero-based row; default selects the minimum objective")
    args = parser.parse_args()

    with args.input_csv.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if args.row is None:
        row_index, row = min(
            enumerate(rows), key=lambda item: float(item[1]["greedy_obj"])
        )
    else:
        if not 0 <= args.row < len(rows):
            raise ValueError(f"row {args.row} is outside 0..{len(rows) - 1}")
        row_index, row = args.row, rows[args.row]
    n = int(row["n"])
    label2 = parse_nodes(row["label2_nodes"])
    label1 = parse_nodes(row["label1_nodes"])
    if label2 & label1:
        raise ValueError("label-1 and label-2 sets overlap")
    if not label1 | label2 <= set(range(1, n + 1)):
        raise ValueError("node ID outside 1..n")
    labels = [2 if i in label2 else 1 if i in label1 else 0 for i in range(1, n + 1)]
    objective = 2 * len(label2) + len(label1)
    recorded = int(float(row["greedy_obj"]))
    if objective != recorded:
        raise ValueError(f"computed objective {objective} != recorded {recorded}")
    args.output_labels.parent.mkdir(parents=True, exist_ok=True)
    args.output_labels.write_text(" ".join(map(str, labels)) + "\n", encoding="ascii")
    print(f"wrote row={row_index}, n={n}, objective={objective} to {args.output_labels}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
