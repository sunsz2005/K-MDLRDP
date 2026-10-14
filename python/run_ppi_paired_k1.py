"""Complete Greedy/VNS seed pairs without changing archived Hybrid results.

Commands: --prepare, --run, --status. Outputs are separate JSON certificates.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import sys
import time
import traceback

from exact_ilp import load_instance, layer_violations
from run_experiments import run_one

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "formal_experiments" / "ppi_paired_k1_protocol.json"
OUT = ROOT / "Results" / "paired_ppi" / "k1_10seeds_20261013"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def key(row):
    return (Path(row["instance"]).name, row["method"], int(row["seed"]))


def validate(row, graph):
    labels = row["labels"]
    if isinstance(labels, str):
        labels = [int(x) for x in labels.split()]
    if row.get("status") != "ok" or row.get("feasible") not in (True, "True"):
        raise ValueError("unsuccessful solver record")
    counts = layer_violations(graph, labels)
    valid = sum(x == 0 for x in counts)
    if sum(labels) != int(row["objective"]) or valid != int(row["valid_layers"]) or valid < 1:
        raise ValueError("independent label replay failed")
    return {"labels": labels, "objective": sum(labels), "valid_layers": valid,
            "layer_violations": counts, "validation_errors": 0}


def prepare():
    if (OUT / "manifest.json").exists():
        raise RuntimeError("Already prepared; do not overwrite frozen study")
    protocol = json.loads(PROTOCOL.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "runs").mkdir(exist_ok=True)
    (OUT / "bin").mkdir(exist_ok=True)
    artifacts = {}
    for relative in ("python/run_ppi_paired_k1.py", "python/run_experiments.py",
                     "python/exact_ilp.py", "formal_experiments/ppi_paired_k1_protocol.json",
                     "src/algorithms.cpp", "src/evaluator.cpp", "src/graph.cpp",
                     "include/kmdrdp/algorithms.hpp", "include/kmdrdp/evaluator.hpp",
                     "include/kmdrdp/graph.hpp", "apps/vns_main.cpp", "apps/greedy_main.cpp",
                     "apps/cli_util.hpp"):
        target = OUT / "source_snapshot" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
        artifacts[relative] = sha(target)
    for name in ("kmdrdp_greedy.exe", "kmdrdp_vns.exe"):
        shutil.copyfile(ROOT / "bin" / name, OUT / "bin" / name)
        artifacts["bin/" + name] = sha(OUT / "bin" / name)
    inputs = {}
    for instance in protocol["instances"]:
        directory = ROOT / "Instances" / "PPI" / instance
        inputs[instance] = {"directory": str(directory),
                            "layer_sha256": {p.name: sha(p) for p in directory.glob("layer_*.txt")}}
    reused = []
    source_hashes = {}
    for source_name, relative in protocol["reuse"].items():
        path = ROOT / relative
        source_hashes[relative] = sha(path)
        with path.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                if int(row["k"]) != 1:
                    continue
                if Path(row["instance"]).name not in protocol["instances"]:
                    raise ValueError("unexpected archived instance")
                row = dict(row, provenance=relative)
                for field in ("n", "layers", "k", "seed", "objective", "valid_layers"):
                    row[field] = int(row[field])
                row["elapsed_seconds"] = float(row["elapsed_seconds"])
                reused.append(row)
    expected_reused = {(i, "ga-g-hybrid", s) for i in protocol["instances"] for s in protocol["seeds"]}
    expected_reused |= {(i, m, 13579) for i in protocol["instances"] for m in protocol["new_methods"]}
    if len(reused) != 36 or {key(r) for r in reused} != expected_reused:
        raise ValueError("archived key set mismatch")
    # Validate by graph to avoid retaining several large PPI graphs at once.
    for instance in protocol["instances"]:
        graph = load_instance(inputs[instance]["directory"])
        for row in reused:
            if Path(row["instance"]).name == instance:
                checked = validate(row, graph)
                row["labels"] = checked.pop("labels")
                row["validation"] = checked
        print("ARCHIVE VALIDATED", instance, flush=True)
        del graph
    schedule = [{"instance": i, "seed": s, "method": m}
                for i in protocol["instances"] for s in protocol["new_seeds"]
                for m in protocol["new_methods"]]
    manifest = {"protocol": protocol, "protocol_sha256": sha(PROTOCOL),
                "inputs": inputs, "source_ledger_sha256": source_hashes,
                "artifact_sha256": artifacts, "schedule": schedule,
                "prepared_unix_time": time.time()}
    save(OUT / "reused_rows.json", reused)
    manifest["reused_rows_sha256"] = sha(OUT / "reused_rows.json")
    save(OUT / "manifest.json", manifest)
    save(OUT / "status.json", {"state": "prepared", "new_completed": 0,
                              "new_expected": 54, "reused_validated": 36})
    print("PREPARED: 36 archived certificates replayed; 54 new runs scheduled.", flush=True)


def completed(manifest):
    rows = []
    for item in manifest["schedule"]:
        path = OUT / "runs" / f'{item["instance"]}_{item["method"]}_{item["seed"]}.json'
        if path.exists():
            row = json.loads(path.read_text())
            if row["protocol_sha256"] != manifest["protocol_sha256"] or row["validation"]["validation_errors"]:
                raise ValueError("invalid completed record")
            rows.append(row)
    return rows


def summarize(manifest, new_rows):
    rows = json.loads((OUT / "reused_rows.json").read_text()) + new_rows
    indexed = {key(r): r for r in rows}
    assert len(indexed) == len(rows)
    summaries = []
    for instance in manifest["protocol"]["instances"]:
        methods = {}
        for method in ("greedy", "vns-g-fi", "ga-g-hybrid"):
            selected = [r for r in rows if key(r)[:2] == (instance, method)]
            objectives = [r["objective"] for r in selected]
            methods[method] = {"runs": len(selected), "mean_objective": statistics.mean(objectives),
                               "sample_sd": statistics.stdev(objectives) if len(objectives) > 1 else None,
                               "range": [min(objectives), max(objectives)],
                               "mean_elapsed_seconds": statistics.mean(r["elapsed_seconds"] for r in selected)}
        comparisons = {}
        for first, second in (("ga-g-hybrid", "greedy"), ("ga-g-hybrid", "vns-g-fi"),
                              ("vns-g-fi", "greedy")):
            records = []
            for seed in manifest["protocol"]["seeds"]:
                a = indexed.get((instance, first, seed))
                b = indexed.get((instance, second, seed))
                if a and b:
                    records.append({"seed": seed, "difference": b["objective"] - a["objective"],
                                    "relative_improvement_percent": 100 * (b["objective"]-a["objective"]) / b["objective"]})
            def describe(paired):
                return {"paired_runs": len(paired),
                        "first_wins_ties_losses": [sum(r["difference"] > 0 for r in paired),
                                                  sum(r["difference"] == 0 for r in paired),
                                                  sum(r["difference"] < 0 for r in paired)],
                        "mean_paired_relative_improvement_percent": statistics.mean(
                            r["relative_improvement_percent"] for r in paired) if paired else None}
            comparisons[first + "_vs_" + second] = {"all_ten": describe(records),
                                                      "gate_seed_excluded": describe([r for r in records if r["seed"] != 13579])}
        summaries.append({"instance": instance, "methods": methods, "comparisons": comparisons})
    summary = {"new_completed": len(new_rows), "new_expected": 54, "combined_rows": len(rows),
               "combined_expected": 90, "validation_errors": 0,
               "complete": len(rows) == 90, "summary": summaries,
               "interpretation": "Matched algorithm seeds on three nested reconstructed PPI networks; gate reuse disclosed; descriptive runtime comparisons"}
    save(OUT / "summary.json", summary)
    if len(rows) == 90:
        save(OUT / "combined_90run_ledger.json", [indexed[k] for k in sorted(indexed)])
    return summary


def run():
    manifest = json.loads((OUT / "manifest.json").read_text())
    assert sha(PROTOCOL) == manifest["protocol_sha256"]
    assert sha(OUT / "reused_rows.json") == manifest["reused_rows_sha256"]
    for name in ("kmdrdp_greedy.exe", "kmdrdp_vns.exe"):
        assert sha(OUT / "bin" / name) == manifest["artifact_sha256"]["bin/" + name]
    for info in manifest["inputs"].values():
        assert all(sha(Path(info["directory"]) / name) == expected
                   for name, expected in info["layer_sha256"].items())
    started = time.time()
    graph = None
    last_instance = None
    try:
        for item in manifest["schedule"]:
            path = OUT / "runs" / f'{item["instance"]}_{item["method"]}_{item["seed"]}.json'
            if path.exists():
                continue
            count = len(completed(manifest))
            save(OUT / "status.json", {"state": "running", "pid": os.getpid(),
                                      "new_completed": count, "new_expected": 54,
                                      "reused_validated": 36, "current": item,
                                      "started_unix_time": started, "updated_unix_time": time.time()})
            directory = Path(manifest["inputs"][item["instance"]]["directory"])
            if item["instance"] != last_instance:
                del graph
                graph = load_instance(directory)
                last_instance = item["instance"]
            print("RUN", item, flush=True)
            record = run_one(OUT, sys.executable, item["method"], directory, 1,
                             item["seed"], manifest["protocol"]["time_limit_seconds"])
            row = {**record, "instance": str(directory), "method": item["method"],
                   "seed": item["seed"], "n": graph.n, "layers": graph.m, "k": 1,
                   "protocol_sha256": manifest["protocol_sha256"],
                   "finished_unix_time": time.time()}
            checked = validate(row, graph)
            row["labels"] = checked.pop("labels")
            row["validation"] = checked
            save(path, row)
            print("VALIDATED", item, "objective", row["objective"], "seconds", row["elapsed_seconds"], flush=True)
            summarize(manifest, completed(manifest))
        rows = completed(manifest)
        assert len(rows) == 54
        summary = summarize(manifest, rows)
        assert summary["complete"]
        save(OUT / "status.json", {"state": "complete", "new_completed": 54,
                                  "new_expected": 54, "combined_rows": 90,
                                  "validation_errors": 0, "pid": os.getpid(),
                                  "started_unix_time": started, "finished_unix_time": time.time()})
        print("COMPLETE: 90 unique validated paired records.", flush=True)
    except Exception as error:
        save(OUT / "status.json", {"state": "failed", "pid": os.getpid(),
                                  "new_completed": len(completed(manifest)), "new_expected": 54,
                                  "error": str(error), "traceback": traceback.format_exc()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()
    assert sum((args.prepare, args.run, args.status)) == 1
    if args.prepare:
        prepare()
    elif args.run:
        run()
    else:
        print((OUT / "status.json").read_text())
