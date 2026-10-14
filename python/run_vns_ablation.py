"""Frozen exploratory ablation; separate artifacts, never alters formal results.

Usage: python/run_vns_ablation.py --prepare | --run | --status
Requires Python 3. All data below are outputs of actual executable runs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "formal_experiments" / "vns_ablation_protocol.json"
OUTPUT = ROOT / "Results" / "ablation" / "vns_small_20261013"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2)
    temporary.replace(path)


def load_graph(directory):
    layers = []
    n = None
    files = sorted(directory.glob("layer_*.txt"),
                   key=lambda path: int(path.stem.split("_")[1]))
    for path in files:
        with path.open() as stream:
            size, declared = map(int, stream.readline().split())
            if n is None:
                n = size
            assert size == n
            adjacency = [set() for _ in range(n)]
            count = 0
            for line in stream:
                if not line.strip():
                    continue
                u, v = (int(x) - 1 for x in line.split())
                assert 0 <= u < n and 0 <= v < n
                if u == v:
                    continue
                adjacency[u].add(v)
                adjacency[v].add(u)
                count += 1
            assert count == declared
            layers.append(adjacency)
    assert layers and n == 1000
    return layers


def replay(labels, layers, k):
    assert len(labels) == len(layers[0]) and all(x in (0, 1, 2, 3) for x in labels)
    violations = []
    for adjacency in layers:
        invalid = 0
        for i, neighbors in enumerate(adjacency):
            if labels[i] not in (0, 1):
                continue
            d2 = sum(labels[j] == 2 for j in neighbors)
            d3 = sum(labels[j] == 3 for j in neighbors)
            invalid += (labels[i] == 0 and d3 == 0 and d2 < 2) or (
                labels[i] == 1 and d2 + d3 == 0)
        violations.append(invalid)
    assert sum(x == 0 for x in violations) >= k
    return {"objective": sum(labels), "valid_layers": sum(x == 0 for x in violations),
            "layer_violations": violations, "validation_errors": 0}


def prepare():
    assert not (OUTPUT / "manifest.json").exists(), "Already prepared; do not overwrite frozen protocol"
    protocol = json.loads(PROTOCOL.read_text())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "runs").mkdir(exist_ok=True)
    (OUTPUT / "initial").mkdir(exist_ok=True)
    source_paths = ["include/kmdrdp/algorithms.hpp", "include/kmdrdp/evaluator.hpp",
                    "include/kmdrdp/graph.hpp", "src/graph.cpp", "src/evaluator.cpp",
                    "src/algorithms.cpp", "apps/vns_ablation_main.cpp", "apps/cli_util.hpp",
                    "apps/evaluation_benchmark_main.cpp", "tests/test_ablation.cpp",
                    "python/run_vns_ablation.py", "formal_experiments/vns_ablation_protocol.json",
                    "formal_experiments/build_vns_ablation.ps1"]
    hashes = {}
    for relative in source_paths:
        target = OUTPUT / "source_snapshot" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
        hashes[relative] = digest(target)
    for name in ("vns_ablation.exe", "evaluation_benchmark.exe"):
        shutil.copyfile(ROOT / "bin" / name, OUTPUT / name)
        hashes[name] = digest(OUTPUT / name)
    with (ROOT / protocol["initial_ledger"]).open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    configurations = []
    for index, config in enumerate(protocol["configurations"]):
        directory = ROOT.parent / "k-MLRDP-main" / "Instances" / config["family"] / config["instance"]
        matches = [r for r in rows if Path(r["instance"]).name == config["instance"]
                   and int(r["k"]) == config["k"] and r["method"] == "greedy"
                   and int(r["seed"]) == 24680]
        assert len(matches) == 1
        labels = [int(x) for x in matches[0]["labels"].split()]
        layers = load_graph(directory)
        validation = replay(labels, layers, config["k"])
        assert validation["objective"] == int(matches[0]["objective"])
        initial_path = OUTPUT / "initial" / (str(index) + ".txt")
        initial_path.write_text(" ".join(map(str, labels)) + "\n", encoding="ascii")
        configurations.append({**config, "index": index, "directory": str(directory),
                               "layers": len(layers), "initial_file": str(initial_path),
                               "initial_sha256": digest(initial_path),
                               "initial_validation": validation,
                               "layer_sha256": {p.name: digest(p) for p in directory.glob("layer_*.txt")}})
    blocks = [(c["index"], seed) for c in configurations for seed in protocol["seeds"]]
    random.Random(protocol["schedule_seed"]).shuffle(blocks)
    schedule = []
    for i, (index, seed) in enumerate(blocks):
        variants = protocol["variants"]
        order = variants[i % 3:] + variants[:i % 3]
        schedule.extend({"config_index": index, "seed": seed, "variant": variant}
                        for variant in order)
    manifest = {"protocol": protocol, "protocol_sha256": digest(PROTOCOL),
                "initial_ledger_sha256": digest(ROOT / protocol["initial_ledger"]),
                "configurations": configurations, "schedule": schedule, "artifact_sha256": hashes,
                "created_local_date": "2026-10-13", "prepared_unix_time": time.time()}
    write_json(OUTPUT / "manifest.json", manifest)
    write_json(OUTPUT / "status.json", {"state": "prepared", "completed": 0, "expected": 54})
    print("Prepared immutable manifest: 6 replayed initial vectors, 54 scheduled searches.", flush=True)


def completed_rows(manifest):
    rows = []
    for item in manifest["schedule"]:
        path = OUTPUT / "runs" / f'{item["config_index"]}_{item["seed"]}_{item["variant"]}.json'
        if path.exists():
            row = json.loads(path.read_text())
            assert row["protocol_sha256"] == manifest["protocol_sha256"]
            assert row["validation"]["validation_errors"] == 0
            rows.append(row)
    return rows


def summarize(manifest, rows):
    variants = {}
    for variant in manifest["protocol"]["variants"]:
        selected = [r for r in rows if r["variant"] == variant]
        if not selected:
            continue
        variants[variant] = {
            "runs": len(selected),
            "mean_objective": statistics.mean(r["objective"] for r in selected),
            "mean_elapsed_seconds": statistics.mean(r["elapsed_seconds"] for r in selected),
            "repair_calls": sum(r["repair_calls"] for r in selected),
            "repair_steps": sum(r["repair_steps"] for r in selected),
            "fallback_count": sum(r["fallback_count"] for r in selected),
            "repairs_without_fallback": sum(r["repairs_without_fallback"] for r in selected),
            "accepted_ties": sum(r["accepted_ties"] for r in selected),
            "mean_evaluations_per_second": statistics.mean(
                (r["repair_evaluations"] + r["local_search_evaluations"]) / r["elapsed_seconds"]
                for r in selected),
        }
    lookup = {(r["config_index"], r["seed"], r["variant"]): r for r in rows}
    comparisons = {}
    for variant in ("no-repair-deficit", "fixed-tie"):
        wins = ties = losses = 0
        improvements = []
        config_records = []
        for config in manifest["configurations"]:
            per_config = []
            for seed in manifest["protocol"]["seeds"]:
                baseline = lookup.get((config["index"], seed, "full"))
                ablated = lookup.get((config["index"], seed, variant))
                if not baseline or not ablated:
                    continue
                difference = ablated["objective"] - baseline["objective"]
                wins += difference > 0
                ties += difference == 0
                losses += difference < 0
                improvement = 100 * difference / ablated["objective"]
                improvements.append(improvement)
                per_config.append(improvement)
            if per_config:
                config_records.append({"instance": config["instance"], "k": config["k"],
                                       "paired_seeds": len(per_config),
                                       "mean_full_relative_improvement_percent": statistics.mean(per_config)})
        comparisons[variant] = {"full_wins_ties_losses": [wins, ties, losses],
                                "paired_runs": len(improvements),
                                "mean_full_relative_improvement_percent": statistics.mean(improvements) if improvements else None,
                                "by_configuration": config_records}
    summary = {"completed": len(rows), "expected": 54, "validation_errors": 0,
               "interpretation": "Exploratory six-graph screen with three seeds; no confirmatory significance claims",
               "variants": variants, "paired_comparisons": comparisons}
    write_json(OUTPUT / "summary.json", summary)
    return summary


def run():
    manifest = json.loads((OUTPUT / "manifest.json").read_text())
    assert digest(PROTOCOL) == manifest["protocol_sha256"], "Protocol changed"
    for name in ("vns_ablation.exe", "evaluation_benchmark.exe"):
        assert digest(OUTPUT / name) == manifest["artifact_sha256"][name]
    graphs = {}
    for config in manifest["configurations"]:
        directory = Path(config["directory"])
        assert all(digest(directory / name) == expected for name, expected in config["layer_sha256"].items())
        assert digest(Path(config["initial_file"])) == config["initial_sha256"]
        graphs[config["index"]] = load_graph(directory)
    started = time.time()
    try:
        write_json(OUTPUT / "status.json", {"state": "running_benchmark", "pid": os.getpid(),
                                          "completed": len(completed_rows(manifest)), "expected": 54,
                                          "started_unix_time": started})
        benchmark_path = OUTPUT / "evaluation_benchmark.json"
        if not benchmark_path.exists():
            benchmarks = []
            for config in manifest["configurations"]:
                settings = manifest["protocol"]["incremental_benchmark"]
                command = [str(OUTPUT / "evaluation_benchmark.exe"), config["directory"],
                           str(config["k"]), str(settings["seed"]), str(settings["probes"]),
                           str(settings["repeats"]), config["initial_file"]]
                result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=180)
                benchmarks.append({"instance": config["instance"], **json.loads(result.stdout)})
                print("BENCHMARK", config["instance"], benchmarks[-1]["speedup_full_over_incremental"], flush=True)
            write_json(benchmark_path, {"results": benchmarks,
                                       "scope": "Candidate evaluator throughput, not full VNS acceleration",
                                       "total_verified_comparisons": sum(r["verified_comparisons"] for r in benchmarks),
                                       "mismatches": sum(r["mismatches"] for r in benchmarks)})
        for item in manifest["schedule"]:
            config = manifest["configurations"][item["config_index"]]
            key = f'{item["config_index"]}_{item["seed"]}_{item["variant"]}'
            path = OUTPUT / "runs" / (key + ".json")
            if path.exists():
                continue
            rows = completed_rows(manifest)
            write_json(OUTPUT / "status.json", {"state": "running_search", "pid": os.getpid(),
                                              "completed": len(rows), "expected": 54,
                                              "current": item, "started_unix_time": started,
                                              "updated_unix_time": time.time()})
            command = [str(OUTPUT / "vns_ablation.exe"), config["directory"], str(config["k"]),
                       str(item["seed"]), str(manifest["protocol"]["time_limit_seconds"]),
                       item["variant"], config["initial_file"]]
            result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=300)
            row = json.loads(result.stdout)
            assert (row["variant"], row["seed"], row["k"], row["n"], row["layers"]) == (
                item["variant"], item["seed"], config["k"], 1000, config["layers"])
            validation = replay(row["labels"], graphs[config["index"]], config["k"])
            assert validation["objective"] == row["objective"]
            assert validation["valid_layers"] == row["valid_layers"]
            assert row["objective"] <= config["initial_validation"]["objective"]
            assert row["initial_objective"] == config["initial_validation"]["objective"]
            if item["variant"] == "fixed-tie":
                assert row["final_tie_probability"] == 0.4
            row.update({"config_index": config["index"], "instance": config["instance"],
                        "family": config["family"], "protocol_sha256": manifest["protocol_sha256"],
                        "initial_sha256": config["initial_sha256"], "validation": validation,
                        "finished_unix_time": time.time()})
            write_json(path, row)
            print("RUN", key, "objective", row["objective"], "seconds", row["elapsed_seconds"], flush=True)
            summarize(manifest, completed_rows(manifest))
        rows = completed_rows(manifest)
        assert len(rows) == 54
        summary = summarize(manifest, rows)
        write_json(OUTPUT / "status.json", {"state": "complete", "completed": 54, "expected": 54,
                                          "validation_errors": 0, "pid": os.getpid(),
                                          "started_unix_time": started, "finished_unix_time": time.time()})
        print("COMPLETE", json.dumps(summary["paired_comparisons"]), flush=True)
    except Exception as error:
        write_json(OUTPUT / "status.json", {"state": "failed", "pid": os.getpid(),
                                          "completed": len(completed_rows(manifest)), "expected": 54,
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
        print((OUTPUT / "status.json").read_text())
