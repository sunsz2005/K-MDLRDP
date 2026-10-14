"""Derive validated heuristic-to-certified-optimum gaps; no new solver runs."""
import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path

from exact_ilp import load_instance, layer_violations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference-instances",
        type=Path,
        help="Path to the reference k-MLRDP-main/Instances directory",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    reference_instances = args.reference_instances or root.parent / "k-MLRDP-main" / "Instances"
    source = root / "Results/formal/n30_runs.csv"
    rows = list(csv.DictReader(source.open(encoding="utf-8-sig", newline="")))
    methods = ["greedy", "vns-g-fi", "ga-g-hybrid"]
    selected = [r for r in rows if r["method"] in methods + ["dr-new1", "dr-new2"]]
    keys = [(r["instance"], int(r["k"]), r["method"], int(r["seed"])) for r in selected]
    assert len(keys) == len(set(keys)), "duplicate record"
    optima = {}
    exact = {}
    for r in selected:
        assert r["status"] == "ok" and r["feasible"] == "True"
        key = (r["instance"], int(r["k"]))
        if r["method"].startswith("dr-new"):
            assert r["optimal"] == "True" and float(r["relative_gap"]) == 0
            exact.setdefault(key, {})[r["method"]] = int(r["objective"])
    assert len(exact) == 24
    for key, pair in exact.items():
        assert set(pair) == {"dr-new1", "dr-new2"}
        assert pair["dr-new1"] == pair["dr-new2"] > 0
        optima[key] = pair["dr-new1"]
    for instance in sorted({r["instance"] for r in selected}):
        normalized = instance.replace("\\", "/")
        if normalized.startswith("REFERENCE_INSTANCES/"):
            instance_path = reference_instances / normalized.removeprefix("REFERENCE_INSTANCES/")
        else:
            instance_path = Path(instance)
        if not instance_path.is_dir():
            raise FileNotFoundError(
                f"Missing {instance_path}; pass --reference-instances to the reference instance directory"
            )
        graph = load_instance(instance_path)
        for r in (r for r in selected if r["instance"] == instance):
            labels = [int(x) for x in r["labels"].split()]
            assert len(labels) == graph.n and all(x in (0, 1, 2, 3) for x in labels)
            assert sum(labels) == int(r["objective"])
            valid = sum(x == 0 for x in layer_violations(graph, labels))
            assert valid == int(r["valid_layers"]) and valid >= int(r["k"])
        print("VALIDATED", Path(instance).name, flush=True)
    summaries = []
    seeds = {13579, 24680, 97531, 86420, 112358, 271828, 314159, 161803, 424242, 987654}
    for method in methods:
        method_rows = [r for r in selected if r["method"] == method]
        assert len(method_rows) == 240
        gaps = []
        configurations = []
        for key, optimum in sorted(optima.items()):
            group = [r for r in method_rows if (r["instance"], int(r["k"])) == key]
            assert {int(r["seed"]) for r in group} == seeds and len(group) == 10
            group_gaps = [(int(r["objective"]) - optimum) / optimum * 100 for r in group]
            assert min(group_gaps) >= 0
            gaps.extend(group_gaps)
            configurations.append({"instance": Path(key[0]).name, "k": key[1],
                                   "optimum": optimum, "mean_gap_percent": statistics.mean(group_gaps),
                                   "optimal_hits": group_gaps.count(0)})
        summaries.append({"method": method, "runs": len(gaps),
                          "mean_gap_percent": statistics.mean(gaps),
                          "max_gap_percent": max(gaps), "optimal_hits": gaps.count(0),
                          "configs_with_optimal_hit": sum(c["optimal_hits"] > 0 for c in configurations),
                          "configurations": configurations})
    report = {"source": str(source.relative_to(root)),
              "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "certified_configurations": 24, "validated_records": len(selected),
              "validation_errors": 0, "gap_definition": "100*(heuristic weight - optimum)/optimum",
              "scope": "Retrospective descriptive analysis of screening runs; not post-selection replication",
              "methods": summaries}
    output = root / "Results/formal/n30_optimal_gap_summary.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps([{k: v for k, v in r.items() if k != "configurations"} for r in summaries], indent=2))


if __name__ == "__main__":
    main()
