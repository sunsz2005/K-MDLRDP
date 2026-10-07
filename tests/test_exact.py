from itertools import product
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from exact_ilp import is_feasible, load_instance, solve_exact  # noqa: E402


def brute_force(graph, k):
    best = None
    for labels in product(range(4), repeat=graph.n):
        value = sum(labels)
        if best is not None and value >= best:
            continue
        if is_feasible(graph, labels, k):
            best = value
    return best


def main():
    graph = load_instance(ROOT / "Instances" / "sample-2layer")
    for k in (1, 2):
        expected = brute_force(graph, k)
        direct = solve_exact(graph, k, "direct", 30)
        cumulative = solve_exact(graph, k, "cumulative", 30)
        assert direct["objective"] == expected
        assert cumulative["objective"] == expected
        assert is_feasible(graph, direct["labels"], k)
        assert is_feasible(graph, cumulative["labels"], k)
        for model, result in (("direct", direct), ("cumulative", cumulative)):
            fixed = list(result["labels"])
            replay = solve_exact(graph, k, model, 30, fixed_labels=fixed)
            assert replay["objective"] == result["objective"]
            fixed[0] = None
            neighborhood = solve_exact(graph, k, model, 30, fixed_labels=fixed)
            assert neighborhood["objective"] <= result["objective"]
            assert is_feasible(graph, neighborhood["labels"], k)
        loose = [2] * graph.n
        improved = solve_exact(
            graph, k, "direct", 30,
            fixed_labels=[None] * graph.n,
            objective_upper_bound=sum(loose) - 1,
        )
        assert improved["objective"] < sum(loose)
    print("exact model tests passed")


if __name__ == "__main__":
    main()

