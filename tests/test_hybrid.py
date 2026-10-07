from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from exact_ilp import is_feasible, load_instance  # noqa: E402
from ga_ilp_hybrid import hybrid_refine  # noqa: E402


def main():
    graph = load_instance(ROOT / "Instances" / "sample-2layer")
    incumbent = [2] * graph.n
    assert is_feasible(graph, incumbent, 1)
    refined, history = hybrid_refine(
        graph=graph,
        k=1,
        incumbent=incumbent,
        rng=random.Random(13579),
        deadline=time.monotonic() + 10,
        subproblem_time=5,
        size_fractions=[1.0],
        max_subproblems=1,
    )
    assert history and history[0]["accepted"]
    assert sum(refined) < sum(incumbent)
    assert is_feasible(graph, refined, 1)
    print("hybrid refinement test passed")


if __name__ == "__main__":
    main()
