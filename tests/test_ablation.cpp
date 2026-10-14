#include "kmdrdp/algorithms.hpp"
#include <cassert>
#include <iostream>

int main() {
    kmdrdp::MultiplexGraph graph;
    graph.n = 7;
    graph.name = "ablation-test";
    graph.adjacency = {{{1}, {0,2,3}, {1,4}, {1,4}, {2,3,5}, {4,6}, {5}},
                       {{2}, {2,3}, {0,1,4}, {1,5}, {2,6}, {3,6}, {4,5}}};
    for (bool deficit : {false, true}) {
        for (int trial = 0; trial < 40; ++trial) {
            std::mt19937_64 rng(17 + trial);
            std::vector<std::uint8_t> labels(graph.n);
            for (auto& value : labels) value = rng() % 4;
            kmdrdp::EvaluationState state(graph, labels, 2);
            kmdrdp::SearchDiagnostics diagnostic;
            kmdrdp::repair(state, rng, deficit, &diagnostic);
            assert(kmdrdp::full_is_feasible(graph, state.labels(), 2));
            assert(diagnostic.repair_calls == 1);
            assert(diagnostic.repairs_without_fallback + diagnostic.fallback_count == 1);
        }
    }
    const std::vector<std::uint8_t> initial(graph.n, 2);
    for (bool deficit : {false, true}) for (bool adaptive : {false, true}) {
        kmdrdp::VNSConfig config;
        config.initial_labels = initial;
        config.max_iterations = 30;
        config.max_no_improve = 1000;
        config.time_limit_seconds = 10;
        config.seed = 654;
        config.deficit_guided_repair = deficit;
        config.adaptive_tie_probability = adaptive;
        kmdrdp::SearchDiagnostics diagnostic;
        config.diagnostics = &diagnostic;
        const auto result = kmdrdp::solve_vns(graph, 2, config);
        assert(kmdrdp::full_is_feasible(graph, result.labels, 2));
        assert(result.objective <= 2 * graph.n);
        assert(diagnostic.shaking_attempts >= 30);
        if (!adaptive) assert(diagnostic.final_tie_probability == 0.4);
    }
    kmdrdp::VNSConfig config;
    config.initial_labels = initial;
    config.max_iterations = 30;
    config.max_no_improve = 1000;
    config.time_limit_seconds = 10;
    config.seed = 654;
    const auto uninstrumented = kmdrdp::solve_vns(graph, 2, config);
    kmdrdp::SearchDiagnostics diagnostic;
    config.diagnostics = &diagnostic;
    const auto instrumented = kmdrdp::solve_vns(graph, 2, config);
    assert(uninstrumented.labels == instrumented.labels);
    assert(uninstrumented.iterations == instrumented.iterations);
    std::cout << "ablation tests passed: 80 repair replays, 4 variants, instrumentation equivalence\n";
}
