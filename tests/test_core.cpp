#include "kmdrdp/algorithms.hpp"
#include "kmdrdp/evaluator.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

namespace {

kmdrdp::MultiplexGraph graph_from_edges(
    int n, const std::vector<std::vector<std::pair<int, int>>>& layers) {
    kmdrdp::MultiplexGraph graph;
    graph.n = n;
    graph.name = "test";
    for (const auto& edges : layers) {
        std::vector<std::vector<int>> adjacency(n);
        for (const auto [u, v] : edges) {
            adjacency[u].push_back(v);
            adjacency[v].push_back(u);
        }
        graph.adjacency.push_back(std::move(adjacency));
    }
    return graph;
}

}  // namespace

int main() {
    const auto path3 = graph_from_edges(3, {{{0, 1}, {1, 2}}});
    assert(kmdrdp::full_is_feasible(path3, {0, 3, 0}, 1));
    assert(!kmdrdp::full_is_feasible(path3, {0, 2, 0}, 1));
    assert(kmdrdp::full_is_feasible(path3, {2, 0, 2}, 1));
    assert(!kmdrdp::full_is_feasible(path3, {1, 0, 0}, 1));
    assert(kmdrdp::full_is_feasible(
        path3, {0, 2, 0}, 1, kmdrdp::DominationMode::Roman));
    assert(!kmdrdp::full_is_feasible(
        path3, {0, 1, 0}, 1, kmdrdp::DominationMode::Roman));
    assert(!kmdrdp::full_is_feasible(
        path3, {0, 3, 0}, 1, kmdrdp::DominationMode::Roman));

    const auto multiplex = graph_from_edges(
        5, {{{0, 1}, {1, 2}, {2, 3}, {3, 4}},
            {{0, 2}, {1, 2}, {2, 4}, {3, 4}}});
    std::mt19937_64 rng(20261001);
    std::uniform_int_distribution<int> label(0, 3);
    std::uniform_int_distribution<int> vertex(0, multiplex.n - 1);
    std::vector<std::uint8_t> labels(multiplex.n, 0);
    kmdrdp::EvaluationState state(multiplex, labels, 1);
    for (int step = 0; step < 1000; ++step) {
        const int v = vertex(rng);
        const auto new_label = static_cast<std::uint8_t>(label(rng));
        state.apply_move(v, new_label);
        labels[v] = new_label;
        assert(state.feasible() == kmdrdp::full_is_feasible(multiplex, labels, 1));
        assert(state.weight() == static_cast<int>(
            labels[0] + labels[1] + labels[2] + labels[3] + labels[4]));
    }

    auto greedy = kmdrdp::solve_greedy(multiplex, 2, 7, 5);
    assert(greedy.feasible);
    assert(kmdrdp::full_is_feasible(multiplex, greedy.labels, 2));

    kmdrdp::VNSConfig roman_config;
    roman_config.mode = kmdrdp::DominationMode::Roman;
    roman_config.seed = 17;
    roman_config.max_iterations = 20;
    roman_config.max_no_improve = 10;
    roman_config.time_limit_seconds = 1.0;
    auto roman = kmdrdp::solve_vns(multiplex, 2, roman_config);
    assert(roman.feasible);
    assert(kmdrdp::full_is_feasible(
        multiplex, roman.labels, 2, kmdrdp::DominationMode::Roman));
    for (const auto value : roman.labels) {
        assert(value <= 2);
        (void)value;
    }
    std::cout << "core tests passed\n";
    return 0;
}

