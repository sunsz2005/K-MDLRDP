#include "kmdrdp/evaluator.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

using Labels = std::vector<std::uint8_t>;
using Clock = std::chrono::steady_clock;

// Independent full scan: no incremental cache and no early layer stopping.
kmdrdp::EvalKey full_key(const kmdrdp::MultiplexGraph& graph, const Labels& labels, int k) {
    std::vector<std::pair<int,int>> layers;
    for (int s = 0; s < graph.layers(); ++s) {
        int violations = 0, deficit = 0;
        for (int i = 0; i < graph.n; ++i) {
            int d2 = 0, d3 = 0;
            for (int j : graph.adjacency[s][i]) {
                d2 += labels[j] == 2;
                d3 += labels[j] == 3;
            }
            if (labels[i] == 0) {
                violations += d3 == 0 && d2 < 2;
                deficit += std::max(0, 2-d2-2*d3);
            } else if (labels[i] == 1) {
                violations += d2+d3 == 0;
                deficit += d2+d3 == 0;
            }
        }
        layers.emplace_back(violations, deficit);
    }
    std::sort(layers.begin(), layers.end());
    kmdrdp::EvalKey key;
    key.weight = std::accumulate(labels.begin(), labels.end(), 0);
    for (int s = 0; s < k; ++s) {
        key.violations += layers[s].first;
        key.deficit += layers[s].second;
    }
    return key;
}

struct Move { int vertex; std::uint8_t label; bool commit; };
struct Trace { double seconds; std::vector<kmdrdp::EvalKey> keys; };

Trace execute(const kmdrdp::MultiplexGraph& graph, const Labels& initial, int k,
              const std::vector<Move>& moves, bool incremental) {
    Labels labels = initial;
    kmdrdp::EvaluationState state(graph, initial, k);
    Trace trace;
    trace.keys.reserve(moves.size());
    const auto start = Clock::now();
    for (const auto& move : moves) {
        if (incremental) {
            trace.keys.push_back(state.key_after_move(move.vertex, move.label));
            if (move.commit) state.apply_move(move.vertex, move.label);
        } else {
            const auto old = labels[move.vertex];
            labels[move.vertex] = move.label;
            trace.keys.push_back(full_key(graph, labels, k));
            if (!move.commit) labels[move.vertex] = old;
        }
    }
    trace.seconds = std::chrono::duration<double>(Clock::now()-start).count();
    return trace;
}

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: evaluation_benchmark INSTANCE K SEED PROBES REPEATS INITIAL_FILE\n";
        return 2;
    }
    try {
        const auto graph = kmdrdp::MultiplexGraph::load_directory(argv[1]);
        const int k = std::stoi(argv[2]), probes = std::stoi(argv[4]), repeats = std::stoi(argv[5]);
        if (probes < 1 || repeats < 1) throw std::invalid_argument("positive counts required");
        Labels initial;
        std::ifstream input(argv[6]);
        int value;
        while (input >> value) {
            if (value < 0 || value > 3) throw std::invalid_argument("invalid label");
            initial.push_back(static_cast<std::uint8_t>(value));
        }
        if (!kmdrdp::full_is_feasible(graph, initial, k)) throw std::runtime_error("invalid initial");
        std::mt19937_64 rng(std::stoull(argv[3]));
        Labels simulated = initial;
        std::vector<Move> moves;
        for (int i = 0; i < probes; ++i) {
            const int v = rng() % graph.n;
            const auto label = static_cast<std::uint8_t>((simulated[v]+1+rng()%3)%4);
            moves.push_back({v, label, i % 5 == 0});
            if (moves.back().commit) simulated[v] = label;
        }
        std::uint64_t checksum = 0;
        std::vector<double> inc_times, full_times;
        // Warm up both paths; exclude startup and validation from timing.
        execute(graph, initial, k, moves, true);
        execute(graph, initial, k, moves, false);
        for (int r = 0; r < repeats; ++r) {
            Trace incremental, full;
            if (r % 2 == 0) {
                incremental = execute(graph, initial, k, moves, true);
                full = execute(graph, initial, k, moves, false);
            } else {
                full = execute(graph, initial, k, moves, false);
                incremental = execute(graph, initial, k, moves, true);
            }
            for (int i = 0; i < probes; ++i) {
                if (!(incremental.keys[i] == full.keys[i]))
                    throw std::runtime_error("incremental/full key mismatch at probe " + std::to_string(i));
                checksum += static_cast<std::uint64_t>(full.keys[i].violations)*1000000
                          + full.keys[i].weight*1000 + full.keys[i].deficit;
            }
            inc_times.push_back(incremental.seconds);
            full_times.push_back(full.seconds);
        }
        std::sort(inc_times.begin(), inc_times.end());
        std::sort(full_times.begin(), full_times.end());
        const double inc = inc_times[repeats/2], full = full_times[repeats/2];
        std::cout << std::setprecision(12) << "{\"n\":" << graph.n
                  << ",\"layers\":" << graph.layers() << ",\"k\":" << k
                  << ",\"probes\":" << probes << ",\"repeats\":" << repeats
                  << ",\"verified_comparisons\":" << probes*repeats
                  << ",\"mismatches\":0,\"checksum\":" << checksum
                  << ",\"incremental_median_seconds\":" << inc
                  << ",\"full_scan_median_seconds\":" << full
                  << ",\"speedup_full_over_incremental\":" << full/inc << "}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
