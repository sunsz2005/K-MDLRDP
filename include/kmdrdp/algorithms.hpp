#pragma once

#include "kmdrdp/evaluator.hpp"

#include <chrono>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace kmdrdp {

enum class LocalSearchType { FirstImprovement = 1, BestImprovement = 2 };

struct SearchDiagnostics {
    std::uint64_t repair_calls = 0;
    std::uint64_t repair_steps = 0;
    std::uint64_t repair_evaluations = 0;
    std::uint64_t repairs_without_fallback = 0;
    std::uint64_t fallback_count = 0;
    std::uint64_t local_search_evaluations = 0;
    std::uint64_t shaking_attempts = 0;
    std::uint64_t accepted_ties = 0;
    std::uint64_t accepted_improvements = 0;
    double repair_seconds = 0.0;
    double final_tie_probability = 0.0;
};

struct SolveResult {
    std::vector<std::uint8_t> labels;
    int objective = 0;
    int valid_layers = 0;
    int iterations = 0;
    double elapsed_seconds = 0.0;
    double time_to_best_seconds = 0.0;
    bool feasible = false;
};

struct GAConfig {
    int population_size = 50;
    int generations = 100;
    int elites = 2;
    double mutation_percentage = 10.0;
    double crossover_rate = 1.0;
    int tournament_size = 5;
    LocalSearchType local_search = LocalSearchType::FirstImprovement;
    bool greedy_initialization = true;
    int greedy_candidates = 10;
    double time_limit_seconds = 3600.0;
    std::uint64_t seed = 13579;
    // Optional incumbent injected as the first population member.  This is
    // used by the epochal GA--ILP hybrid to preserve improvements.
    std::vector<std::uint8_t> initial_labels;
};

struct VNSConfig {
    int kappa_min = 2;
    int kappa_max = 12;
    double tie_probability = 0.4;
    LocalSearchType local_search = LocalSearchType::FirstImprovement;
    bool greedy_initialization = true;
    int greedy_candidates = 10;
    int max_iterations = 10000;
    int max_no_improve = 5000;
    double time_limit_seconds = 3600.0;
    std::uint64_t seed = 13579;
    DominationMode mode = DominationMode::DoubleRoman;
    std::vector<std::uint8_t> initial_labels;
    // Defaults preserve the original algorithm. Ablation changes repair
    // ranking only; layer selection, Greedy fallback and search keys stay fixed.
    bool deficit_guided_repair = true;
    bool adaptive_tie_probability = true;
    SearchDiagnostics* diagnostics = nullptr;
};

std::vector<std::uint8_t> greedy_construct(const MultiplexGraph& graph,
                                           int k,
                                           std::mt19937_64& rng,
                                           DominationMode mode = DominationMode::DoubleRoman);
void repair(EvaluationState& state, std::mt19937_64& rng,
            bool deficit_guided = true, SearchDiagnostics* diagnostics = nullptr);
void reduce_feasible(EvaluationState& state, std::mt19937_64& rng);
void local_search(EvaluationState& state,
                  LocalSearchType type,
                  std::mt19937_64& rng,
                  const std::chrono::steady_clock::time_point& deadline,
                  SearchDiagnostics* diagnostics = nullptr);

SolveResult solve_greedy(const MultiplexGraph& graph,
                         int k,
                         std::uint64_t seed,
                         int constructions);
SolveResult solve_ga(const MultiplexGraph& graph, int k, const GAConfig& config);
SolveResult solve_vns(const MultiplexGraph& graph, int k, const VNSConfig& config);

std::string labels_to_string(const std::vector<std::uint8_t>& labels);
std::string indices_with_label(const std::vector<std::uint8_t>& labels, std::uint8_t label);

}  // namespace kmdrdp

