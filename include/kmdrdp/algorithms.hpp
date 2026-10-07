#pragma once

#include "kmdrdp/evaluator.hpp"

#include <chrono>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace kmdrdp {

enum class LocalSearchType { FirstImprovement = 1, BestImprovement = 2 };

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
};

std::vector<std::uint8_t> greedy_construct(const MultiplexGraph& graph,
                                           int k,
                                           std::mt19937_64& rng,
                                           DominationMode mode = DominationMode::DoubleRoman);
void repair(EvaluationState& state, std::mt19937_64& rng);
void reduce_feasible(EvaluationState& state, std::mt19937_64& rng);
void local_search(EvaluationState& state,
                  LocalSearchType type,
                  std::mt19937_64& rng,
                  const std::chrono::steady_clock::time_point& deadline);

SolveResult solve_greedy(const MultiplexGraph& graph,
                         int k,
                         std::uint64_t seed,
                         int constructions);
SolveResult solve_ga(const MultiplexGraph& graph, int k, const GAConfig& config);
SolveResult solve_vns(const MultiplexGraph& graph, int k, const VNSConfig& config);

std::string labels_to_string(const std::vector<std::uint8_t>& labels);
std::string indices_with_label(const std::vector<std::uint8_t>& labels, std::uint8_t label);

}  // namespace kmdrdp

