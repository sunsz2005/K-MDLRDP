#include "kmdrdp/algorithms.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>

namespace kmdrdp {
namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

bool before_deadline(const Clock::time_point& deadline) {
    return Clock::now() < deadline;
}

std::vector<std::uint8_t> random_labels(int n,
                                        std::mt19937_64& rng,
                                        DominationMode mode) {
    const int max_label = mode == DominationMode::Roman ? 2 : 3;
    std::uniform_int_distribution<int> label_distribution(0, max_label);
    std::vector<std::uint8_t> labels(n);
    for (auto& label : labels) {
        label = static_cast<std::uint8_t>(label_distribution(rng));
    }
    return labels;
}

std::vector<std::uint8_t> best_greedy(const MultiplexGraph& graph,
                                      int k,
                                      int candidates,
                                      std::mt19937_64& rng,
                                      DominationMode mode = DominationMode::DoubleRoman) {
    candidates = std::max(1, candidates);
    std::vector<std::uint8_t> best;
    int best_weight = std::numeric_limits<int>::max();
    for (int i = 0; i < candidates; ++i) {
        auto labels = greedy_construct(graph, k, rng, mode);
        const int weight = std::accumulate(labels.begin(), labels.end(), 0);
        if (weight < best_weight) {
            best_weight = weight;
            best = std::move(labels);
        }
    }
    return best;
}

struct Individual {
    std::vector<std::uint8_t> labels;
    EvalKey key;
};

Individual improve_individual(const MultiplexGraph& graph,
                              int k,
                              std::vector<std::uint8_t> labels,
                              LocalSearchType type,
                              std::mt19937_64& rng,
                              const Clock::time_point& deadline) {
    EvaluationState state(graph, std::move(labels), k);
    if (!state.feasible()) {
        repair(state, rng);
    }
    if (before_deadline(deadline)) {
        local_search(state, type, rng, deadline);
    }
    return {state.labels(), state.key()};
}

bool individual_less(const Individual& left, const Individual& right) {
    return left.key < right.key;
}

}  // namespace

std::vector<std::uint8_t> greedy_construct(const MultiplexGraph& graph,
                                           int k,
                                           std::mt19937_64& rng,
                                           DominationMode mode) {
    EvaluationState state(graph, std::vector<std::uint8_t>(graph.n, 0), k, mode);
    while (!state.feasible()) {
        bool found = false;
        RepairKey best_key;
        std::vector<std::pair<int, std::uint8_t>> best_moves;
        for (int vertex = 0; vertex < graph.n; ++vertex) {
            const auto old_label = state.labels()[vertex];
            for (std::uint8_t label = old_label + 1; label <= state.max_label(); ++label) {
                const auto candidate = state.repair_key_after_move(vertex, label);
                if (!found || candidate < best_key) {
                    found = true;
                    best_key = candidate;
                    best_moves.assign(1, {vertex, label});
                } else if (candidate.tuple() == best_key.tuple()) {
                    best_moves.push_back({vertex, label});
                }
            }
        }
        if (!found) {
            break;
        }
        std::uniform_int_distribution<std::size_t> choice(0, best_moves.size() - 1);
        const auto move = best_moves[choice(rng)];
        state.apply_move(move.first, move.second);
    }
    if (!state.feasible()) {
        for (int vertex = 0; vertex < graph.n; ++vertex) {
            state.apply_move(vertex, 2);
        }
    }
    reduce_feasible(state, rng);
    return state.labels();
}

void repair(EvaluationState& state, std::mt19937_64& rng,
            bool deficit_guided, SearchDiagnostics* diagnostics) {
    const auto repair_start = diagnostics ? Clock::now() : Clock::time_point{};
    if (diagnostics) ++diagnostics->repair_calls;
    const auto ranking_key = [deficit_guided](RepairKey key) {
        if (!deficit_guided) key.deficit = 0;
        return key;
    };
    const int max_steps = std::max(20, 8 * state.graph().n);
    for (int step = 0; step < max_steps && !state.feasible(); ++step) {
        const RepairKey current = ranking_key(state.repair_key());
        bool found = false;
        RepairKey best_key;
        std::vector<std::pair<int, std::uint8_t>> best_moves;
        for (int vertex = 0; vertex < state.graph().n; ++vertex) {
            const auto old_label = state.labels()[vertex];
            for (std::uint8_t label = 0; label <= state.max_label(); ++label) {
                if (label == old_label) {
                    continue;
                }
                const auto candidate = ranking_key(state.repair_key_after_move(vertex, label));
                if (diagnostics) ++diagnostics->repair_evaluations;
                if (!(candidate < current)) {
                    continue;
                }
                if (!found || candidate < best_key) {
                    found = true;
                    best_key = candidate;
                    best_moves.assign(1, {vertex, label});
                } else if (candidate.tuple() == best_key.tuple()) {
                    best_moves.push_back({vertex, label});
                }
            }
        }
        if (!found) {
            break;
        }
        std::uniform_int_distribution<std::size_t> choice(0, best_moves.size() - 1);
        const auto move = best_moves[choice(rng)];
        state.apply_move(move.first, move.second);
        if (diagnostics) ++diagnostics->repair_steps;
    }

    if (!state.feasible()) {
        if (diagnostics) ++diagnostics->fallback_count;
        auto fallback = greedy_construct(state.graph(), state.k(), rng, state.mode());
        for (int vertex = 0; vertex < state.graph().n; ++vertex) {
            state.apply_move(vertex, fallback[vertex]);
        }
    } else if (diagnostics) {
        ++diagnostics->repairs_without_fallback;
    }
    if (diagnostics) diagnostics->repair_seconds += seconds_since(repair_start);
}

void reduce_feasible(EvaluationState& state, std::mt19937_64& rng) {
    if (!state.feasible()) {
        return;
    }
    std::vector<int> vertices(state.graph().n);
    std::iota(vertices.begin(), vertices.end(), 0);
    bool improved = true;
    while (improved) {
        improved = false;
        std::shuffle(vertices.begin(), vertices.end(), rng);
        for (int vertex : vertices) {
            const auto old_label = state.labels()[vertex];
            for (std::uint8_t label = 0; label < old_label; ++label) {
                const auto candidate = state.key_after_move(vertex, label);
                if (candidate.violations == 0 && candidate.weight < state.weight()) {
                    state.apply_move(vertex, label);
                    improved = true;
                    break;
                }
            }
        }
    }
}

void local_search(EvaluationState& state,
                  LocalSearchType type,
                  std::mt19937_64& rng,
                  const Clock::time_point& deadline,
                  SearchDiagnostics* diagnostics) {
    if (!state.feasible()) {
        repair(state, rng);
    }
    std::vector<int> vertices(state.graph().n);
    std::iota(vertices.begin(), vertices.end(), 0);

    while (before_deadline(deadline)) {
        const EvalKey current = state.key();
        bool found = false;
        int best_vertex = -1;
        std::uint8_t best_label = 0;
        EvalKey best_key = current;
        std::shuffle(vertices.begin(), vertices.end(), rng);

        for (int vertex : vertices) {
            std::vector<std::uint8_t> labels;
            for (std::uint8_t label = 0; label <= state.max_label(); ++label) {
                labels.push_back(label);
            }
            std::shuffle(labels.begin(), labels.end(), rng);
            for (const auto label : labels) {
                if (label == state.labels()[vertex]) {
                    continue;
                }
                const auto candidate = state.key_after_move(vertex, label);
                if (diagnostics) ++diagnostics->local_search_evaluations;
                if (candidate < best_key) {
                    found = true;
                    best_vertex = vertex;
                    best_label = label;
                    best_key = candidate;
                    if (type == LocalSearchType::FirstImprovement) {
                        break;
                    }
                }
            }
            if (found && type == LocalSearchType::FirstImprovement) {
                break;
            }
            if (!before_deadline(deadline)) {
                return;
            }
        }
        if (!found) {
            return;
        }
        state.apply_move(best_vertex, best_label);
    }
}

SolveResult solve_greedy(const MultiplexGraph& graph,
                         int k,
                         std::uint64_t seed,
                         int constructions) {
    const auto start = Clock::now();
    std::mt19937_64 rng(seed);
    auto labels = best_greedy(graph, k, constructions, rng);
    EvaluationState state(graph, labels, k);
    const double elapsed = seconds_since(start);
    return {labels, state.weight(), state.valid_layer_count(), constructions,
            elapsed, elapsed, state.feasible()};
}

SolveResult solve_ga(const MultiplexGraph& graph, int k, const GAConfig& config) {
    const auto start = Clock::now();
    const auto deadline = start + std::chrono::duration_cast<Clock::duration>(
                                      std::chrono::duration<double>(config.time_limit_seconds));
    std::mt19937_64 rng(config.seed);
    std::vector<Individual> population;
    population.reserve(config.population_size);

    if (!config.initial_labels.empty()) {
        if (static_cast<int>(config.initial_labels.size()) != graph.n ||
            std::any_of(config.initial_labels.begin(), config.initial_labels.end(),
                        [](std::uint8_t label) { return label > 3; })) {
            throw std::invalid_argument("invalid injected GA solution");
        }
        population.push_back(improve_individual(graph, k, config.initial_labels,
                                                config.local_search, rng, deadline));
    }

    for (int i = static_cast<int>(population.size());
         i < config.population_size && before_deadline(deadline); ++i) {
        std::vector<std::uint8_t> labels;
        if (config.greedy_initialization && i < config.greedy_candidates) {
            labels = greedy_construct(graph, k, rng);
        } else {
            labels = random_labels(graph.n, rng, DominationMode::DoubleRoman);
        }
        population.push_back(improve_individual(graph, k, std::move(labels),
                                                config.local_search, rng, deadline));
    }
    if (population.empty()) {
        population.push_back(improve_individual(graph, k, greedy_construct(graph, k, rng),
                                                config.local_search, rng, deadline));
    }
    std::sort(population.begin(), population.end(), individual_less);
    Individual best = population.front();
    double time_to_best = seconds_since(start);
    int completed_generations = 0;

    auto tournament = [&]() -> const Individual& {
        std::uniform_int_distribution<std::size_t> pick(0, population.size() - 1);
        std::size_t winner = pick(rng);
        for (int t = 1; t < std::max(1, config.tournament_size); ++t) {
            const std::size_t candidate = pick(rng);
            if (population[candidate].key < population[winner].key) {
                winner = candidate;
            }
        }
        return population[winner];
    };

    for (int generation = 0;
         generation < config.generations && before_deadline(deadline);
         ++generation) {
        std::sort(population.begin(), population.end(), individual_less);
        std::vector<Individual> next;
        const int elite_count = std::min({config.elites, config.population_size,
                                          static_cast<int>(population.size())});
        next.insert(next.end(), population.begin(), population.begin() + elite_count);

        while (static_cast<int>(next.size()) < config.population_size &&
               before_deadline(deadline)) {
            auto child1 = tournament().labels;
            auto child2 = tournament().labels;
            std::uniform_real_distribution<double> probability(0.0, 1.0);
            if (graph.n > 1 && probability(rng) < config.crossover_rate) {
                std::uniform_int_distribution<int> point_distribution(1, graph.n - 1);
                const int point = point_distribution(rng);
                for (int i = point; i < graph.n; ++i) {
                    std::swap(child1[i], child2[i]);
                }
            }

            auto mutate = [&](std::vector<std::uint8_t>& child) {
                int count = std::max(1, static_cast<int>(std::ceil(
                    graph.n * config.mutation_percentage / 100.0)));
                count = std::min(count, graph.n);
                std::vector<int> positions(graph.n);
                std::iota(positions.begin(), positions.end(), 0);
                std::shuffle(positions.begin(), positions.end(), rng);
                std::uniform_int_distribution<int> offset(1, 3);
                for (int j = 0; j < count; ++j) {
                    const int v = positions[j];
                    child[v] = static_cast<std::uint8_t>((child[v] + offset(rng)) % 4);
                }
            };
            mutate(child1);
            mutate(child2);
            next.push_back(improve_individual(graph, k, std::move(child1),
                                              config.local_search, rng, deadline));
            if (static_cast<int>(next.size()) < config.population_size &&
                before_deadline(deadline)) {
                next.push_back(improve_individual(graph, k, std::move(child2),
                                                  config.local_search, rng, deadline));
            }
        }
        population = std::move(next);
        if (population.empty()) {
            break;
        }
        std::sort(population.begin(), population.end(), individual_less);
        if (population.front().key < best.key) {
            best = population.front();
            time_to_best = seconds_since(start);
        }
        completed_generations = generation + 1;
    }

    EvaluationState final_state(graph, best.labels, k);
    const double elapsed = seconds_since(start);
    return {best.labels, final_state.weight(), final_state.valid_layer_count(),
            completed_generations, elapsed, time_to_best, final_state.feasible()};
}

SolveResult solve_vns(const MultiplexGraph& graph, int k, const VNSConfig& config) {
    const auto start = Clock::now();
    const auto deadline = start + std::chrono::duration_cast<Clock::duration>(
                                      std::chrono::duration<double>(config.time_limit_seconds));
    std::mt19937_64 rng(config.seed);

    std::vector<std::uint8_t> initial;
    if (!config.initial_labels.empty()) {
        initial = config.initial_labels;
    } else {
        initial = config.greedy_initialization
            ? best_greedy(graph, k, config.greedy_candidates, rng, config.mode)
            : random_labels(graph.n, rng, config.mode);
    }
    EvaluationState current(graph, std::move(initial), k, config.mode);
    if (!current.feasible()) {
        repair(current, rng, config.deficit_guided_repair, config.diagnostics);
    }
    local_search(current, config.local_search, rng, deadline, config.diagnostics);
    auto best_labels = current.labels();
    EvalKey best_key = current.key();
    double time_to_best = seconds_since(start);

    int no_improve = 0;
    int iterations = 0;
    double tie_probability = std::clamp(config.tie_probability, 0.0, 1.0);
    std::uniform_real_distribution<double> probability(0.0, 1.0);

    while (iterations < config.max_iterations && no_improve < config.max_no_improve &&
           before_deadline(deadline)) {
        bool accepted = false;
        for (int kappa = std::max(1, config.kappa_min);
             kappa <= std::min(config.kappa_max, graph.n) && before_deadline(deadline);
             ++kappa) {
            if (config.diagnostics) ++config.diagnostics->shaking_attempts;
            auto candidate_labels = current.labels();
            std::vector<int> positions(graph.n);
            std::iota(positions.begin(), positions.end(), 0);
            std::shuffle(positions.begin(), positions.end(), rng);
            const int max_label = config.mode == DominationMode::Roman ? 2 : 3;
            std::uniform_int_distribution<int> offset(1, max_label);
            for (int j = 0; j < kappa; ++j) {
                const int v = positions[j];
                candidate_labels[v] = static_cast<std::uint8_t>(
                    (candidate_labels[v] + offset(rng)) % (max_label + 1));
            }

            EvaluationState candidate(graph, std::move(candidate_labels), k, config.mode);
            if (!candidate.feasible()) {
                repair(candidate, rng, config.deficit_guided_repair, config.diagnostics);
            }
            local_search(candidate, config.local_search, rng, deadline, config.diagnostics);
            const EvalKey candidate_key = candidate.key();
            const EvalKey current_key = current.key();
            if (candidate_key < current_key ||
                (candidate_key == current_key && probability(rng) < tie_probability)) {
                if (config.diagnostics) {
                    if (candidate_key == current_key) ++config.diagnostics->accepted_ties;
                    else ++config.diagnostics->accepted_improvements;
                }
                for (int v = 0; v < graph.n; ++v) {
                    current.apply_move(v, candidate.labels()[v]);
                }
                accepted = true;
                if (candidate_key < best_key) {
                    best_key = candidate_key;
                    best_labels = candidate.labels();
                    time_to_best = seconds_since(start);
                    no_improve = 0;
                    if (config.adaptive_tie_probability)
                        tie_probability = std::max(0.2, tie_probability * 0.9);
                } else {
                    ++no_improve;
                    if (config.adaptive_tie_probability)
                        tie_probability = std::min(0.95, tie_probability + 0.01);
                }
                break;
            }
            ++no_improve;
        }
        if (!accepted && config.adaptive_tie_probability) {
            tie_probability = std::min(0.95, tie_probability + 0.01);
        }
        ++iterations;
    }

    EvaluationState final_state(graph, best_labels, k, config.mode);
    if (config.diagnostics) config.diagnostics->final_tie_probability = tie_probability;
    const double elapsed = seconds_since(start);
    return {best_labels, final_state.weight(), final_state.valid_layer_count(),
            iterations, elapsed, time_to_best, final_state.feasible()};
}

std::string labels_to_string(const std::vector<std::uint8_t>& labels) {
    std::ostringstream output;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (i) {
            output << ' ';
        }
        output << static_cast<int>(labels[i]);
    }
    return output.str();
}

std::string indices_with_label(const std::vector<std::uint8_t>& labels, std::uint8_t label) {
    std::ostringstream output;
    bool first = true;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (labels[i] == label) {
            if (!first) {
                output << ' ';
            }
            output << i + 1;
            first = false;
        }
    }
    return output.str();
}

}  // namespace kmdrdp

