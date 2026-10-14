#include "cli_util.hpp"
#include "kmdrdp/graph.hpp"
#include <exception>
#include <iostream>
#include <numeric>

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: vns_ablation INSTANCE K SEED SECONDS "
                     "full|no-repair-deficit|fixed-tie INITIAL_LABEL_FILE\n";
        return 2;
    }
    try {
        const auto graph = kmdrdp::MultiplexGraph::load_directory(argv[1]);
        const int k = std::stoi(argv[2]);
        const std::string variant = argv[5];
        kmdrdp::VNSConfig config;
        config.seed = std::stoull(argv[3]);
        config.time_limit_seconds = std::stod(argv[4]);
        config.deficit_guided_repair = variant != "no-repair-deficit";
        config.adaptive_tie_probability = variant != "fixed-tie";
        if (variant != "full" && variant != "no-repair-deficit" && variant != "fixed-tie")
            throw std::invalid_argument("unknown ablation variant");
        std::ifstream initial(argv[6]);
        int value;
        while (initial >> value) {
            if (value < 0 || value > 3) throw std::invalid_argument("invalid initial label");
            config.initial_labels.push_back(static_cast<std::uint8_t>(value));
        }
        if (!initial.eof() || !kmdrdp::full_is_feasible(graph, config.initial_labels, k))
            throw std::invalid_argument("initial file missing, malformed or infeasible");
        kmdrdp::SearchDiagnostics diagnostics;
        config.diagnostics = &diagnostics;
        const auto result = kmdrdp::solve_vns(graph, k, config);
        if (!kmdrdp::full_is_feasible(graph, result.labels, k))
            throw std::runtime_error("independent C++ feasibility replay failed");
        std::cout << std::setprecision(12) << std::boolalpha
                  << "{\"variant\":\"" << variant << "\",\"seed\":" << config.seed
                  << ",\"n\":" << graph.n << ",\"layers\":" << graph.layers()
                  << ",\"k\":" << k << ",\"objective\":" << result.objective
                  << ",\"initial_objective\":" << std::accumulate(config.initial_labels.begin(), config.initial_labels.end(), 0)
                  << ",\"feasible\":" << result.feasible
                  << ",\"valid_layers\":" << result.valid_layers
                  << ",\"iterations\":" << result.iterations
                  << ",\"elapsed_seconds\":" << result.elapsed_seconds
                  << ",\"time_to_best_seconds\":" << result.time_to_best_seconds
                  << ",\"repair_calls\":" << diagnostics.repair_calls
                  << ",\"repair_steps\":" << diagnostics.repair_steps
                  << ",\"repair_evaluations\":" << diagnostics.repair_evaluations
                  << ",\"repairs_without_fallback\":" << diagnostics.repairs_without_fallback
                  << ",\"fallback_count\":" << diagnostics.fallback_count
                  << ",\"repair_seconds\":" << diagnostics.repair_seconds
                  << ",\"local_search_evaluations\":" << diagnostics.local_search_evaluations
                  << ",\"shaking_attempts\":" << diagnostics.shaking_attempts
                  << ",\"accepted_ties\":" << diagnostics.accepted_ties
                  << ",\"accepted_improvements\":" << diagnostics.accepted_improvements
                  << ",\"final_tie_probability\":" << diagnostics.final_tie_probability
                  << ",\"labels\":[";
        for (std::size_t i = 0; i < result.labels.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << static_cast<int>(result.labels[i]);
        }
        std::cout << "]}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
