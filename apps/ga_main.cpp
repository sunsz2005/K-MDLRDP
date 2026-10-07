#include "cli_util.hpp"
#include "kmdrdp/graph.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 12) {
        std::cerr << "Usage: kmdrdp_ga INSTANCE_DIR K POP GENERATIONS ELITES MUTATION_PCT "
                     "SELECTION CROSSOVER_RATE SEED LOCAL_SEARCH TIME_LIMIT "
                     "[--init greedy|random] [--initial-solution FILE]\n";
        return 2;
    }
    try {
        const std::string directory = argv[1];
        const int k = std::stoi(argv[2]);
        kmdrdp::GAConfig config;
        config.population_size = std::stoi(argv[3]);
        config.generations = std::stoi(argv[4]);
        config.elites = std::stoi(argv[5]);
        config.mutation_percentage = std::stod(argv[6]);
        const std::string selection = argv[7];
        if (selection != "TOURNAMENT" && selection != "tournament") {
            throw std::invalid_argument("only TOURNAMENT selection is supported");
        }
        config.crossover_rate = std::stod(argv[8]);
        config.seed = std::stoull(argv[9]);
        config.local_search = kmdrdp::cli::local_search_type(argv[10]);
        config.time_limit_seconds = std::stod(argv[11]);
        config.greedy_initialization =
            kmdrdp::cli::parse_initialization(argc, argv, 12, true);
        config.initial_labels =
            kmdrdp::cli::parse_initial_solution(argc, argv, 12);

        const auto graph = kmdrdp::MultiplexGraph::load_directory(directory);
        const auto result = kmdrdp::solve_ga(graph, k, config);
        kmdrdp::cli::print_result(config.greedy_initialization ? "GA-G" : "GA-R",
                                  directory, k, config.seed, result);
        return result.feasible ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

