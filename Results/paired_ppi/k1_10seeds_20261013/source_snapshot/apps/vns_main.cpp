#include "cli_util.hpp"
#include "kmdrdp/graph.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 12) {
        std::cerr << "Usage: kmdrdp_vns vns INSTANCE_DIR K SEED KAPPA_MIN KAPPA_MAX "
                     "TIE_PROB LOCAL_SEARCH ITER_MAX NO_IMPROVE_MAX TIME_LIMIT "
                     "[--init greedy|random] [--mode double|roman] "
                     "[--initial-solution FILE]\n";
        return 2;
    }
    try {
        if (std::string(argv[1]) != "vns") {
            throw std::invalid_argument("first argument must be vns for interface compatibility");
        }
        const std::string directory = argv[2];
        const int k = std::stoi(argv[3]);
        kmdrdp::VNSConfig config;
        config.seed = std::stoull(argv[4]);
        config.kappa_min = std::stoi(argv[5]);
        config.kappa_max = std::stoi(argv[6]);
        config.tie_probability = std::stod(argv[7]);
        config.local_search = kmdrdp::cli::local_search_type(argv[8]);
        config.max_iterations = std::stoi(argv[9]);
        config.max_no_improve = std::stoi(argv[10]);
        config.time_limit_seconds = std::stod(argv[11]);
        config.greedy_initialization =
            kmdrdp::cli::parse_initialization(argc, argv, 12, true);
        config.mode = kmdrdp::cli::parse_domination_mode(argc, argv, 12);
        config.initial_labels = kmdrdp::cli::parse_initial_solution(argc, argv, 12);

        const auto graph = kmdrdp::MultiplexGraph::load_directory(directory);
        const auto result = kmdrdp::solve_vns(graph, k, config);
        std::string algorithm = config.greedy_initialization ? "VNS-G" : "VNS-R";
        if (config.mode == kmdrdp::DominationMode::Roman) {
            algorithm += "-Roman";
        }
        kmdrdp::cli::print_result(algorithm,
                                  directory, k, config.seed, result);
        return result.feasible ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

