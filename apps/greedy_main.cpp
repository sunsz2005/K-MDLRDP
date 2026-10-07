#include "cli_util.hpp"
#include "kmdrdp/graph.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: kmdrdp_greedy INSTANCE_DIR K [SEED=13579] [CONSTRUCTIONS=10]\n";
        return 2;
    }
    try {
        const std::string directory = argv[1];
        const int k = std::stoi(argv[2]);
        const std::uint64_t seed = argc >= 4 ? std::stoull(argv[3]) : 13579;
        const int constructions = argc >= 5 ? std::stoi(argv[4]) : 10;
        const auto graph = kmdrdp::MultiplexGraph::load_directory(directory);
        const auto result = kmdrdp::solve_greedy(graph, k, seed, constructions);
        kmdrdp::cli::print_result("GREEDY", directory, k, seed, result);
        return result.feasible ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

