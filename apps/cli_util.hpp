#pragma once

#include "kmdrdp/algorithms.hpp"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace kmdrdp::cli {

inline LocalSearchType local_search_type(const std::string& value) {
    const int parsed = std::stoi(value);
    if (parsed == 1) {
        return LocalSearchType::FirstImprovement;
    }
    if (parsed == 2) {
        return LocalSearchType::BestImprovement;
    }
    throw std::invalid_argument("local search must be 1 (first) or 2 (best)");
}

inline bool parse_initialization(int argc, char** argv, int start, bool default_value) {
    bool greedy = default_value;
    for (int i = start; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == "--init") {
            const std::string value = argv[i + 1];
            if (value == "greedy") {
                greedy = true;
            } else if (value == "random") {
                greedy = false;
            } else {
                throw std::invalid_argument("--init must be greedy or random");
            }
        }
    }
    return greedy;
}

inline DominationMode parse_domination_mode(int argc, char** argv, int start) {
    for (int i = start; i + 1 < argc; ++i) {
        if (std::string(argv[i]) != "--mode") {
            continue;
        }
        const std::string value = argv[i + 1];
        if (value == "roman") {
            return DominationMode::Roman;
        }
        if (value == "double") {
            return DominationMode::DoubleRoman;
        }
        throw std::invalid_argument("--mode must be roman or double");
    }
    return DominationMode::DoubleRoman;
}

inline std::vector<std::uint8_t> parse_initial_solution(int argc, char** argv, int start) {
    for (int i = start; i + 1 < argc; ++i) {
        if (std::string(argv[i]) != "--initial-solution") {
            continue;
        }
        std::ifstream input(argv[i + 1]);
        if (!input) {
            throw std::runtime_error("cannot open --initial-solution file");
        }
        std::vector<std::uint8_t> labels;
        int value = 0;
        while (input >> value) {
            if (value < 0 || value > 3) {
                throw std::invalid_argument("initial labels must be in {0,1,2,3}");
            }
            labels.push_back(static_cast<std::uint8_t>(value));
        }
        return labels;
    }
    return {};
}

inline void print_result(const std::string& algorithm,
                         const std::string& directory,
                         int k,
                         std::uint64_t seed,
                         const SolveResult& result) {
    std::cout << "==================================================\n"
              << "ALGORITHM: " << algorithm << '\n'
              << "DIRECTORY: " << directory << '\n'
              << "K: " << k << '\n'
              << "SEED: " << seed << '\n'
              << "FEASIBLE: " << (result.feasible ? "YES" : "NO") << '\n'
              << "VALID_LAYERS: " << result.valid_layers << '\n'
              << "ITERATIONS: " << result.iterations << '\n'
              << "OBJECTIVE: " << result.objective << '\n'
              << std::fixed << std::setprecision(6)
              << "ELAPSED_SECONDS: " << result.elapsed_seconds << '\n'
              << "TIME_TO_BEST_SECONDS: " << result.time_to_best_seconds << '\n'
              << "SOLUTION:\n  " << labels_to_string(result.labels) << '\n'
              << "LABEL_1_VERTICES: " << indices_with_label(result.labels, 1) << '\n'
              << "LABEL_2_VERTICES: " << indices_with_label(result.labels, 2) << '\n'
              << "LABEL_3_VERTICES: " << indices_with_label(result.labels, 3) << '\n'
              << "==================================================\n";
}

}  // namespace kmdrdp::cli

