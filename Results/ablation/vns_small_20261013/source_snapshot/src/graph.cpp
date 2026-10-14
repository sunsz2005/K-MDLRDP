#include "kmdrdp/graph.hpp"

#include <algorithm>
#include <dirent.h>
#include <fstream>
#include <regex>
#include <stdexcept>

namespace kmdrdp {
namespace {

int layer_number(const std::string& filename) {
    static const std::regex pattern(R"(^layer_([0-9]+)\.txt$)", std::regex::icase);
    std::smatch match;
    if (!std::regex_match(filename, match, pattern)) {
        return -1;
    }
    return std::stoi(match[1].str());
}

}  // namespace

MultiplexGraph MultiplexGraph::load_directory(const std::string& directory) {
    DIR* handle = opendir(directory.c_str());
    if (!handle) {
        throw std::runtime_error("instance directory does not exist: " + directory);
    }
    std::vector<std::string> files;
    while (const dirent* entry = readdir(handle)) {
        const std::string filename = entry->d_name;
        if (layer_number(filename) > 0) {
            files.push_back(filename);
        }
    }
    closedir(handle);
    std::sort(files.begin(), files.end(), [](const auto& left, const auto& right) {
        return layer_number(left) < layer_number(right);
    });
    if (files.empty()) {
        throw std::runtime_error("no layer_*.txt files found in " + directory);
    }

    MultiplexGraph graph;
    const auto slash = directory.find_last_of("/\\");
    graph.name = slash == std::string::npos ? directory : directory.substr(slash + 1);
    for (const auto& filename : files) {
        const std::string file = directory + "/" + filename;
        std::ifstream input(file);
        if (!input) {
            throw std::runtime_error("cannot open " + file);
        }
        int n = 0;
        long long declared_edges = 0;
        if (!(input >> n >> declared_edges) || n <= 0 || declared_edges < 0) {
            throw std::runtime_error("invalid header in " + file);
        }
        if (graph.n == 0) {
            graph.n = n;
        } else if (graph.n != n) {
            throw std::runtime_error("layers do not share the same number of vertices");
        }

        std::vector<std::vector<int>> adjacency(n);
        int u = 0;
        int v = 0;
        long long read_edges = 0;
        while (input >> u >> v) {
            --u;
            --v;
            if (u < 0 || v < 0 || u >= n || v >= n) {
                throw std::runtime_error("edge endpoint outside 1..n in " + file);
            }
            if (u == v) {
                continue;
            }
            adjacency[u].push_back(v);
            adjacency[v].push_back(u);
            ++read_edges;
        }
        for (auto& neighbors : adjacency) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }
        if (read_edges != declared_edges) {
            throw std::runtime_error("edge count does not match header in " + file);
        }
        graph.adjacency.push_back(std::move(adjacency));
    }
    return graph;
}

}  // namespace kmdrdp

