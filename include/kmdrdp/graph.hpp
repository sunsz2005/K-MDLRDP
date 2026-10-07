#pragma once

#include <string>
#include <vector>

namespace kmdrdp {

struct MultiplexGraph {
    int n = 0;
    std::vector<std::vector<std::vector<int>>> adjacency;
    std::string name;

    int layers() const { return static_cast<int>(adjacency.size()); }
    static MultiplexGraph load_directory(const std::string& directory);
};

}  // namespace kmdrdp

