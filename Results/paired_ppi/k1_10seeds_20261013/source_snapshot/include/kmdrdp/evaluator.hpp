#pragma once

#include "kmdrdp/graph.hpp"

#include <cstdint>
#include <tuple>
#include <vector>

namespace kmdrdp {

enum class DominationMode { DoubleRoman, Roman };

struct EvalKey {
    int violations = 0;
    int weight = 0;
    int deficit = 0;

    auto tuple() const { return std::tie(violations, weight, deficit); }
    bool operator<(const EvalKey& other) const { return tuple() < other.tuple(); }
    bool operator==(const EvalKey& other) const { return tuple() == other.tuple(); }
};

struct RepairKey {
    int violations = 0;
    int deficit = 0;
    int weight = 0;

    auto tuple() const { return std::tie(violations, deficit, weight); }
    bool operator<(const RepairKey& other) const { return tuple() < other.tuple(); }
};

class EvaluationState {
public:
    EvaluationState(const MultiplexGraph& graph,
                    std::vector<std::uint8_t> labels,
                    int k,
                    DominationMode mode = DominationMode::DoubleRoman);

    const MultiplexGraph& graph() const { return graph_; }
    const std::vector<std::uint8_t>& labels() const { return labels_; }
    int k() const { return k_; }
    DominationMode mode() const { return mode_; }
    std::uint8_t max_label() const {
        return mode_ == DominationMode::Roman ? 2 : 3;
    }
    int weight() const { return weight_; }
    bool feasible() const;
    int valid_layer_count() const;
    EvalKey key() const;
    RepairKey repair_key() const;
    double fitness() const;

    void apply_move(int vertex, std::uint8_t new_label);
    EvalKey key_after_move(int vertex, std::uint8_t new_label);
    RepairKey repair_key_after_move(int vertex, std::uint8_t new_label);

    std::vector<int> valid_layers() const;
    std::vector<int> layer_violations() const { return layer_violations_; }
    std::vector<int> layer_deficits() const { return layer_deficits_; }

private:
    const MultiplexGraph& graph_;
    int k_;
    DominationMode mode_;
    std::vector<std::uint8_t> labels_;
    std::vector<std::vector<int>> count2_;
    std::vector<std::vector<int>> count3_;
    std::vector<std::vector<std::uint8_t>> violating_;
    std::vector<std::vector<int>> deficit_;
    std::vector<int> layer_violations_;
    std::vector<int> layer_deficits_;
    int weight_ = 0;

    bool is_violating(std::uint8_t label, int count2, int count3) const;
    int protection_deficit(std::uint8_t label, int count2, int count3) const;
    std::pair<int, int> best_k_layer_sums() const;
};

bool full_is_feasible(const MultiplexGraph& graph,
                      const std::vector<std::uint8_t>& labels,
                      int k,
                      DominationMode mode = DominationMode::DoubleRoman);

}  // namespace kmdrdp

