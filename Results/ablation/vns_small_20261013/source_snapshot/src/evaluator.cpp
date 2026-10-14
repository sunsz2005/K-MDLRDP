#include "kmdrdp/evaluator.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace kmdrdp {

bool EvaluationState::is_violating(std::uint8_t label, int count2, int count3) const {
    if (mode_ == DominationMode::Roman) {
        return label == 0 && count2 == 0;
    }
    if (label == 0) {
        return count3 == 0 && count2 < 2;
    }
    if (label == 1) {
        return count2 + count3 == 0;
    }
    return false;
}

int EvaluationState::protection_deficit(std::uint8_t label,
                                        int count2,
                                        int count3) const {
    if (mode_ == DominationMode::Roman) {
        return label == 0 && count2 == 0 ? 1 : 0;
    }
    if (label == 0) {
        return std::max(0, 2 - count2 - 2 * count3);
    }
    if (label == 1) {
        return count2 + count3 == 0 ? 1 : 0;
    }
    return 0;
}

EvaluationState::EvaluationState(const MultiplexGraph& graph,
                                 std::vector<std::uint8_t> labels,
                                 int k,
                                 DominationMode mode)
    : graph_(graph), k_(k), mode_(mode), labels_(std::move(labels)) {
    if (k < 1 || k > graph.layers()) {
        throw std::invalid_argument("k must be between 1 and the number of layers");
    }
    if (static_cast<int>(labels_.size()) != graph.n) {
        throw std::invalid_argument("label vector has the wrong length");
    }
    if (std::any_of(labels_.begin(), labels_.end(),
                    [this](auto label) { return label > max_label(); })) {
        throw std::invalid_argument(mode_ == DominationMode::Roman
                                        ? "Roman labels must be in {0,1,2}"
                                        : "double-Roman labels must be in {0,1,2,3}");
    }

    const int m = graph.layers();
    count2_.assign(m, std::vector<int>(graph.n, 0));
    count3_.assign(m, std::vector<int>(graph.n, 0));
    violating_.assign(m, std::vector<std::uint8_t>(graph.n, 0));
    deficit_.assign(m, std::vector<int>(graph.n, 0));
    layer_violations_.assign(m, 0);
    layer_deficits_.assign(m, 0);
    weight_ = std::accumulate(labels_.begin(), labels_.end(), 0);

    for (int s = 0; s < m; ++s) {
        for (int i = 0; i < graph.n; ++i) {
            for (int neighbor : graph.adjacency[s][i]) {
                count2_[s][i] += labels_[neighbor] == 2;
                count3_[s][i] += labels_[neighbor] == 3;
            }
            violating_[s][i] = is_violating(labels_[i], count2_[s][i], count3_[s][i]);
            deficit_[s][i] = protection_deficit(labels_[i], count2_[s][i], count3_[s][i]);
            layer_violations_[s] += violating_[s][i];
            layer_deficits_[s] += deficit_[s][i];
        }
    }
}

std::pair<int, int> EvaluationState::best_k_layer_sums() const {
    std::vector<std::pair<int, int>> values;
    values.reserve(graph_.layers());
    for (int s = 0; s < graph_.layers(); ++s) {
        values.emplace_back(layer_violations_[s], layer_deficits_[s]);
    }
    std::sort(values.begin(), values.end());
    int violations = 0;
    int deficits = 0;
    for (int i = 0; i < k_; ++i) {
        violations += values[i].first;
        deficits += values[i].second;
    }
    return {violations, deficits};
}

EvalKey EvaluationState::key() const {
    const auto sums = best_k_layer_sums();
    return {sums.first, weight_, sums.second};
}

RepairKey EvaluationState::repair_key() const {
    const auto sums = best_k_layer_sums();
    return {sums.first, sums.second, weight_};
}

bool EvaluationState::feasible() const {
    return key().violations == 0;
}

int EvaluationState::valid_layer_count() const {
    return static_cast<int>(std::count(layer_violations_.begin(), layer_violations_.end(), 0));
}

double EvaluationState::fitness() const {
    const auto current = key();
    return current.violations +
           static_cast<double>(weight_) / (max_label() * graph_.n + 1.0);
}

std::vector<int> EvaluationState::valid_layers() const {
    std::vector<int> result;
    for (int s = 0; s < graph_.layers(); ++s) {
        if (layer_violations_[s] == 0) {
            result.push_back(s + 1);
        }
    }
    return result;
}

void EvaluationState::apply_move(int vertex, std::uint8_t new_label) {
    if (vertex < 0 || vertex >= graph_.n || new_label > max_label()) {
        throw std::invalid_argument("invalid move");
    }
    const std::uint8_t old_label = labels_[vertex];
    if (old_label == new_label) {
        return;
    }

    for (int s = 0; s < graph_.layers(); ++s) {
        std::vector<int> affected;
        affected.reserve(graph_.adjacency[s][vertex].size() + 1);
        affected.push_back(vertex);
        affected.insert(affected.end(), graph_.adjacency[s][vertex].begin(),
                        graph_.adjacency[s][vertex].end());

        for (int i : affected) {
            layer_violations_[s] -= violating_[s][i];
            layer_deficits_[s] -= deficit_[s][i];
        }
        for (int neighbor : graph_.adjacency[s][vertex]) {
            count2_[s][neighbor] += (new_label == 2) - (old_label == 2);
            count3_[s][neighbor] += (new_label == 3) - (old_label == 3);
        }

        labels_[vertex] = new_label;
        for (int i : affected) {
            violating_[s][i] = is_violating(labels_[i], count2_[s][i], count3_[s][i]);
            deficit_[s][i] = protection_deficit(labels_[i], count2_[s][i], count3_[s][i]);
            layer_violations_[s] += violating_[s][i];
            layer_deficits_[s] += deficit_[s][i];
        }
        labels_[vertex] = old_label;
    }
    labels_[vertex] = new_label;
    weight_ += static_cast<int>(new_label) - static_cast<int>(old_label);
}

EvalKey EvaluationState::key_after_move(int vertex, std::uint8_t new_label) {
    const std::uint8_t old_label = labels_.at(vertex);
    apply_move(vertex, new_label);
    const EvalKey result = key();
    apply_move(vertex, old_label);
    return result;
}

RepairKey EvaluationState::repair_key_after_move(int vertex, std::uint8_t new_label) {
    const std::uint8_t old_label = labels_.at(vertex);
    apply_move(vertex, new_label);
    const RepairKey result = repair_key();
    apply_move(vertex, old_label);
    return result;
}

bool full_is_feasible(const MultiplexGraph& graph,
                      const std::vector<std::uint8_t>& labels,
                      int k,
                      DominationMode mode) {
    if (static_cast<int>(labels.size()) != graph.n || k < 1 || k > graph.layers()) {
        return false;
    }
    const std::uint8_t max_label = mode == DominationMode::Roman ? 2 : 3;
    if (std::any_of(labels.begin(), labels.end(),
                    [max_label](auto label) { return label > max_label; })) {
        return false;
    }
    int valid = 0;
    for (int s = 0; s < graph.layers(); ++s) {
        bool layer_valid = true;
        for (int i = 0; i < graph.n && layer_valid; ++i) {
            int count2 = 0;
            int count3 = 0;
            for (int neighbor : graph.adjacency[s][i]) {
                count2 += labels[neighbor] == 2;
                count3 += labels[neighbor] == 3;
            }
            const bool violating = mode == DominationMode::Roman
                ? labels[i] == 0 && count2 == 0
                : (labels[i] == 0 && count3 == 0 && count2 < 2) ||
                      (labels[i] == 1 && count2 + count3 == 0);
            if (violating) {
                layer_valid = false;
            }
        }
        valid += layer_valid;
    }
    return valid >= k;
}

}  // namespace kmdrdp

