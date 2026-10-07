#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

struct Graph {
    int n = 0;
    long long m = 0;
    std::vector<std::vector<int>> adj;
};

Graph load_graph(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path);
    Graph graph;
    input >> graph.n >> graph.m;
    graph.adj.assign(graph.n, {});
    int u, v;
    long long edges = 0;
    while (input >> u >> v) {
        --u; --v;
        if (u < 0 || v < 0 || u >= graph.n || v >= graph.n || u == v)
            throw std::runtime_error("invalid edge in " + path);
        graph.adj[u].push_back(v);
        graph.adj[v].push_back(u);
        ++edges;
    }
    if (edges != graph.m) throw std::runtime_error("edge-count mismatch in " + path);
    for (auto& neighbors : graph.adj) {
        std::sort(neighbors.begin(), neighbors.end());
        if (std::adjacent_find(neighbors.begin(), neighbors.end()) != neighbors.end())
            throw std::runtime_error("duplicate edge in " + path);
    }
    return graph;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: ppi_exact_centrality LAYER_FILE OUTPUT_CSV\n";
        return 2;
    }
    try {
        const Graph graph = load_graph(argv[1]);
        const int n = graph.n;
        std::vector<double> betweenness(n, 0.0), harmonic(n, 0.0);

#pragma omp parallel
        {
            std::vector<double> local_bc(n, 0.0);
            std::vector<std::vector<int>> pred(n);
            std::vector<int> dist(n), stack;
            std::vector<double> sigma(n), delta(n);
            std::queue<int> queue;
            stack.reserve(n);

#pragma omp for schedule(dynamic, 1)
            for (int source = 0; source < n; ++source) {
                for (auto& values : pred) values.clear();
                std::fill(dist.begin(), dist.end(), -1);
                std::fill(sigma.begin(), sigma.end(), 0.0);
                std::fill(delta.begin(), delta.end(), 0.0);
                stack.clear();
                while (!queue.empty()) queue.pop();
                dist[source] = 0;
                sigma[source] = 1.0;
                queue.push(source);
                double harmonic_sum = 0.0;

                while (!queue.empty()) {
                    const int v = queue.front();
                    queue.pop();
                    stack.push_back(v);
                    for (const int w : graph.adj[v]) {
                        if (dist[w] < 0) {
                            dist[w] = dist[v] + 1;
                            harmonic_sum += 1.0 / static_cast<double>(dist[w]);
                            queue.push(w);
                        }
                        if (dist[w] == dist[v] + 1) {
                            sigma[w] += sigma[v];
                            pred[w].push_back(v);
                        }
                    }
                }
                harmonic[source] = harmonic_sum / static_cast<double>(n - 1);
                while (!stack.empty()) {
                    const int w = stack.back();
                    stack.pop_back();
                    for (const int v : pred[w])
                        delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                    if (w != source) local_bc[w] += delta[w];
                }
            }
#pragma omp critical
            for (int i = 0; i < n; ++i) betweenness[i] += local_bc[i];
        }

        const double bc_scale = n > 2 ? 1.0 / (static_cast<double>(n - 1) * (n - 2)) : 0.0;
        std::ofstream output(argv[2]);
        if (!output) throw std::runtime_error("cannot create output file");
        output << "internal_id,degree,harmonic_closeness,betweenness\n";
        output << std::setprecision(17);
        for (int vertex = 0; vertex < n; ++vertex) {
            output << vertex + 1 << ',' << graph.adj[vertex].size() << ','
                   << harmonic[vertex] << ',' << betweenness[vertex] * bc_scale << '\n';
        }
        std::cout << "OK n=" << n << " m=" << graph.m << " output=" << argv[2] << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
