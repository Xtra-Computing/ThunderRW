#include <omp.h>
#include <queue>
#include "util/graph/graph.h"
#include "util/command_parser.h"
#include "util/log/log.h"

int main(int argc, char *argv[]) {
    InputParser cmd_parser(argc, argv);
    // Load graphs.
    Graph graph;
    graph.is_edge_weighted_ = true;

    graph.load_csr(cmd_parser.get_cmd_option("-f"));

    graph.print_metadata();

    log_info("Start generate alias table...");

    auto start = std::chrono::high_resolution_clock::now();

    graph.is_edge_weight_alias_generated_ = true;
    graph.edge_weight_alias_table_ = new AliasSlot[graph.num_edges()];

#pragma omp parallel default(none) shared(graph)
    {
        std::queue<intT> small;
        std::queue<intT> large;

        #pragma omp for schedule(dynamic, 1000)
        for (int i = 0; i < graph.num_vertices_; ++i) {
            auto src = graph.edge_weight_ + graph.offset_[i];
            auto dst_table = graph.edge_weight_alias_table_ + graph.offset_[i];
            auto adj_table = graph.adj_ + graph.offset_[i];
            auto edge_num = graph.degree(i);
            double weight_sum = 0, weight_avg;
            int small_idx, large_idx, index = 0;

            for (int j = 0; j < edge_num; j++) {
                weight_sum += src[j];
            }
            weight_avg = weight_sum / edge_num;
            for (int j = 0; j < edge_num; j++) {
                src[j] /= weight_avg;
                if (src[j] < 1)
                    small.push(j);
                else
                    large.push(j);
            }
            while ((!small.empty()) && (!large.empty())) {
                small_idx = small.front();
                large_idx = large.front();
                dst_table[index] = {src[small_idx], adj_table[small_idx], adj_table[large_idx]};
                index++;
                small.pop();
                src[large_idx] = src[large_idx] + src[small_idx] - 1;
                if (src[large_idx] < 1) {
                    large.pop();
                    small.push(large_idx);
                }
            }
            while (!small.empty()) {
                small_idx = small.front();
                small.pop();
                dst_table[index] = {1.0, adj_table[small_idx], adj_table[small_idx]};
                index++;
            }
            while (!large.empty()) {
                large_idx = large.front();
                large.pop();
                dst_table[index] = {1.0, adj_table[large_idx], adj_table[large_idx]};
                index++;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    log_info("Preprocess time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);

    graph.is_edge_weighted_ = false;
    graph.store_csr(cmd_parser.get_cmd_option("-f"));
    log_info("Done...");

    return 0;
}

