//
// Created by Shixuan Sun on 11/20/20.
//

#include <omp.h>
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

    // Compute prefix sum.
    log_info("Start compute prefix sum...");

    auto start = std::chrono::high_resolution_clock::now();

    graph.is_edge_weight_prefix_summed_ = true;
    graph.edge_weight_prefix_sum_ = new double[graph.num_edges()];

#pragma omp parallel for schedule(dynamic, 1000)
    for (int i = 0; i < graph.num_vertices_; ++i) {
        auto length = (uint32_t)(graph.offset_[i + 1] - graph.offset_[i]);
        auto src = graph.edge_weight_ + graph.offset_[i];
        auto dst = graph.edge_weight_prefix_sum_ + graph.offset_[i];
        Utility::sequential_prefix_sum(src, dst, length);
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Preprocess time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);

    graph.is_edge_weighted_ = false;
    graph.store_csr(cmd_parser.get_cmd_option("-f"));
    log_info("Done...");

    return 0;
}

