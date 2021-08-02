//
// Created by Shixuan Sun on 2020/6/30.
//

#include "util/log/log.h"
#include "util/graph/graph.h"
#include "util/utility.h"

int main(int argc, char *argv[]) {
    std::string input_graph_folder(argv[1]);
    std::string output_graph_folder(argv[2]);
    int target_number_edges_per_partition = std::stoi(argv[3]);
    // By default, do not shuffle
    bool is_shuffle = false;
    if(argc > 4){
        is_shuffle = std::stoi(argv[4]) == 1;
    }

    Graph graph;
    graph.load_csr(input_graph_folder);
    graph.print_metadata();

    if (is_shuffle) {
        std::vector<intT> new_to_old(graph.num_vertices());
        for (auto i = 0; i < graph.num_vertices(); ++i) {
            new_to_old[i] = i;
        }

        Utility::shuffle(new_to_old.data(), new_to_old.size());
        graph.relabel(new_to_old.data());
    }

    graph.sequential_partition(target_number_edges_per_partition);
    graph.print_metadata();
    graph.store_partition_csr(output_graph_folder);
    log_info("Done.");
    return 0;
}
