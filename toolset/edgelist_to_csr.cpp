//
// Created by Shixuan Sun on 2020/6/30.
//
#include <iostream>
#include "util/log/log.h"
#include "util/graph/graph.h"

int main(int argc, char *argv[]) {
    std::string input_graph_folder(argv[1]);
    std::string output_graph_folder(argv[2]);
    char skip_character = '#';
    bool is_vertex_labeled = false, is_vertex_weighted = false;
    bool is_edge_labeled = false, is_edge_weighted = false;
    if (argc > 3) {
        skip_character = argv[3][0];
    }
    if (argc > 4){
        is_vertex_labeled = std::stoi(argv[4]) == 1;
    }
    if (argc > 5){
        is_vertex_weighted = std::stoi(argv[5]) == 1;
    }
    if (argc > 6){
        is_edge_labeled = std::stoi(argv[6]) == 1;
    }
    if (argc > 7){
        is_edge_weighted = std::stoi(argv[7]) == 1;
    }

    log_info("Skip character is %c", skip_character);

    Graph graph(is_vertex_labeled, is_vertex_weighted, is_edge_labeled, is_edge_weighted);
    graph.load_edge_list(input_graph_folder, skip_character);
    graph.print_metadata();
    graph.store_csr(output_graph_folder);
    log_info("Done.");

    return 0;
}
