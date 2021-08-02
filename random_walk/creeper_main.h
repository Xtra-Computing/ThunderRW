//
// Created by Shixuan Sun on 12/3/20.
//

#ifndef XTRAGRAPHCOMPUTING_CREEPER_MAIN_H
#define XTRAGRAPHCOMPUTING_CREEPER_MAIN_H

#include "creeper.h"

void execute(Graph &graph, InputParser &cmd_parser);

int g_num_threads = 1;

int main(int argc, char *argv[]) {
    assert(RING_SIZE % SMALL_RING_SIZE == 0);
#ifdef ENABLE_INTERLEAVING
    log_info("Step interleaving enabled.");
#else
    log_info("Step interleaving disabled");
#endif
    log_info("Task ring size %d, Small ring size %d, Search ring size %d", RING_SIZE, SMALL_RING_SIZE, SEARCH_RING_SIZE);

    InputParser cmd_parser(argc, argv);

    // Load graphs.
    Graph graph;

    if (cmd_parser.check_cmd_option_exists("-vl")) {
        graph.is_vertex_labeled_ = true;
    }

    if (cmd_parser.check_cmd_option_exists("-el")) {
        graph.is_edge_labeled_ = true;
    }

    if (cmd_parser.check_cmd_option_exists("-vw")) {
        graph.is_vertex_weighted_ = true;
    }

    if (cmd_parser.check_cmd_option_exists("-ew")) {
        graph.is_edge_weighted_ = true;
    }

    if (cmd_parser.check_cmd_option_exists("-ep")) {
        graph.is_edge_weight_prefix_summed_ = true;
    }

    if (cmd_parser.check_cmd_option_exists("-al")) {
        graph.is_edge_weight_alias_generated_ = true;
    }

    if (cmd_parser.check_cmd_option_exists("-rj")) {
        graph.is_edge_weight_rejection_generated_ = true;
    }

    g_num_threads = 1;
    if (!cmd_parser.get_cmd_option("-n").empty()) {
        g_num_threads = std::stoi(cmd_parser.get_cmd_option("-n"));
    }

    graph.is_offset_pair_ = true;

    graph.load_csr(cmd_parser.get_cmd_option("-f"));

    graph.print_metadata();

    execute(graph, cmd_parser);

    log_info("Done.");

    return 0;
}

#endif //XTRAGRAPHCOMPUTING_CREEPER_MAIN_H
