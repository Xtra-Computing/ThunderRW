//
// Created by Shixuan Sun on 11/3/20.
//

#ifndef XTRAGRAPHCOMPUTING_WALKER_GENERATOR_H
#define XTRAGRAPHCOMPUTING_WALKER_GENERATOR_H

#include "util/graph/graph.h"
#include "walker.h"

void generate_all_source(Graph &graph, std::vector<WalkerMeta> &walkers);
void generate_random_source(Graph &graph, std::vector<WalkerMeta> &walkers, int target_num_vertex, int seed);
void generate_single_source(Graph &graph, std::vector<WalkerMeta> &walkers, int target_num_vertex, int seed);

#endif //XTRAGRAPHCOMPUTING_WALKER_GENERATOR_H
