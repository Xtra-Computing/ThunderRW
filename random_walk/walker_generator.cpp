//
// Created by Shixuan Sun on 11/3/20.
//

#include "walker_generator.h"
#include "util/utility.h"
#include "SFMT.h"

void generate_all_source(Graph &graph, std::vector<WalkerMeta> &walkers) {
    auto num_vertex =graph.num_vertices();
    std::vector<intT> source;

    for (intT i = 0; i < num_vertex; ++i) {
        if (graph.degree(i) != 0) {
            source.push_back(i);
        }
    }

    Utility::shuffle(source.data(), source.size());

    walkers.reserve(2048);
    for (intT i = 0; i < source.size(); ++i) {
        walkers.push_back({i, source[i], source[i], 1, nullptr});
    }
}

void generate_random_source(Graph &graph, std::vector<WalkerMeta> &walkers, int target_num_vertex, int seed) {
    auto num_vertex =graph.num_vertices();
    sfmt_t sfmt;
    sfmt_init_gen_rand(&sfmt, seed);

    walkers.reserve(2048);
    while (walkers.size() < target_num_vertex) {
        intT u = sfmt_genrand_uint32(&sfmt) % num_vertex;

        if (graph.degree(u) != 0) {
            walkers.push_back({(int)walkers.size(), u, u, 1, nullptr});
        }
    }
}

void generate_single_source(Graph &graph, std::vector<WalkerMeta> &walkers, int target_num_vertex, int seed) {
    bool is_set = false;
    auto num_vertex =graph.num_vertices();
    sfmt_t sfmt;
    sfmt_init_gen_rand(&sfmt, seed);

    intT u = 0;
    intT d = graph.degree(u);
    for (intT v = 1; v < graph.num_vertices(); ++v) {
        intT cur_d = graph.degree(v);
        if (cur_d > d) {
            d = cur_d;
            u = v;
        }
    }

    std::cout << u << ' ' << graph.degree(u) << std::endl;
    walkers.reserve(2048);
    for (int i = 0; i < target_num_vertex; ++i) {
        walkers.push_back({i, u, u, 1, nullptr});
    }
}
