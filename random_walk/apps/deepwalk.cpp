//
// Created by Shixuan Sun on 11/3/20.
//

#include "creeper_main.h"

HyperParameter g_para = {ExecutionMode::Static, SamplingMethod::AliasSampling, 0, {nullptr, 0}};

struct DeepWalk {
    int l_;
    Graph* g_;

    DeepWalk() {}

    DeepWalk(int l, Graph* g) : l_(l), g_(g) {}

    DeepWalk(const DeepWalk& dw) : DeepWalk(dw.l_, dw.g_) {}

    inline double weight(WalkerMeta& w, intT begin, intT end, int64_t eid) {
        return g_->edge_weight_[eid];
    }

    inline bool update(WalkerMeta& w, intT begin, intT end, int64_t eid) {
        return w.length_ == l_;
    }

    inline double max_weight(WalkerMeta& w) {
        return 5.0;
    }
};

void execute(Graph &graph, InputParser &cmd_parser) {
    /**
     * Extract the parameters.
     *  1. -l: the length.
     */
    int l = 80;
    if (!cmd_parser.get_cmd_option("-l").empty()) {
        l = std::stoi(cmd_parser.get_cmd_option("-l"));
    }

    g_para.length_ = l;

    int em = static_cast<int>(ExecutionMode::Static);
    if (!cmd_parser.get_cmd_option("-em").empty()) {
        em = std::stoi(cmd_parser.get_cmd_option("-em"));
    }

    g_para.execution_ = static_cast<ExecutionMode>(em);

    int sm = static_cast<int>(SamplingMethod::AliasSampling);
    if (!cmd_parser.get_cmd_option("-sm").empty()) {
        sm = std::stoi(cmd_parser.get_cmd_option("-sm"));
    }

    g_para.sample_ = static_cast<SamplingMethod>(sm);

    /**
     * Initialize the walkers.
     */
    std::vector<WalkerMeta> walkers;
    generate_all_source(graph, walkers);


    /**
     * Execute the walkers.
     */
    compute(graph, walkers, DeepWalk(l, &graph));

    /**
     * Finalize the walkers as necessary.
     */
}