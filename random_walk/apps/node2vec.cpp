//
// Created by Shixuan Sun on 11/5/20.
//

#include "creeper_main.h"

HyperParameter g_para = {ExecutionMode::Dynamic, SamplingMethod::MaxWeightRejectionSampling, 0, {nullptr, 0}};

struct Node2Vec {
    int l_;
    double p_;
    double q_;
    double max_;
    Graph* g_;

    Node2Vec() {}

    Node2Vec(int l, double p, double q, Graph* g) : l_(l), p_(p), q_(q), g_(g) {
        max_ = std::max(1.0, std::max(q_, p_));
    }

    Node2Vec(const Node2Vec& n2v) : Node2Vec(n2v.l_, n2v.p_, n2v.q_, n2v.g_) {}

    inline double weight(WalkerMeta& w, intT begin, intT end, int64_t eid) {
        if (w.length_ == 1) {
            return max_;
        }
        else {
            int prev = w.seq_[w.length_ - 2];
            if (end == prev) {
                return p_;
            } else if (g_->is_neighbor(prev, end)) {
                return 1.0;
            } else {
                return q_;
            }
        }
    }

    inline bool update(WalkerMeta& w, intT begin, intT end, int64_t eid) {
        return w.length_ == l_;
    }

    inline double max_weight(WalkerMeta& w) {
        return max_;
    }
};

void execute(Graph &graph, InputParser &cmd_parser) {
    /**
     * Extract the parameters.
     *  1. -l: the target length.
     *  2. -p: the return parameter.
     *  3. -q: the in-out parameter.
     */

    int l = 80;
    if (!cmd_parser.get_cmd_option("-l").empty()) {
        l = std::stoi(cmd_parser.get_cmd_option("-l"));
    }

    g_para.length_ = l;

    double p = 2;
    if (!cmd_parser.get_cmd_option("-p").empty()) {
        p = std::stof(cmd_parser.get_cmd_option("-p"));
    }

    p = 1.0 / p;

    double q = 0.5;
    if (!cmd_parser.get_cmd_option("-q").empty()) {
        q = std::stof(cmd_parser.get_cmd_option("-q"));
    }

    q = 1.0 / q;

    int em = static_cast<int>(ExecutionMode::Dynamic);
    if (!cmd_parser.get_cmd_option("-em").empty()) {
        em = std::stoi(cmd_parser.get_cmd_option("-em"));
    }

    g_para.execution_ = static_cast<ExecutionMode>(em);

    int sm = static_cast<int>(SamplingMethod::MaxWeightRejectionSampling);
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

    compute(graph, walkers, Node2Vec(l, p, q, &graph));

    /**
     * Finalize the walkers as necessary.
     */
}