//
// Created by Shixuan Sun on 12/8/20.
//

#include "creeper_main.h"

HyperParameter g_para = {ExecutionMode::Static, SamplingMethod::AliasSampling, 0, {nullptr, 0}};

struct GenericWalk {
    int l_;
    Graph* g_;

    GenericWalk() {}

    GenericWalk(int l, Graph* g) : l_(l), g_(g) {}

    GenericWalk(const GenericWalk& gw) : GenericWalk(gw.l_, gw.g_) {}

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
     *  2. -st: the source type.
     *  3. -wn: the number of walks
     *  4. -em: the execution mode.
     *  5. -sm: the sampling method.
     */
    int l = 80;
    if (!cmd_parser.get_cmd_option("-l").empty()) {
        l = std::stoi(cmd_parser.get_cmd_option("-l"));
    }

    g_para.length_ = l;

    SourceType st = SourceType::All;
    if (!cmd_parser.get_cmd_option("-st").empty()) {
        st = static_cast<SourceType>(std::stoi(cmd_parser.get_cmd_option("-st")));
    }

    intT wn = graph.num_vertices();
    if (!cmd_parser.get_cmd_option("-wn").empty()) {
        intT temp = std::stoi(cmd_parser.get_cmd_option("-wn"));
        wn = temp == 0 ? wn : temp;
    }

    int em = static_cast<int>(ExecutionMode::Static);
    if (!cmd_parser.get_cmd_option("-em").empty()) {
        em = std::stoi(cmd_parser.get_cmd_option("-em"));
    }

    g_para.execution_ = static_cast<ExecutionMode>(em);

    int sm = static_cast<int>(SamplingMethod::InverseTransformationSampling);
    if (!cmd_parser.get_cmd_option("-sm").empty()) {
        sm = std::stoi(cmd_parser.get_cmd_option("-sm"));
    }

    g_para.sample_ = static_cast<SamplingMethod>(sm);

    /**
     * Initialize the walkers.
     */
    std::vector<WalkerMeta> walkers;

    if (st == SourceType::All) {
        generate_all_source(graph, walkers);
    }
    else if (st == SourceType::Single) {
        generate_single_source(graph, walkers, wn, 0);
    }
    else {
        generate_random_source(graph, walkers, wn, 0);
    }

    /**
     * Execute the walkers.
     */
    compute(graph, walkers, GenericWalk(l, &graph));

    /**
     * Finalize the walkers as necessary.
     */
}
