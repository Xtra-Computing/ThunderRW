//
// Created by Shixuan Sun on 11/5/20.
//

#include "creeper_main.h"

HyperParameter g_para = {ExecutionMode::Dynamic, SamplingMethod::InverseTransformationSampling, 0, {nullptr, 0}};

struct MetaPath {
    int l_;
    int schema_length_;
    int* schema_;
    Graph* g_;

    MetaPath() {}

    MetaPath(int l, int schema_length, int* schema, Graph* g) : l_(l), schema_length_(schema_length),
            schema_(schema), g_(g) {}

    MetaPath(const MetaPath& mp) : MetaPath(mp.l_, mp.schema_length_, mp.schema_, mp.g_) {}

    inline double weight(WalkerMeta& w, intT begin, intT end, int64_t eid) {
        auto label = g_->edge_label(eid);
        auto pos = (w.length_ - 1) % schema_length_;

        if (label != schema_[pos])
            return 0.0;
        else
            return 1.0;
    }

    inline bool update(WalkerMeta& w, intT begin, intT end, int64_t eid) {
        return w.length_ == l_;
    }

    inline double max_weight(WalkerMeta& w) {
        return 1.0;
    }
};

void execute(Graph &graph, InputParser &cmd_parser) {
    /**
     * Extract the parameters.
     *  1. -l: the target length.
     *  2. -s: the schema.
     */

    int l = 80;
    if (!cmd_parser.get_cmd_option("-l").empty()) {
        l = std::stoi(cmd_parser.get_cmd_option("-l"));
    }

    g_para.length_ = l;

    int em = 2;
    if (!cmd_parser.get_cmd_option("-em").empty()) {
        em = std::stoi(cmd_parser.get_cmd_option("-em"));
    }

    g_para.execution_ = static_cast<ExecutionMode>(em);

    int sm = 1;
    if (!cmd_parser.get_cmd_option("-sm").empty()) {
        sm = std::stoi(cmd_parser.get_cmd_option("-sm"));
    }

    g_para.sample_ = static_cast<SamplingMethod>(sm);

    int* schema = nullptr;
    int schema_length = 0;
    if (!cmd_parser.get_cmd_option("-s").empty()) {

        std::vector<std::string> tokens;
        Utility::split(cmd_parser.get_cmd_option("-s"), ',', tokens);
        schema_length = tokens.size();
        schema = new int[schema_length];

        for (int i = 0; i < tokens.size(); ++i) {
            schema[i] = std::stoi(tokens[i]);
            std::cout << schema[i] << ' ';
        }
        std::cout << std::endl;
    }
    assert(schema != nullptr);

    /**
     * Initialize the walkers.
     */
    std::vector<WalkerMeta> walkers;
    generate_all_source(graph, walkers);


    /**
     * Execute the walkers.
     */

    compute(graph, walkers, MetaPath(l, schema_length, schema, &graph));

    delete[] schema;

    /**
     * Finalize the walkers as necessary.
     */
}