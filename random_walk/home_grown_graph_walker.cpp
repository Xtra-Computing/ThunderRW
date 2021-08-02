#include <iostream>
#include <numeric>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <cassert>
#include <omp.h>
#include <queue>
#include "util/command_parser.h"
#include "config/type_config.h"
#include "util/graph/graph.h"
#include "util/log/log.h"
#include "util/io/io.h"
#include "util/utility.h"
#include "SFMT.h"
#include "creeper.h"
#include "linux-perf-events.h"

enum Application {
    DeepWalk = 0,
    PPR = 1,
    Node2Vec = 2,
    MetaPath = 3
};

const char* applications[] = {"DeepWalk", "PPR", "Node2Vec", "MetaPath"};

struct Walk {
    intT current_;
    int length_;
};

struct InputParameter {
    ExecutionMode execute_mode;
    SamplingMethod sampling_method;
    Application application;
    double stop_probability;
    int length;
    double p;
    double q;
    int* schema;
    int schema_length;
};

InputParameter g_input_parameter;

inline void deep_walk_update_weight(Graph& graph, Walk& w, double* weight) {
    auto weight_base = graph.edge_weight_ + graph.offset_[w.current_];
    auto degree = graph.offset_[w.current_ + 1] - graph.offset_[w.current_];
    for (auto i = 0; i < degree; ++i) {
        weight[i] = weight_base[i];
    }
}

inline void ppr_update_weight(Graph& graph, Walk& w, double* weight) {
    auto weight_base = graph.edge_weight_ + graph.offset_[w.current_];
    auto degree = graph.offset_[w.current_ + 1] - graph.offset_[w.current_];
    for (auto i = 0; i < degree; ++i) {
        weight[i] = weight_base[i];
    }
}

inline double node2vec_weight_function(Graph& graph, Walk& w, intT next, intT prev, double p, double q) {
    if (w.length_ == 1) {
        return std::max(1.0, std::max(p, q));
    }
    else if (next == prev) {
        return p;
    }
    else if (graph.is_neighbor(next, prev)) {
        return 1.0;
    }
    else {
        return q;
    }
}

inline void node2vec_update_weight(Graph& graph, Walk& w, double* weight, intT prev, const InputParameter& para) {
    auto neighbors = graph.neighbors(w.current_);
    for (auto i = 0; i < neighbors.second; ++i) {
        auto v = neighbors.first[i];
        weight[i] = node2vec_weight_function(graph, w, v, prev, para.p, para.q);
    }
}

inline void metapath_update_weight(Graph& graph, Walk& w, double* weight, const int* schema, const int schema_length) {
    int label = schema[(w.length_ - 1) % schema_length];
    auto offset = graph.offset_pair_[w.current_];

    int j = 0;
    for (auto i = offset.first; i < offset.second; ++i) {
        weight[j++] = graph.edge_label_[i] == label ? 1 : 0;
    }
}

void uniform_move(Graph &graph, Walk& w, sfmt_t *sfmt) {
    auto neighbors = graph.neighbors(w.current_);
    auto random_value = sfmt_genrand_uint32(sfmt);
    auto selected_position = random_value % neighbors.second;
    w.current_ = neighbors.first[selected_position];
}

void its_move(Graph &graph, Walk &w, sfmt_t *sfmt, double *weight) {
    auto neighbor_base = graph.adj_ + graph.offset_[w.current_];
    auto degree = graph.offset_[w.current_ + 1] - graph.offset_[w.current_];

    double max_weight = weight[degree - 1];
    auto target_weight = sfmt_genrand_real2(sfmt) * max_weight;

    int selected_position = Utility::upper_bound_branchless(weight, 0, degree, target_weight);

    // assert(selected_position < degree);

    w.current_ = neighbor_base[selected_position];
}

void alias_move(Graph &graph, Walk &w, sfmt_t *sfmt, AliasSlot* alias_table) {
    auto degree = graph.offset_[w.current_ + 1] - graph.offset_[w.current_];

    // Select a position.
    auto r0 = sfmt_genrand_uint32(sfmt);
    auto p0 = r0 % degree;

    // Select the vertex.
    auto r1 = sfmt_genrand_real2(sfmt);
    auto vertex = r1 <= alias_table[p0].alias_value_ ? alias_table[p0].first_ : alias_table[p0].second_;

    // Find the neighbor and update.
    w.current_ = vertex;
}

void rj_move(Graph &graph, Walk &w, sfmt_t *sfmt, double *edge_weight, double max_weight) {
    auto degree = graph.offset_[w.current_ + 1] - graph.offset_[w.current_];
    auto neighbors = graph.adj_ + graph.offset_[w.current_];

    // Select a neighbor.
    double r;
    uint32_t p0;
    do {
        // Generate x position.
        auto r0 = sfmt_genrand_uint32(sfmt);
        p0 = r0 % degree;

        // Generate y position.
        r = sfmt_genrand_real2(sfmt) * max_weight;
    } while (r > edge_weight[p0]);

    // Find the neighbor and update.
    w.current_ = neighbors[p0];
}

void max_rj_move_for_n2v(Graph &graph, Walk &w, sfmt_t *sfmt, intT prev, double max_weight, const InputParameter& para) {
    auto degree = graph.offset_[w.current_ + 1] - graph.offset_[w.current_];
    auto neighbors = graph.adj_ + graph.offset_[w.current_];

    // Select a neighbor.
    double r, selected_r;
    uint32_t p0;
    intT v;
    do {
        // Generate x position.
        auto r0 = sfmt_genrand_uint32(sfmt);
        p0 = r0 % degree;
        v = neighbors[p0];

        // Generate y position.
        r = sfmt_genrand_real2(sfmt) * max_weight;
        selected_r = node2vec_weight_function(graph, w, v, prev, para.p, para.q);
    } while (r > selected_r);

    // Find the neighbor and update.
    w.current_ = v;
}


void static_execute(Graph& graph, std::vector<intT>& sources) {
    std::vector<Walk> walks;
    std::vector<intT> last_node;
    for (auto u : sources) {
        walks.push_back({u, 1});
        last_node.push_back(u);
    }

    int num_walks = walks.size();
    uint64_t total_step = 0;

    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel
    {
        auto local_start = std::chrono::high_resolution_clock::now();
        intT* seq = nullptr;
        sfmt_t sfmt;
        sfmt_init_gen_rand(&sfmt, 0);
        const InputParameter para = g_input_parameter;

        if (para.length > 0) {
            seq = new intT[para.length];
        }

        uint64_t local_step = 0;
        uint64_t local_walks = 0;

#pragma omp for
        for (int i = 0; i < num_walks; ++i) {
            auto w = walks[i];
            local_walks += 1;

            bool flag;
            do {
                if (w.length_ < para.length) {
                    seq[w.length_ - 1] = w.current_;
                }

                switch (para.sampling_method) {
                    case UniformSampling:
                        uniform_move(graph, w, &sfmt);
                        break;
                    case InverseTransformationSampling:
                        its_move(graph, w, &sfmt, graph.edge_weight_prefix_sum_ + graph.offset_[w.current_]);
                        break;
                    case AliasSampling:
                        alias_move(graph, w, &sfmt, graph.edge_weight_alias_table_ + graph.offset_[w.current_]);
                        break;
                    case RejectionSampling:
                        rj_move(graph, w, &sfmt, graph.edge_weight_ + graph.offset_[w.current_], graph.edge_weight_rejection_max_[w.current_]);
                        break;
                }

                w.length_ += 1;

                local_step += 1;

                switch (para.application) {
                    case MetaPath:
                    case DeepWalk:
                    case Node2Vec:
                        flag = w.length_ < para.length;
                        break;
                    case PPR:
                        flag = sfmt_genrand_real2(&sfmt) >= para.stop_probability;
                        break;
                }
            } while (flag);

            if (para.length > 0) {
                seq[para.length - 1] = w.current_;
            }
        }

        delete[] seq;
        auto local_end = std::chrono::high_resolution_clock::now();
        double local_execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(local_end - local_start).count() / 1000000000.0;
        log_info("Execution time (seconds): %.6lf", local_execution_time);
        log_info("Num of steps: %zu", local_step);
        log_info("Num of walks: %zu", local_walks);
#pragma omp atomic
        total_step += local_step;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Execution time: %.6lf seconds", execution_time);
    log_info("Num of steps: %zu", total_step);
    log_info("Throughput (steps per second): %.2lf", total_step / execution_time);

}

void dynamic_execute(Graph& graph, std::vector<intT>& sources) {
    std::vector<Walk> walks;
    std::vector<intT> last_node;
    for (auto u : sources) {
        walks.push_back({u, 1});
        last_node.push_back(u);
    }

    int num_walks = walks.size();

    auto start = std::chrono::high_resolution_clock::now();

    uint64_t total_step = 0;

#pragma omp parallel
    {
        auto local_start = std::chrono::high_resolution_clock::now();
        std::queue<intT> small;
        std::queue<intT> large;

        intT *seq = nullptr;
        double *weight = new double[graph.max_degree()];
        AliasSlot *alias_table = new AliasSlot[graph.max_degree()];

        sfmt_t sfmt;
        sfmt_init_gen_rand(&sfmt, 0);
        InputParameter para = g_input_parameter;

        if (g_input_parameter.schema != nullptr) {
            para.schema_length = g_input_parameter.schema_length;
            para.schema = new int[para.schema_length];
            std::copy(g_input_parameter.schema, g_input_parameter.schema + g_input_parameter.schema_length,
                      para.schema);
        }

        if (para.length > 0) {
            seq = new intT[para.length];
        }
        int tid = omp_get_thread_num();

        uint64_t local_step = 0;
        uint64_t local_walks = 0;
        uint64_t local_terminated_walks = 0;
#pragma omp for
        for (int i = 0; i < num_walks; ++i) {
            auto w = walks[i];
            local_walks += 1;

            bool flag;
            do {
                flag = true;
                if (w.length_ < para.length) {
                    seq[w.length_ - 1] = w.current_;
                }

                if (para.sampling_method != MaxWeightRejectionSampling) {
                    switch (para.application) {
                        case DeepWalk: {
                            deep_walk_update_weight(graph, w, weight);
                            break;
                        }
                        case PPR: {
                            ppr_update_weight(graph, w, weight);
                            break;
                        }
                        case Node2Vec: {
                            node2vec_update_weight(graph, w, weight,
                                                   w.length_ == 1 ? 0 : seq[w.length_ - 2], para);
                            break;
                        }
                        case MetaPath: {
                            metapath_update_weight(graph, w, weight, para.schema, para.schema_length);
                            break;
                        }
                    }
                }

                switch (para.sampling_method) {
                    case InverseTransformationSampling: {
                        its_in_place_initialization(weight, graph.degree(w.current_));
                        if (weight[graph.degree(w.current_) - 1] != 0) {
                            its_move(graph, w, &sfmt, weight);
                        } else {
                            local_terminated_walks += 1;
                            flag = false;
                        }
                        break;
                    }
                    case AliasSampling: {
                        if (alias_initialization(weight, graph.adj_ + graph.offset_[w.current_],
                                                 alias_table, graph.degree(w.current_), large, small)) {
                            alias_move(graph, w, &sfmt, alias_table);
                        } else {
                            local_terminated_walks += 1;
                            flag = false;
                        }
                        break;
                    }
                    case RejectionSampling: {
                        double max_weight = *std::max_element(weight, weight + graph.offset_[w.current_ + 1] -
                                                                      graph.offset_[w.current_]);
                        rj_move(graph, w, &sfmt, weight, max_weight);
                        break;
                    }
                    case MaxWeightRejectionSampling: {
                        double max_weight = std::max(1.0, std::max(para.p, para.q));
                        max_rj_move_for_n2v(graph, w, &sfmt, w.length_ == 1 ? 0 : seq[w.length_ - 2], max_weight, para);
                        break;
                    }
                }

                // By default, reset the flag.
                if (flag) {
                    w.length_ += 1;
                    local_step += 1;

                    switch (para.application) {
                        case MetaPath:
                        case DeepWalk:
                        case Node2Vec:
                            flag = w.length_ < para.length;
                            break;
                        case PPR:
                            flag = sfmt_genrand_real2(&sfmt) >= para.stop_probability;
                            break;
                    }
                }

            } while (flag);

            if (para.length > 0) {
                seq[para.length - 1] = w.current_;
            }

            //if (tid == 0 && local_walks % 1000 == 0) {
            //    log_info("%d, %d, %d", local_walks, local_terminated_walks, i);
            //}
        }

        delete[] para.schema;
        delete[] seq;
        delete[] weight;
        delete[] alias_table;

        auto local_end = std::chrono::high_resolution_clock::now();
        double local_execution_time =
                std::chrono::duration_cast<std::chrono::nanoseconds>(local_end - local_start).count() / 1000000000.0;

        log_info("Thread %d Execution time (seconds): %.6lf", tid, local_execution_time);
        log_info("Thread %d Num of steps: %zu", tid, local_step);
        log_info("Thread %d Num of walks: %zu", tid, local_walks);
        log_info("Thread %d Num of terminated walks: %zu", tid, local_terminated_walks);
#pragma omp atomic
        total_step += local_step;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Execution time: %.6lf seconds", execution_time);
    log_info("Num of steps: %zu", total_step);
    log_info("Throughput (steps per second): %.2lf", total_step / execution_time);
}

void generate_single_source(Graph &graph, std::vector<intT> &source, int target_num_vertex) {
    bool is_set = false;

    while (!is_set) {
        intT u = rand() % graph.num_vertices();
        if (graph.degree(u) != 0) {
            source.resize(target_num_vertex, u);
            is_set = true;
        }
    }
}

void generate_all_source(Graph &graph, std::vector<intT> &source, intT num_vertex) {
    for (intT i = 0; i < num_vertex; ++i) {
        if (graph.degree(i) != 0) {
            source.push_back(i);
        }
    }
    Utility::shuffle(source.data(), source.size());
}


void parse_parameter(InputParser& cmd_parser, Graph& graph) {
    ExecutionMode execute_mode = Static;
    SamplingMethod sampling_method = AliasSampling;
    Application application = DeepWalk;

    // PPR parameter.
    double stop_probability = 0.2;

    // DeepWalk and Node2Vec parameters.
    int length = 80;
    double p = 2.0;
    double q = 0.5;

    // Set application and execution.
    if (!cmd_parser.get_cmd_option("-em").empty()) {
        execute_mode = static_cast<ExecutionMode>(std::stoi(cmd_parser.get_cmd_option("-em")));
    }

    if (!cmd_parser.get_cmd_option("-sm").empty()) {
        sampling_method = static_cast<SamplingMethod>(std::stoi(cmd_parser.get_cmd_option("-sm")));
    }

    if (!cmd_parser.get_cmd_option("-app").empty()) {
        application = static_cast<Application>(std::stoi(cmd_parser.get_cmd_option("-app")));
    }

    // Set parameter.
    if (!cmd_parser.get_cmd_option("-p").empty()) {
        p = std::stof(cmd_parser.get_cmd_option("-p"));
    }

    if (!cmd_parser.get_cmd_option("-q").empty()) {
        q = std::stof(cmd_parser.get_cmd_option("-q"));
    }

    if (!cmd_parser.get_cmd_option("-l").empty()) {
        length = std::stoi(cmd_parser.get_cmd_option("-l"));
    }

    if (!cmd_parser.get_cmd_option("-sp").empty()) {
        stop_probability = std::stof(cmd_parser.get_cmd_option("-sp"));
    }

    int* schema = nullptr;
    int schema_length = 0;

    if (!cmd_parser.get_cmd_option("-s").empty()) {

        std::vector<std::string> tokens;
        Utility::split(cmd_parser.get_cmd_option("-s"), ',', tokens);
        schema_length = tokens.size();
        schema = new int[schema_length];

        for (int i = 0; i < tokens.size(); ++i) {
            schema[i] = std::stoi(tokens[i]);
        }
    }

    g_input_parameter = {execute_mode, sampling_method, application,
                         stop_probability, length, 1.0 / p, 1.0 / q, schema, schema_length};

    // Set the graph property.
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



    graph.is_offset_pair_ = true;
}


int main(int argc, char *argv[]) {
    InputParser cmd_parser(argc, argv);
    Graph graph;

    parse_parameter(cmd_parser, graph);

    log_info("Execution Mode: %s, Sampling Method: %s, Application: %s, Stop Probability: %lf, Length: %d, p: %d, q: %d",
            execution[g_input_parameter.execute_mode], sampling[g_input_parameter.sampling_method],
            applications[g_input_parameter.application], g_input_parameter.stop_probability,
            g_input_parameter.length, g_input_parameter.p, g_input_parameter.q);

    graph.load_csr(cmd_parser.get_cmd_option("-f"));

    graph.print_metadata();

    // compute memory cost
//    size_t memory_cost = 0;
//    for (intT u = 0; u < graph.num_vertices(); ++u) {
//        size_t d = graph.degree(u);
//        memory_cost += d*d;
//    }
//
//    std::cout << 4 * memory_cost / 1024.0 / 1024.0 / 1024.0 / 1024.0 << "TB" << std::endl;
//    exit(0);

    std::vector<intT> source;


    log_info("Generate source...");


    if (g_input_parameter.application != PPR) {
        generate_all_source(graph, source, graph.num_vertices());
    }
    else {
        generate_single_source(graph, source, graph.num_vertices());
        std::string my_source = cmd_parser.get_cmd_option("-source");
        int m_source = std::stoi(my_source);
        for(int i = 0; i < source.size(); i++){
            source[i] = m_source;
        }
    }

    log_info("Start walk...");

    if (g_input_parameter.execute_mode != Dynamic) {
        std::vector<int> evts;
        evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
        evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
        evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
        LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);

        std::vector<unsigned long long> creeper_results;
        creeper_results.resize(evts.size());

        unified.start();
        static_execute(graph, source);
        unified.end(creeper_results);
        log_info("%d cycles, %d instructions, %d LLC misses\n", creeper_results[0], creeper_results[1], creeper_results[2]);
    }
    else {
        dynamic_execute(graph, source);
    }

    log_info("Done.");

    return 0;
}

