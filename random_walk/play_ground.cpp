//
// Created by Shixuan Sun on 2020/10/18.
//

#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <cassert>
// #include <omp.h>
#include "config/type_config.h"
#include "util/graph/graph.h"
#include "util/log/log.h"
#include "util/io/io.h"
#include "util/utility.h"
#include "SFMT.h"

sfmt_t sfmt;

bool ITS = false;
// #pragma omp threadprivate(sfmt)

#define COLLECT_METRICS

#if defined(COLLECT_METRICS)

/**
 * Global
 */
std::vector<int64_t> vertex_access_distribution;

/**
 * BSP
 */
std::vector<int64_t> active_query_per_round;

/**
 * ASP
 */
std::vector<int64_t> total_active_query_per_load;
std::vector<int64_t> active_query_per_load;
std::vector<int64_t> moves_per_load;
std::vector<int64_t> partition_access_distribution;

#else

#endif

uint64_t total_num_moves = 0;
uint64_t num_load = 0;
uint64_t num_round = 0;

char message[256];

struct Walk {
    intT current_;
    int length_;
};

struct PPR : Walk {

};

struct DeepWalk : Walk {

};

struct Node2Vec : Walk {

};

struct MetaPath : Walk {

};

struct PartitionManager {
    std::pair<intT*, int> get_partition(int pid) {
        return {nullptr, 0};
    }
    int num_partitions() {
        return 0;
    }
    bool is_member(intT u, int pid) {
        return true;
    }
    int get_pid(intT u) {
        return 0;
    }
};

inline intT move(Graph &graph, intT u, sfmt_t *lsf) {
    auto neighbors = graph.neighbors(u);
    int selected;
    if(ITS){
        auto lower_position = graph.edge_weight_prefix_sum_ + graph.offset_[u];
        auto upper_position = graph.edge_weight_prefix_sum_ + graph.offset_[u+1]-1;
        double prefix_upper_bound = graph.edge_weight_prefix_sum_[graph.offset_[u+1]-1];
        double number = prefix_upper_bound*sfmt_genrand_uint32(&sfmt)/UINT32_MAX;
        auto position = std::upper_bound(lower_position, upper_position, number) - graph.edge_weight_prefix_sum_;
        selected = graph.adj_[position];
    }
    else{
        selected = sfmt_genrand_uint32(lsf) % neighbors.second;
    }
    std::cout << std::endl;
    return neighbors.first[selected];
}

int select_partition(std::vector<std::vector<Walk>>& partition_states) {
    int selected_pid = 0;
    int max_value = partition_states[0].size();

    for (int i = 1; i < partition_states.size(); ++i) {
        int cur_value = partition_states[i].size();

        if (cur_value > max_value) {
            max_value = cur_value;
            selected_pid = i;
        }
    }

    return selected_pid;
}


void gb_parallel_bsp(Graph& graph, std::vector<intT>& sources, int length) {
    std::vector<Walk> walks;
    for (auto u : sources) {
        walks.push_back({u, 1});
    }

    const int step_length = 16;
    std::pair<intT*, intT> degree_state[step_length];
    
    auto start = std::chrono::high_resolution_clock::now();

    int num_walks = sources.size();
    num_walks = num_walks - num_walks % step_length;
    
    for (int i = 0; i < num_walks; i += step_length) {
        for (int cur_length = 1; cur_length < length; ++cur_length) {
            degree_state[0] = graph.neighbors(walks[i].current_);
            walks[i].current_ = degree_state[0].first[sfmt_genrand_uint32(&sfmt) % degree_state[0].second];
            walks[i].length_ += 1;

            degree_state[1] = graph.neighbors(walks[i + 1].current_);
            walks[i + 1].current_ = degree_state[1].first[sfmt_genrand_uint32(&sfmt) % degree_state[1].second];
            walks[i + 1].length_ += 1;

            degree_state[2] = graph.neighbors(walks[i + 2].current_);
            walks[i + 2].current_ = degree_state[2].first[sfmt_genrand_uint32(&sfmt) % degree_state[2].second];
            walks[i + 2].length_ += 1;

            degree_state[3] = graph.neighbors(walks[i + 3].current_);
            walks[i + 3].current_ = degree_state[3].first[sfmt_genrand_uint32(&sfmt) % degree_state[3].second];
            walks[i + 3].length_ += 1;

            degree_state[4] = graph.neighbors(walks[i + 4].current_);
            walks[i + 4].current_ = degree_state[4].first[sfmt_genrand_uint32(&sfmt) % degree_state[4].second];
            walks[i + 4].length_ += 1;

            degree_state[5] = graph.neighbors(walks[i + 5].current_);
            walks[i + 5].current_ = degree_state[5].first[sfmt_genrand_uint32(&sfmt) % degree_state[5].second];
            walks[i + 5].length_ += 1;

            degree_state[6] = graph.neighbors(walks[i + 6].current_);
            walks[i + 6].current_ = degree_state[6].first[sfmt_genrand_uint32(&sfmt) % degree_state[6].second];
            walks[i + 6].length_ += 1;

            degree_state[7] = graph.neighbors(walks[i + 7].current_);
            walks[i + 7].current_ = degree_state[7].first[sfmt_genrand_uint32(&sfmt) % degree_state[7].second];
            walks[i + 7].length_ += 1;

            degree_state[8] = graph.neighbors(walks[i + 8].current_);
            walks[i + 8].current_ = degree_state[8].first[sfmt_genrand_uint32(&sfmt) % degree_state[8].second];
            walks[i + 8].length_ += 1;

            degree_state[9] = graph.neighbors(walks[i + 9].current_);
            walks[i + 9].current_ = degree_state[9].first[sfmt_genrand_uint32(&sfmt) % degree_state[9].second];
            walks[i + 9].length_ += 1;

            degree_state[10] = graph.neighbors(walks[i + 10].current_);
            walks[i + 10].current_ = degree_state[10].first[sfmt_genrand_uint32(&sfmt) %
                                                                degree_state[10].second];
            walks[i + 10].length_ += 1;

            degree_state[11] = graph.neighbors(walks[i + 11].current_);
            walks[i + 11].current_ = degree_state[11].first[sfmt_genrand_uint32(&sfmt) %
                                                                degree_state[11].second];
            walks[i + 11].length_ += 1;

            degree_state[12] = graph.neighbors(walks[i + 12].current_);
            walks[i + 12].current_ = degree_state[12].first[sfmt_genrand_uint32(&sfmt) %
                                                                degree_state[12].second];
            walks[i + 12].length_ += 1;

            degree_state[13] = graph.neighbors(walks[i + 13].current_);
            walks[i + 13].current_ = degree_state[13].first[sfmt_genrand_uint32(&sfmt) %
                                                                degree_state[13].second];
            walks[i + 13].length_ += 1;

            degree_state[14] = graph.neighbors(walks[i + 14].current_);
            walks[i + 14].current_ = degree_state[14].first[sfmt_genrand_uint32(&sfmt) %
                                                                degree_state[14].second];
            walks[i + 14].length_ += 1;

            degree_state[15] = graph.neighbors(walks[i + 15].current_);
            walks[i + 15].current_ = degree_state[15].first[sfmt_genrand_uint32(&sfmt) %
                                                                degree_state[15].second];
            walks[i + 15].length_ += 1;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Execution time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}


void parallel_bsp(Graph& graph, std::vector<intT>& sources, int length) {
    std::vector<Walk> walks;
    for (auto u : sources) {
        walks.push_back({u, 1});
    }

    auto start = std::chrono::high_resolution_clock::now();

    int num_walks = sources.size();
    for (auto cur_length = 1; cur_length < length; ++cur_length) {
// #pragma omp parallel for
        for (int i = 0; i < num_walks; ++i) {
            // Move one step
            Walk &w = walks[i];
            w.current_ = move(graph, w.current_, &sfmt);
            w.length_ += 1;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Execution time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}

void parallel_asp(Graph& graph, std::vector<intT>& sources, int length) {
    std::vector<Walk> walks;
    for (auto u : sources) {
        walks.push_back({u, 1});
    }

    int num_walks = walks.size();
    auto start = std::chrono::high_resolution_clock::now();

// #pragma omp parallel for
    for (int i = 0; i < num_walks; ++i) {
        Walk& w = walks[i];
        while (w.length_ < length) {
            w.current_ = move(graph, w.current_, &sfmt);
            w.length_ += 1;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Execution time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}

void parallel_partition_asp(Graph& graph, std::vector<intT>& sources, int length) {
//    int num_partition = graph.num_partitions();
//    int num_threads = omp_get_max_threads();
//
//    // Maintain the global states.
//    std::vector<std::vector<Walk>> partition_states(num_partition);
//
//    // Maintain the local states for each thread
//    std::vector<std::vector<std::vector<Walk>>> local_partition_states(num_threads);
//    for (int i = 0; i < num_threads; ++i) {
//        local_partition_states[i].resize(num_partition);
//    }
//
//    // Initialize the states of each partition.
//    for (auto u : sources) {
//        int pid = graph.get_pid(u);
//        partition_states[pid].push_back({u, 1});
//    }
//
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // Execute each query in the selected partition.
//    int num_queries = sources.size();
//    int num_completed_queries = 0;
//
//    while (num_completed_queries < num_queries) {
//        // Select the next partition to be executed.
//        int pid = select_partition(partition_states);
//        std::vector<Walk>& states = partition_states[pid];
//        int num_walks = states.size();
//
//// #pragma omp parallel
//        {
//            int tid = omp_get_thread_num();
//
//// #pragma omp for reduction(+:num_completed_queries)
//            for (int i = 0; i < num_walks; ++i) {
//                Walk &w = states[i];
//
//                while (true) {
//                    w.current_ = move(graph, w.current_, &sfmt);
//                    w.degree_ += 1;
//
//                    if (w.degree_ == length) {
//                        // The query is completed.
//                        num_completed_queries += 1;
//                        break;
//                    } else {
//                        if (!graph.is_member(w.current_, pid)) {
//                            // The query moves out current partition.
//                            int target_pid = graph.get_pid(w.current_);
//                            local_partition_states[tid][target_pid].push_back(w);
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//
//        // Synchronize
//        for (int i = 0; i < num_threads; ++i) {
//            for (int j = 0; j < num_partition; ++j) {
//                if (!local_partition_states[i][j].empty()) {
//                    partition_states[j].insert(partition_states[j].end(), local_partition_states[i][j].begin(),
//                            local_partition_states[i][j].end());
//
//                    local_partition_states[i][j].clear();
//                }
//            }
//        }
//        states.clear();
//    }
//
//    auto end = std::chrono::high_resolution_clock::now();
//    log_info("Execution time: %.6lf seconds",
//             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}

void BSP_DeepWalk_execute(Graph& graph, std::vector<intT>& sources, int length) {
    std::vector<Walk> walks;
    for (auto u : sources) {
        walks.push_back({u, 1});

#if defined(COLLECT_METRICS)
        vertex_access_distribution[u] += 1;
#endif

    }

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Walk> next_round_walks;

    for (auto cur_length = 1; cur_length < length; ++cur_length) {
#if defined(COLLECT_METRICS)
        active_query_per_round.emplace_back(walks.size());
#endif

        for (auto w : walks) {
            // Move one step
            w.current_ = move(graph, w.current_, nullptr);
            w.length_ += 1;
            next_round_walks.emplace_back(w);

#if defined(COLLECT_METRICS)
            vertex_access_distribution[w.current_] += 1;
#endif
        }

        walks.swap(next_round_walks);
        next_round_walks.clear();

        num_round += 1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Execution time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}

void ASP_DeepWalk_execute(Graph& graph, std::vector<intT>& sources, int length) {
    std::vector<Walk> walks;
    for (auto u : sources) {
        walks.push_back({u, 1});

#if defined(COLLECT_METRICS)
        vertex_access_distribution[u] += 1;
#endif
    }

    auto start = std::chrono::high_resolution_clock::now();

    for (auto w : walks) {
        while (w.length_ < length) {
            w.current_ = move(graph, w.current_, nullptr);
            w.length_ += 1;

#if defined(COLLECT_METRICS)
            vertex_access_distribution[w.current_] += 1;
#endif
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Execution time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}

void Partition_ASP_DeepWalk_execute(Graph& graph, std::vector<intT>& sources, int length) {
    int num_partition = graph.num_partitions();
    std::vector<std::vector<Walk>> partition_states(num_partition);

    // Initialize the states of each partition.
    for (auto u : sources) {
        int pid = graph.get_pid(u);
        partition_states[pid].push_back({u, 1});
    }

    auto start = std::chrono::high_resolution_clock::now();

    // Execute each query in the selected partition.
    int num_queries = sources.size();
    int num_completed_queries = 0;

    while (num_completed_queries < num_queries) {
        // Select the next partition to be executed.
        int pid = select_partition(partition_states);
        std::vector<Walk>& states = partition_states[pid];

        num_load += 1;

#if defined(COLLECT_METRICS)
        uint64_t num_moves = 0;
#endif

        // Execute the queries in this parition.
        while (!states.empty()) {
            Walk w = states.back();
            states.pop_back();

            while (true) {
                w.current_ = move(graph, w.current_, nullptr);
                w.length_ += 1;

#if defined(COLLECT_METRICS)
                num_moves += 1;
#endif
                if (w.length_ == length) {
                    // The query is completed.
                    num_completed_queries += 1;
                    break;
                }
                else {
                    if (!graph.is_member(w.current_, pid)) {
                        // The query moves out current partition.
                        partition_states[graph.get_pid(w.current_)].emplace_back(w);
                        break;
                    }
                }
            }
        }

#if defined(COLLECT_METRICS)
        moves_per_load.push_back(num_moves);
#endif
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Execution time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}

void generate_all_source(Graph &graph, std::vector<intT> &source, intT num_vertex) {
    for (intT i = 0; i < num_vertex; ++i) {
        if (graph.degree(i) != 0) {
            source.push_back(i);
        }
    }

    Utility::shuffle(source.data(), source.size());
}

void generate_random_source(Graph &graph, std::vector<intT> &source, int target_num_vertex, intT num_vertex) {
    while (source.size() < target_num_vertex) {
        intT u = rand() % num_vertex;

        if (graph.degree(u) != 0) {
            source.push_back(u);
        }
    }
}

void generate_single_source(Graph &graph, std::vector<intT> &source, int target_num_vertex, intT num_vertex) {
    bool is_set = false;

    while (!is_set) {
        intT u = rand() % num_vertex;
        if (graph.degree(u) != 0) {
            source.resize(target_num_vertex, u);
            is_set = true;
        }
    }
}

int main(int argc, char *argv[]) {
    std::string input_graph_dir(argv[1]);
    std::string input_method_type(argv[2]);
    std::string input_source_num(argv[3]);
    std::string input_source_type(argv[4]);
    std::string input_length(argv[5]);
    std::string input_sampling_type(argv[6]);
    std::string output_file_dir(argv[7]);

    log_info("Input graph path: %s", input_graph_dir.c_str());
    log_info("Input method type: %s", input_method_type.c_str());
    log_info("Input source num: %s", input_source_num.c_str());
    log_info("Input source type: %s", input_source_type.c_str());
    log_info("Input length: %s", input_length.c_str());
    log_info("Input sampling type: %s", input_sampling_type.c_str());
    log_info("Input result path: %s", output_file_dir.c_str());

    bool is_vertex_labeled = false, is_vertex_weighted = false;
    bool is_edge_labeled = false, is_edge_weighted = true;

    uint64_t target_num_vertex = std::stoul(input_source_num);
    int length = std::stoi(input_length);

    Graph graph(is_vertex_labeled, is_vertex_weighted, is_edge_labeled, is_edge_weighted);
    graph.load_partition_csr(input_graph_dir);
    graph.print_metadata();

    intT num_vertices = graph.num_vertices();
    vertex_access_distribution.resize(num_vertices, 0);

    std::vector<intT> source;

// #pragma omp parallel
    {
        sfmt_init_gen_rand(&sfmt, 1234);
    }

    if(input_sampling_type == "Inverse_Transform_Sampling"){
        ITS = true;
        auto prefix_sum_start = std::chrono::high_resolution_clock::now();

        graph.edge_weight_prefix_sum_ = (double*)malloc(sizeof(double)*graph.num_edges());
        for(int i = 0; i < graph.num_vertices();i++){
            Utility::sequential_prefix_sum<double>(graph.edge_weight_+graph.offset_[i], graph.edge_weight_prefix_sum_+graph.offset_[i], graph.degree(i));
        }

        auto prefix_sum_end = std::chrono::high_resolution_clock::now();
        log_info("Generate prefix sum time: %.6lf seconds",
                 std::chrono::duration_cast<std::chrono::nanoseconds>(prefix_sum_end - prefix_sum_start).count() / 1000000000.0);
    }

    log_info("Generate source...");

    auto start = std::chrono::high_resolution_clock::now();

    if (input_source_type == "single") {
        generate_single_source(graph, source, target_num_vertex, num_vertices);
    }
    else if (input_source_type == "multiple") {
        generate_random_source(graph, source, target_num_vertex, num_vertices);
    }
    else if (input_source_type == "all") {
        generate_all_source(graph, source, num_vertices);
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Generate source time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);

    log_info("Start walk...");

    if (input_method_type == "BSP_DeepWalk") {
        BSP_DeepWalk_execute(graph, source, length);
    }
    else if (input_method_type == "ASP_DeepWalk") {
        ASP_DeepWalk_execute(graph, source, length);
    }
    else if (input_method_type == "Partition_ASP_DeepWalk") {
        partition_access_distribution.resize(graph.num_partitions(), 0);
        Partition_ASP_DeepWalk_execute(graph, source, length);
    }
    else if (input_method_type == "PBSP") {
//        parallel_bsp(graph, source, length);
        gb_parallel_bsp(graph, source, length);
    }
    else if (input_method_type == "PASP") {
        parallel_asp(graph, source, length);
    }
    else if (input_method_type == "PPASP") {
        parallel_partition_asp(graph, source, length);
    }
    else {
        log_error("The method type %s cannot be supported.", input_method_type.c_str());
        exit(-1);
    }

    log_info("Num of round: %lu", num_round);
    log_info("Num of load: %lu", num_load);
    log_info("Num of move: %lu", total_num_moves);
    // Output
    log_info("Dump results into file...");
    std::string file_path_prefix = output_file_dir + "/" + input_method_type + "-"
            + input_source_type + "-" + input_source_num + "-" + input_length;
    std::string file_path = file_path_prefix + std::string("-vertex_access_distribution.bin");
    
    IO::write(file_path, vertex_access_distribution);

    file_path = file_path_prefix + std::string("-active_query_per_round.bin");
    IO::write(file_path, active_query_per_round);

    file_path = file_path_prefix + std::string("-total_active_query_per_load.bin");
    IO::write(file_path, total_active_query_per_load);

    file_path = file_path_prefix + std::string("-active_query_per_load.bin");
    IO::write(file_path, active_query_per_load);

    file_path = file_path_prefix + std::string("-moves_per_load.bin");
    IO::write(file_path, moves_per_load);

    file_path = file_path_prefix + std::string("-partition_access_distribution.bin");
    IO::write(file_path, partition_access_distribution);
    log_info("Done.");

    return 0;
}
