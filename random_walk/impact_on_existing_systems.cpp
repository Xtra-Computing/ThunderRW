#include <vector>
#include <cstdlib>
#include <chrono>
#include <cassert>
#include <omp.h>
#include "util/command_parser.h"
#include "config/type_config.h"
#include "util/graph/graph.h"
#include "util/log/log.h"
#include "util/utility.h"
#include "SFMT.h"
#include "creeper.h"

enum ExecutionSetting {
    BSP = 0,        // BSP execution
    IBSP = 1,       // BSP execution with interleaving
    ASP = 2,        // ASP execution
    IASP = 3        // ASP execution with interleaving
};

enum Application {
    DeepWalk = 0,
    PPR = 1,
    Node2Vec = 2,
    MetaPath = 3
};

const char* settings[] = {"BSP", "IBSP", "ASP", "IASP"};
const char* applications[] = {"DeepWalk", "PPR", "Node2Vec", "MetaPath"};

struct InputParameter {
    ExecutionSetting execute_setting;
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

#define LOCAL_RING_SIZE 64
#define LOCAL_PREFETCH_HINT PREFETCH_HINT
InputParameter g_input_parameter;

void interleaving_alias_move(Graph &graph, BufferSlot* ring, sfmt_t *sfmt) {
    for (int i = 0; i < LOCAL_RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            slot.r_ = sfmt_genrand_uint32(sfmt);
            _mm_prefetch((void*)(graph.offset_pair_ + slot.w_.current_), LOCAL_PREFETCH_HINT);
        }
    }

    // Stage 2: generate the position & prefetch the alias slot.
    for (int i = 0; i < LOCAL_RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& offset = graph.offset_pair_[slot.w_.current_];
            slot.r_ = offset.first + (slot.r_ % (offset.second - offset.first));
            slot.dr_ = sfmt_genrand_real2(sfmt);
            _mm_prefetch((void*)(graph.edge_weight_alias_table_ + slot.r_), LOCAL_PREFETCH_HINT);
        }
    }

    // Stage 3: update the walker.
    for (int i = 0; i < LOCAL_RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& alias_slot = graph.edge_weight_alias_table_[slot.r_];
            slot.w_.current_ = slot.dr_ <= alias_slot.alias_value_ ?
                               alias_slot.first_ : alias_slot.second_;
            slot.w_.length_ += 1;
        }
    }
}

// __attribute__ ((noinline))
void alias_move(Graph &graph, WalkerMeta &w, sfmt_t *sfmt, AliasSlot* alias_table) {
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

void interleaving_parallel_asp(Graph& graph, std::vector<intT>& sources) {
    std::vector<WalkerMeta> walks;

    for (auto u : sources) {
        walks.push_back({(int)walks.size(), u, u, 1, nullptr});
    }

    int num_walks = walks.size();
    int num_threads = 1;
    int partition_size = num_walks / num_threads;

    std::vector<std::pair<int, int>> task_range;
    int current_position = 0;
    for (int i = 0; i < num_threads; ++i) {
        int end_position = current_position + partition_size;
        task_range.emplace_back(current_position, end_position);
        current_position = end_position;
    }
    task_range.back().second = num_walks;

    uint64_t total_step = 0;

    auto start = std::chrono::high_resolution_clock::now();

// #pragma omp parallel
    {
        auto local_start = std::chrono::high_resolution_clock::now();
        sfmt_t sfmt;
        int tid = omp_get_thread_num();
        sfmt_init_gen_rand(&sfmt, tid);
        const InputParameter para = g_input_parameter;

        uint64_t local_step = 0;

        int target_length = para.length;
        int next = task_range[tid].first;
        int end_position = task_range[tid].second;
        int local_num_walks = task_range[tid].second - task_range[tid].first;
        int local_completed_num_walks = 0;

        // Initialize ring buffer.
        assert(LOCAL_RING_SIZE < local_num_walks);
        BufferSlot r[LOCAL_RING_SIZE];
        for (int i = 0; i < LOCAL_RING_SIZE; ++i) {
            r[i].empty_ = false;
            r[i].w_ = walks[next++];
        }

// #pragma omp for
        for (int i = 0; i < num_threads; ++i) {
            while (local_completed_num_walks < local_num_walks) {
                interleaving_alias_move(graph, r, &sfmt);


                for (int j = 0; j < LOCAL_RING_SIZE; ++j) {
                    BufferSlot& slot = r[j];
                    if (!slot.empty_) {
                        local_step += 1;
                        if (slot.w_.length_ == target_length) {
                            slot.empty_ = true;
                            local_completed_num_walks += 1;
                        }
                    }

                    if (slot.empty_ && next < end_position) {
                        slot.empty_ = false;
                        slot.w_ = walks[next++];
                    }
                }
            }
        }

        auto local_end = std::chrono::high_resolution_clock::now();
        double local_execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(local_end - local_start).count() / 1000000000.0;
        log_info("Execution time (seconds): %.6lf", local_execution_time);
        log_info("Num of steps: %zu", local_step);
        log_info("Num of walks: %zu", local_completed_num_walks);
// #pragma omp atomic
        total_step += local_step;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Execution time: %.6lf seconds", execution_time);
    log_info("Num of steps: %zu", total_step);
    log_info("Throughput (steps per second): %.2lf", total_step / execution_time);
}

void parallel_asp(Graph& graph, std::vector<intT>& sources) {
    std::vector<WalkerMeta> walks;
    std::vector<intT> last_node;
    for (auto u : sources) {
        walks.push_back({(int)walks.size(), u, u, 1, nullptr});
    }

    int num_walks = walks.size();
    uint64_t total_step = 0;

    auto start = std::chrono::high_resolution_clock::now();

// #pragma omp parallel
    {
        auto local_start = std::chrono::high_resolution_clock::now();
        sfmt_t sfmt;
        int tid = omp_get_thread_num();
        sfmt_init_gen_rand(&sfmt, tid);
        const InputParameter para = g_input_parameter;

        uint64_t local_step = 0;
        uint64_t local_walks = 0;
        uint64_t sum = 0;
        int target_length = para.length;
// #pragma omp for
        for (int i = 0; i < num_walks; ++i) {
            auto w = walks[i];
            local_walks += 1;

            while (w.length_ < target_length) {
                alias_move(graph, w, &sfmt, graph.edge_weight_alias_table_ + graph.offset_[w.current_]);
                w.length_ += 1;
                local_step += 1;
                sum += w.current_;
            }
        }

        auto local_end = std::chrono::high_resolution_clock::now();
        double local_execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(local_end - local_start).count() / 1000000000.0;
        log_info("Sum %zu", sum);
        log_info("Execution time (seconds): %.6lf", local_execution_time);
        log_info("Num of steps: %zu", local_step);
        log_info("Num of walks: %zu", local_walks);
// #pragma omp atomic
        total_step += local_step;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Execution time: %.6lf seconds", execution_time);
    log_info("Num of steps: %zu", total_step);
    log_info("Throughput (steps per second): %.2lf", total_step / execution_time);
}

void interleaving_parallel_bsp(Graph& graph, std::vector<intT>& sources) {
    std::vector<WalkerMeta> walks;
    std::vector<intT> last_node;
    for (auto u : sources) {
        walks.push_back({(int)walks.size(), u, u, 1, nullptr});
    }

    int num_walks = walks.size();
    int num_threads = omp_get_max_threads();
    int partition_size = num_walks / num_threads;

    std::vector<std::pair<int, int>> task_range;
    int current_position = 0;
    for (int i = 0; i < num_threads; ++i) {
        int end_position = current_position + partition_size;
        task_range.emplace_back(current_position, end_position);
        current_position = end_position;
    }
    task_range.back().second = num_walks;

    uint64_t total_step = 0;

    auto start = std::chrono::high_resolution_clock::now();

// #pragma omp parallel
    {
        auto local_start = std::chrono::high_resolution_clock::now();
        sfmt_t sfmt;
        int tid = omp_get_thread_num();
        sfmt_init_gen_rand(&sfmt, tid);
        const InputParameter para = g_input_parameter;

        uint64_t local_step = 0;
        int target_length = para.length;
        int begin = task_range[tid].first;
        int end = task_range[tid].second;
        int round = (end - begin) / LOCAL_RING_SIZE;
        int remain_begin = begin + round * LOCAL_RING_SIZE;
        int remain_size = end - remain_begin;
        assert(LOCAL_RING_SIZE < end - begin);
        BufferSlot r[LOCAL_RING_SIZE];

        for (int current_length = 1; current_length < target_length; ++current_length) {
// #pragma omp for
            for (int i = 0; i < num_threads; ++i) {
                for (int j = begin; j < end; j += LOCAL_RING_SIZE) {
                    for (int k = 0; k < LOCAL_RING_SIZE; ++k) {
                        r[k].empty_ = false;
                        r[k].w_ = walks[j + k];
                    }

                    interleaving_alias_move(graph, r, &sfmt);
                    local_step += LOCAL_RING_SIZE;

                    for (int k = 0; k < LOCAL_RING_SIZE; ++k) {
                        r[k].empty_ = true;
                        walks[j + k] = r[k].w_;
                    }
                }

                for (int j = 0; j < remain_size; ++j) {
                    r[j].empty_ = false;
                    r[j].w_ = walks[remain_begin + j];
                    local_step += 1;
                }

                interleaving_alias_move(graph, r, &sfmt);

                for (int j = 0; j < remain_size; ++j) {
                    r[j].empty_ = true;
                    walks[remain_begin + j] = r[j].w_;
                }
            }
        }

        auto local_end = std::chrono::high_resolution_clock::now();
        double local_execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(local_end - local_start).count() / 1000000000.0;
        log_info("Execution time (seconds): %.6lf", local_execution_time);
        log_info("Num of steps: %zu", local_step);

// #pragma omp atomic
        total_step += local_step;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Execution time: %.6lf seconds", execution_time);
    log_info("Num of steps: %zu", total_step);
    log_info("Throughput (steps per second): %.2lf", total_step / execution_time);
}

void parallel_bsp(Graph& graph, std::vector<intT>& sources) {
    std::vector<WalkerMeta> walks;
    for (auto u : sources) {
        walks.push_back({(int)walks.size(), u, u, 1, nullptr});
    }

    int num_walks = walks.size();
    uint64_t total_step = 0;

    auto start = std::chrono::high_resolution_clock::now();

// #pragma omp parallel
    {
        auto local_start = std::chrono::high_resolution_clock::now();
        sfmt_t sfmt;
        sfmt_init_gen_rand(&sfmt, 0);
        const InputParameter para = g_input_parameter;

        uint64_t local_step = 0;
        uint64_t sum = 0;
        int target_length = para.length;

        for (int current_length = 1; current_length < target_length; ++current_length) {
// #pragma omp for
            for (int i = 0; i < num_walks; ++i) {
                auto w = walks[i];
                alias_move(graph, w, &sfmt, graph.edge_weight_alias_table_ + graph.offset_[w.current_]);
                w.length_ += 1;
                local_step += 1;
                sum += w.current_;
            }
        }

        auto local_end = std::chrono::high_resolution_clock::now();
        double local_execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(local_end - local_start).count() / 1000000000.0;
        log_info("Sum %zu", sum);
        log_info("Execution time (seconds): %.6lf", local_execution_time);
        log_info("Num of steps: %zu", local_step);

// #pragma omp atomic
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
    ExecutionSetting execute_setting = ASP;
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
    if (!cmd_parser.get_cmd_option("-es").empty()) {
        execute_setting = static_cast<ExecutionSetting>(std::stoi(cmd_parser.get_cmd_option("-es")));
    }

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

    g_input_parameter = {execute_setting, execute_mode, sampling_method, application,
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

    log_info("Execution Setting %s, Execution Mode: %s, Sampling Method: %s, Application: %s, Stop Probability: %lf, Length: %d, p: %d, q: %d",
            settings[g_input_parameter.execute_setting],
            execution[g_input_parameter.execute_mode], sampling[g_input_parameter.sampling_method],
            applications[g_input_parameter.application], g_input_parameter.stop_probability,
            g_input_parameter.length, g_input_parameter.p, g_input_parameter.q);

    graph.load_csr(cmd_parser.get_cmd_option("-f"));

    graph.print_metadata();


    std::vector<intT> source;


    log_info("Generate source...");

    generate_all_source(graph, source, graph.num_vertices());

    log_info("Start walk...");

    switch (g_input_parameter.execute_setting) {
        case BSP:
            parallel_bsp(graph, source);
            break;
        case IBSP:
            interleaving_parallel_bsp(graph, source);
            break;
        case ASP:
            log_info("Parallel asp...");
            parallel_asp(graph, source);
            break;
        case IASP:
            log_info("Interleaving parallel asp...");
            interleaving_parallel_asp(graph, source);
            break;
        default:
            std::cerr << "The execution setting cannot be supported." << std::endl;
            break;
    }

    log_info("Done.");

    return 0;
}

