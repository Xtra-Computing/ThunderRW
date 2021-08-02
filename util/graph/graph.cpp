//
// Created by Shixuan Sun on 2020/10/17.
//
#include <fstream>
#include <chrono>
#include <sys/mman.h>
#include <fcntl.h>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <iostream>
#include "graph.h"
#include "util/log/log.h"

const char* degree_file_name = "b_degree.bin";
const char* partition_file_name = "b_partition.bin";
const char* adj_file_name = "b_adj.bin";
const char* vertex_label_file_name = "b_vertex_label.bin";
const char* vertex_weight_file_name = "b_vertex_weight.bin";
const char* edge_label_file_name = "b_edge_label.bin";
const char* edge_weight_file_name = "b_edge_weight.bin";
const char* edge_list_file_name = "b_edge_list.bin";
const char* edge_list_directed_file_name = "b_edge_list_directed.bin";
const char* edge_list_with_weight_file_name = "b_edge_list_with_weight.bin";
const char* edge_weight_prefix_sum_file_name = "b_edge_weight_prefix_sum.bin";
const char* edge_weight_alias_table_file_name = "b_edge_weight_alias_table.bin";
const char* edge_weight_max_file_name = "b_edge_weight_max.bin";

void Graph::load_partition_csr(const std::string &graph_dir){
    log_info("Load partition in CSR from %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();

    std::string partition_file_path = graph_dir + "/" + partition_file_name;
    load_csr_partition_file(partition_file_path, partition_offset_);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load partition time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
    load_csr(graph_dir);
}

void Graph::load_csr(const std::string &graph_dir) {
    log_info("Load graph in CSR from %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();

    std::string degree_file_path = graph_dir + "/" + degree_file_name;
    std::string adj_file_path = graph_dir + "/" + adj_file_name;
    std::string vertex_label_file_path = graph_dir + "/" + vertex_label_file_name;
    std::string vertex_weight_file_path = graph_dir + "/" + vertex_weight_file_name;
    std::string edge_label_file_path = graph_dir + "/" + edge_label_file_name;
    std::string edge_weight_file_path = graph_dir + "/" + edge_weight_file_name;
    std::string edge_weight_prefix_sum_file_path = graph_dir + "/" + edge_weight_prefix_sum_file_name;
    std::string edge_weight_alias_table_path = graph_dir + "/" + edge_weight_alias_table_file_name;
    std::string edge_weight_max_file_path = graph_dir + "/" + edge_weight_max_file_name;

    load_csr_degree_file(degree_file_path, degree_);
    load_csr_adj_file(adj_file_path, degree_, offset_, adj_);

    if (is_vertex_labeled_) {
        load_csr_vertex_label_file(vertex_label_file_path, vertex_label_);
    }

    if (is_vertex_weighted_) {
        load_csr_vertex_weight_file(vertex_weight_file_path, vertex_weight_);
    }

    if (is_edge_labeled_) {
        load_csr_edge_label_file(edge_label_file_path, edge_label_);
    }

    if (is_edge_weighted_) {
        load_csr_edge_weight_file(edge_weight_file_path, edge_weight_);
    }

    if (is_edge_weight_prefix_summed_) {
        load_csr_edge_weight_prefix_sum_file(edge_weight_prefix_sum_file_path, edge_weight_prefix_sum_);
    }

    if (is_edge_weight_alias_generated_) {
        load_csr_edge_weight_alias_table_file(edge_weight_alias_table_path, edge_weight_alias_table_);
    }

    if (is_offset_pair_) {
        offset_pair_ = new std::pair<int64_t, int64_t>[num_vertices_];
        for (auto i = 0; i < num_vertices_; ++i) {
            offset_pair_[i] = {offset_[i], offset_[i + 1]};
        }
    }

    if (is_edge_weight_rejection_generated_) {
        load_csr_edge_weight_max_file(edge_weight_max_file_path, edge_weight_rejection_max_);
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load graph time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
    collect_metadata();
}

void Graph::store_partition_csr(const std::string &graph_dir) {
    store_csr(graph_dir);

    log_info("Store graph partition as CSR into %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    std::string partition_file_path = graph_dir + "/" + partition_file_name;
    {
        log_info("Store partition...");
        std::ofstream partition_file(partition_file_path, std::ios::binary);
        uint32_t element_size = sizeof(intT);
        uint64_t size = (uint64_t)element_size * (num_partitions_ + 1);
        partition_file.write(reinterpret_cast<const char *>(&element_size), 4);
        partition_file.write(reinterpret_cast<const char *>(&num_partitions_), sizeof(int));
        partition_file.write(reinterpret_cast<const char *>(partition_offset_), size);
    }
    auto end = std::chrono::high_resolution_clock::now();
    log_info("Store graph partition time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::store_csr(const std::string &graph_dir) {
    log_info("Store graph as CSR into %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();

    std::string degree_file_path = graph_dir + "/" + degree_file_name;
    std::string adj_file_path = graph_dir + "/" + adj_file_name;
    std::string vertex_label_file_path = graph_dir + "/" + vertex_label_file_name;
    std::string vertex_weight_file_path = graph_dir + "/" + vertex_weight_file_name;
    std::string edge_label_file_path = graph_dir + "/" + edge_label_file_name;
    std::string edge_weight_file_path = graph_dir + "/" + edge_weight_file_name;
    std::string edge_weight_prefix_sum_file_path = graph_dir + "/" + edge_weight_prefix_sum_file_name;
    std::string edge_weight_alias_table_file_path = graph_dir + "/" + edge_weight_alias_table_file_name;
    std::string edge_weight_max_file_path = graph_dir + "/" + edge_weight_max_file_name;
    {
        log_info("Store degree...");
        std::ofstream degree_file(degree_file_path, std::ios::binary);
        uint32_t element_size = sizeof(intT);
        uint64_t size = (uint64_t)element_size * num_vertices_;
        degree_file.write(reinterpret_cast<const char *>(&element_size), 4);
        degree_file.write(reinterpret_cast<const char *>(&num_vertices_), element_size);
        degree_file.write(reinterpret_cast<const char *>(&num_edges_), sizeof(int64_t));
        degree_file.write(reinterpret_cast<const char *>(degree_), size);
    }
    {
        log_info("Store adj...");
        std::ofstream adj_file(adj_file_path, std::ios::binary);
        uint64_t size = sizeof(intT) * num_edges_;
        adj_file.write(reinterpret_cast<const char *>(adj_), size);
    }
    {
        if (is_vertex_labeled_) {
            log_info("Store vertex label...");
            std::ofstream vertex_label_file(vertex_label_file_path, std::ios::binary);
            uint64_t size = sizeof(int) * num_vertices_;
            vertex_label_file.write(reinterpret_cast<const char *>(vertex_label_), size);
        }
    }
    {
        if (is_vertex_weighted_) {
            log_info("Store vertex weight...");
            std::ofstream vertex_weight_file(vertex_weight_file_path, std::ios::binary);
            uint64_t size = sizeof(double) * num_vertices_;
            vertex_weight_file.write(reinterpret_cast<const char *>(vertex_weight_), size);
        }
    }
    {
        if (is_edge_labeled_) {
            log_info("Store edge label...");
            std::ofstream edge_label_file(edge_label_file_path, std::ios::binary);
            uint64_t size = sizeof(int) * num_edges_;
            edge_label_file.write(reinterpret_cast<const char *>(edge_label_), size);
        }
    }
    {
        if (is_edge_weighted_) {
            log_info("Store edge weight...");
            std::ofstream edge_weight_file(edge_weight_file_path, std::ios::binary);
            uint64_t size = sizeof(double) * num_edges_;
            edge_weight_file.write(reinterpret_cast<const char *>(edge_weight_), size);
        }
    }
    {
        if (is_edge_weight_prefix_summed_) {
            log_info("Store edge weight prefix sum...");
            std::ofstream edge_weight_prefix_sum_file(edge_weight_prefix_sum_file_path, std::ios::binary);
            uint64_t size = sizeof(double) * num_edges_;
            edge_weight_prefix_sum_file.write(reinterpret_cast<const char *>(edge_weight_prefix_sum_), size);
        }
    }
    {
        if (is_edge_weight_alias_generated_) {
            log_info("Store edge weight alias table...");
            std::ofstream edge_weight_alias_table_file(edge_weight_alias_table_file_path, std::ios::binary);
            uint64_t size = sizeof(AliasSlot) * num_edges_;
            edge_weight_alias_table_file.write(reinterpret_cast<const char *>(edge_weight_alias_table_), size);
        }
    }
    {
        if (is_edge_weight_rejection_generated_) {
            log_info("Store edge weight max...");
            std::ofstream edge_weight_max_file(edge_weight_max_file_path, std::ios::binary);
            uint64_t size = sizeof(double) * num_vertices_;
            edge_weight_max_file.write(reinterpret_cast<const char *>(edge_weight_rejection_max_), size);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Store graph time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_csr_partition_file(const std::string &file_path, intT *&partition) {
    log_info("Load partition file...");
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream file(file_path, std::ios::binary);
    uint32_t element_size;
    file.read(reinterpret_cast<char *>(&element_size), 4);
    assert(element_size == sizeof(intT) && "The size of the element is not equal to that of intT.");
    file.read(reinterpret_cast<char *>(&num_partitions_), sizeof(int));

    int64_t size = (num_partitions_ + 1) * element_size;
    partition = (intT*)malloc(size);
    file.read(reinterpret_cast<char *>(partition), size);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load partition file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_csr_degree_file(const std::string &file_path, intT *&degree) {
    log_info("Load degree file...");
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream file(file_path, std::ios::binary);
    uint32_t element_size;
    file.read(reinterpret_cast<char *>(&element_size), 4);

    assert( element_size == sizeof(intT) && "The size of the element is not equal to that of intT.");
    file.read(reinterpret_cast<char *>(&num_vertices_), element_size);
    file.read(reinterpret_cast<char *>(&num_edges_), sizeof(int64_t));

    uint64_t size = ((uint64_t) num_vertices_) * element_size;
    degree = (intT*)malloc(size);
    file.read(reinterpret_cast<char *>(degree), size);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load degree file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_csr_adj_file(const std::string &file_path, const intT *degree, int64_t *&offset, intT *&adj) {
    log_info("Load adj file...");
    auto start = std::chrono::high_resolution_clock::now();

    uint64_t size = sizeof(int64_t) * (uint64_t) (num_vertices_ + 1);
    offset = (int64_t*) malloc(size);
    offset[0] = 0;
    for (auto i = 0; i < num_vertices_; ++i) {
        offset[i + 1] = offset[i] + degree[i];
    }

    size = sizeof(intT) * (uint64_t) (num_edges_ + 16);
    adj = (intT*) malloc(size);

    auto dst_v_fd = open(file_path.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
    auto *buffer = (intT *)mmap(0, static_cast<uint64_t >(num_edges_) * sizeof(intT), PROT_READ, MAP_PRIVATE, dst_v_fd, 0);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("malloc and sequential-scan time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

    // load dst vertices into the array
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0; i < num_vertices_; i++) {
        for (auto j = offset[i]; j < offset[i + 1]; j++) {
            adj[j] = buffer[j];
        }
    }
    munmap(buffer, static_cast<uint64_t >(num_edges_) * sizeof(intT));

    auto end2 = std::chrono::high_resolution_clock::now();
    log_info("Load adj file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count() / 1000.0);
}

void Graph::load_csr_vertex_label_file(const std::string &file_path, int *&label) {
    log_info("Load vertex label file...");
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream file(file_path, std::ios::binary);
    uint64_t size = ((uint64_t) num_vertices_) * sizeof(int);
    vertex_label_ = (int*)malloc(size);
    file.read(reinterpret_cast<char *>(vertex_label_), size);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load vertex label file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_csr_vertex_weight_file(const std::string &file_path, double *&weight) {
    log_info("Load vertex weight file...");
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream file(file_path, std::ios::binary);
    uint64_t size = ((uint64_t) num_vertices_) * sizeof(double);
    vertex_weight_ = (double*)malloc(size);
    file.read(reinterpret_cast<char *>(vertex_weight_), size);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load vertex weight file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_csr_edge_label_file(const std::string &file_path, int *&label) {
    log_info("Load edge label file...");
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t size = sizeof(int) * (uint64_t) (num_edges_ + 16);
    label = (int*) malloc(size);

    auto dst_v_fd = open(file_path.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
    auto *buffer = (int *)mmap(0, static_cast<uint64_t >(num_edges_) * sizeof(int), PROT_READ, MAP_PRIVATE, dst_v_fd, 0);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("malloc and sequential-scan time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

    // load dst vertices into the array
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0; i < num_vertices_; i++) {
        for (auto j = offset_[i]; j < offset_[i + 1]; j++) {
            label[j] = buffer[j];
        }
    }
    munmap(buffer, static_cast<uint64_t >(num_edges_) * sizeof(int));

    auto end2 = std::chrono::high_resolution_clock::now();
    log_info("Load edge label file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count() / 1000.0);
}

void Graph::load_csr_edge_weight_file(const std::string &file_path, double *&weight) {
    log_info("Load edge weight file...");
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t size = sizeof(double) * (uint64_t) (num_edges_ + 16);
    weight = (double*) malloc(size);

    auto dst_v_fd = open(file_path.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
    auto *buffer = (double *)mmap(0, static_cast<uint64_t >(num_edges_) * sizeof(double), PROT_READ, MAP_PRIVATE, dst_v_fd, 0);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("malloc and sequential-scan time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);


    // load dst vertices into the array
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0; i < num_vertices_; i++) {
        for (auto j = offset_[i]; j < offset_[i + 1]; j++) {
            weight[j] = buffer[j];
        }
    }
    munmap(buffer, static_cast<uint64_t >(num_edges_) * sizeof(double));

    auto end2 = std::chrono::high_resolution_clock::now();
    log_info("Load edge weight file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count() / 1000.0);
}

void Graph::load_csr_edge_weight_prefix_sum_file(const std::string &file_path, double *&prefix_sum) {
    log_info("Load edge weight prefix sum file...");
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t size = sizeof(double) * (uint64_t) (num_edges_ + 16);
    prefix_sum = (double*)malloc(size);

    auto dst_v_fd = open(file_path.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
    auto *buffer = (double *)mmap(0, static_cast<uint64_t >(num_edges_) * sizeof(double), PROT_READ, MAP_PRIVATE, dst_v_fd, 0);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("malloc and sequential-scan time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

    // load dst vertices into the array
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0; i < num_vertices_; i++) {
        for (auto j = offset_[i]; j < offset_[i + 1]; j++) {
            prefix_sum[j] = buffer[j];
        }
    }
    munmap(buffer, static_cast<uint64_t >(num_edges_) * sizeof(double));

    auto end2 = std::chrono::high_resolution_clock::now();
    log_info("Load edge weight prefix sum file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count() / 1000.0);
}

void Graph::load_csr_edge_weight_alias_table_file(const std::string &alias_table_file_path, AliasSlot *&alias_table) {
    log_info("Load edge weight alias table file...");
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t size = sizeof(AliasSlot) * (uint64_t) (num_edges_ + 16);
    alias_table = (AliasSlot *)malloc(size);

    auto dst_v_fd = open(alias_table_file_path.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
    auto *buffer = (AliasSlot*) mmap(0, static_cast<uint64_t >(num_edges_) * sizeof(AliasSlot),
                                      PROT_READ, MAP_PRIVATE, dst_v_fd, 0);
    auto end = std::chrono::high_resolution_clock::now();
    log_info("malloc and sequential-scan time: %.3lf seconds",
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

        // load dst vertices into the array
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0; i < num_vertices_; i++) {
        for (auto j = offset_[i]; j < offset_[i + 1]; j++) {
            alias_table[j] = buffer[j];
        }
    }

    munmap(buffer, static_cast<uint64_t >(num_edges_) * sizeof(AliasSlot));

    auto end2 = std::chrono::high_resolution_clock::now();
    log_info("Load edge weight alias table time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start).count() / 1000.0);

}

void Graph::load_edge_list(const std::string &file_path, std::vector<std::pair<intT, intT>> &lines, intT &max_vertex_id, char skip) {
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream ifs(file_path);
    max_vertex_id = 0;
    std::string tmp_str;
    intT src, dst;
    uint64_t line_count = 0;
    while (std::getline(ifs, tmp_str)) {
        line_count += 1;

        if (tmp_str[0] != skip) {
            std::stringstream ss(tmp_str);
            if (!(ss >> src >> dst)) {
                log_error("Cannot convert line %lu to edge.", line_count);
                exit(-1);
            }
            if (src > dst)
                std::swap(src, dst);
            lines.emplace_back(src, dst);

            if (src > max_vertex_id) {
                max_vertex_id = src;
            }
            if (dst > max_vertex_id) {
                max_vertex_id = dst;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load edge list file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_edge_list(std::vector<std::pair<intT, intT>> &edge_list, intT max_vertex_id, int64_t *&offset,
                           intT *&degree, intT *&adj) {
    auto start = std::chrono::high_resolution_clock::now();
    std::sort(edge_list.begin(), edge_list.end(),
              [](const std::pair<intT, intT> &left, const std::pair<intT, intT> &right) {
                  if (left.first == right.first) {
                      return left.second < right.second;
                  }
                  return left.first < right.first;
              });

    num_vertices_ = max_vertex_id + 1;
    num_edges_ = 0;
    std::vector<intT> degree_arr(num_vertices_, 0);
    std::vector<std::vector<intT>> adj_arr(num_vertices_);

    std::pair<intT, intT> prev_edge = std::make_pair(num_vertices_, num_vertices_);
    for (const auto &edge : edge_list) {
        // Remove parallel edges.
        if (prev_edge != edge) {
            prev_edge = edge;
            intT src, dst;
            std::tie(src, dst) = edge;
            // Remove self loops.
            if (src != dst) {
                degree_arr[src] += 1;
                adj_arr[src].emplace_back(dst);
                degree_arr[dst] += 1;
                adj_arr[dst].emplace_back(src);
                num_edges_ += 2;
            }
        }
    }

    uint64_t size = sizeof(intT) * num_vertices_;
    degree = (intT*) malloc(size);
    std::copy(degree_arr.begin(), degree_arr.end(), degree);

    size = sizeof(int64_t) * (num_vertices_ + 1);
    offset = (int64_t*) malloc(size);
    offset[0] = 0;
    for (auto i = 0; i < num_vertices_; ++i) {
        offset[i + 1] = offset[i] + degree[i];
    }

    size = sizeof(intT) * (num_edges_ + 16);
    adj = (intT*) malloc(size);

#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0; i < num_vertices_; ++i) {
        std::copy(adj_arr[i].begin(), adj_arr[i].end(), adj + offset[i]);
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Convert edge list to CSR time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::load_edge_list(const std::string &graph_dir, char skip) {
    log_info("Load graph from edge list from %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    std::string file_path = graph_dir + "/" + edge_list_file_name;
    std::vector<std::pair<intT, intT>> edge_list;
    intT max_vertex_id;
    load_edge_list(file_path, edge_list, max_vertex_id, skip);
    load_edge_list(edge_list, max_vertex_id, offset_, degree_, adj_);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("load graph time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

    is_sort();
    collect_metadata();
}

void Graph::store_edge_list(const std::string &graph_dir) {
    log_info("Store graph as edge list into %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    std::string file_path = graph_dir + "/" + edge_list_file_name;
    std::ofstream edge_list_file(file_path);
    edge_list_file << "#t" << num_vertices_ << ' ' << num_edges_/2 << '\n';

    for (auto src = 0; src < num_vertices_; ++src) {
        for (auto i = offset_[src]; i < offset_[src + 1]; ++i) {
            intT dst = adj_[i];
            if (src < dst)
                edge_list_file << src << ' ' << dst << '\n';
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Store graph time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::store_edge_list_directed(const std::string &graph_dir) {
    log_info("Store graph as edge list into %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    std::string file_path = graph_dir + "/" + edge_list_directed_file_name;
    std::ofstream edge_list_file(file_path);
    edge_list_file << "#t" << num_vertices_ << ' ' << num_edges_ << '\n';

    for (auto src = 0; src < num_vertices_; ++src) {
        for (auto i = offset_[src]; i < offset_[src + 1]; ++i) {
            intT dst = adj_[i];
            edge_list_file << src << ' ' << dst << '\n';
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Store graph time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::store_edge_list_with_weight(const std::string &graph_dir) {
    log_info("Store graph as edge list into %s", graph_dir.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    std::string file_path = graph_dir + "/" + edge_list_with_weight_file_name;
    std::ofstream edge_list_file(file_path);

    for (auto src = 0; src < num_vertices_; ++src) {
        for (auto i = offset_[src]; i < offset_[src + 1]; ++i) {
            intT dst = adj_[i];
            if (src < dst)
                edge_list_file << src << ' ' << dst << ' ' << edge_weight_[i] << '\n';
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Store graph time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void Graph::print_metadata() {
    log_info("Graph Metadata:");
    log_info("%s, %s, %s, %s, %s", is_directed_ ? "Directed" : "Undirected",
             is_vertex_labeled_ ? "vertex labeled" : "vertex unlabeled",
             is_vertex_weighted_ ? "vertex weighted" : "vertex unweighted",
             is_edge_labeled_ ? "edge labeled" : "edge unlabeled",
             is_edge_weighted_ ? "edge weighted" : "edge unweighted");
    double csr_memory_cost = (sizeof(intT) * num_edges() + sizeof(int64_t) * (num_vertices() + 1)) / (1024.0 * 1024.0 * 1024.0);
    log_info("|V| = %d, |E| = %lu, d = %.2f, node with max degree = %d, max degree = %d, CSR Memory Cost (GB) = %.6lf",
            num_vertices(), num_edges(), avg_degree(), node_with_max_degree(), max_degree(), csr_memory_cost);
    if (partition_offset_ != nullptr) {
        log_info("Num of partitions = %d, Partition size = %d", num_partitions(), get_partition_size(0));
    }
}

void Graph::collect_metadata() {
    max_degree_ = *std::max_element(degree_, degree_ + num_vertices_);
    node_with_max_degree_ = std::max_element(degree_, degree_ + num_vertices_) - degree_;
}

void Graph::sequential_partition(intT target_num_edges_per_partition) {
    std::vector<intT> partition_offset;
    partition_offset.reserve(512);

    partition_offset.push_back(0);

    int64_t local_sum = 0;
    for (auto u = 0; u < num_vertices(); ++u) {
        local_sum += degree(u);

        // Generate a new partition.
        if (local_sum > target_num_edges_per_partition) {
            partition_offset.push_back(u);
            local_sum = degree(u);
        }
    }
    partition_offset.push_back(num_vertices());

    num_partitions_ = partition_offset.size() - 1;
    partition_offset_ = new intT[partition_offset.size()];
    std::copy(partition_offset.begin(), partition_offset.end(), partition_offset_);
}

void Graph::relabel(const intT *new_to_old) {
    log_info("start relabel...");
    std::vector<intT> old_to_new(num_vertices());

    // Initialize old to new.
    for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
        old_to_new[new_to_old[new_id]] = new_id;
    }

    // Convert csr.
    {
        // Convert vertex.
        {
            std::vector<intT> old_degree(degree_, degree_ + num_vertices());
            std::vector<int> old_vertex_label;
            if (is_vertex_labeled_) {
                old_vertex_label.insert(old_vertex_label.end(), vertex_label_, vertex_label_ + num_vertices());
            }
            std::vector<double> old_vertex_weight;
            if (is_vertex_weighted_) {
                old_vertex_weight.insert(old_vertex_weight.end(), vertex_weight_, vertex_weight_ + num_vertices());
            }

            for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
                auto old_id = new_to_old[new_id];
                degree_[new_id] = old_degree[old_id];
                if (is_vertex_labeled_)
                    vertex_label_[new_id] = old_vertex_label[old_id];
                if (is_vertex_weighted_)
                    vertex_weight_[new_id] = old_vertex_weight[old_id];
            }

            old_degree.clear();
            old_degree.shrink_to_fit();
            old_vertex_label.clear();
            old_vertex_label.shrink_to_fit();
            old_vertex_weight.clear();
            old_vertex_weight.shrink_to_fit();
        }
        {
            // Convert adj.
            std::vector<intT> old_adj(adj_, adj_ + num_edges());

            int64_t offset = 0;
            for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
                auto base = offset_[new_to_old[new_id]];
                for (auto i = base; i < base + degree_[new_id]; ++i) {
                    adj_[offset++] = old_to_new[old_adj[i]];
                }
            }

            old_adj.clear();
            old_adj.shrink_to_fit();
        }
        if (is_edge_labeled_) {
            std::vector<int> old_edge_label(edge_label_, edge_label_ + num_edges());
            int64_t offset = 0;
            for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
                auto base = offset_[new_to_old[new_id]];
                for (auto i = base; i < base + degree_[new_id]; ++i) {
                    edge_label_[offset++] = old_edge_label[i];
                }
            }
            old_edge_label.clear();
            old_edge_label.shrink_to_fit();
        }

        if (is_edge_weighted_) {
            std::vector<double> old_edge_weight(edge_weight_, edge_weight_ + num_edges());
            int64_t offset = 0;
            for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
                auto base = offset_[new_to_old[new_id]];
                for (auto i = base; i < base + degree_[new_id]; ++i) {
                    edge_weight_[offset++] = old_edge_weight[i];
                }
            }
            old_edge_weight.clear();
            old_edge_weight.shrink_to_fit();
        }

        if (is_edge_labeled_ || is_edge_weighted_) {
            std::vector<intT> index(max_degree());
            std::vector<int> old_edge_label(max_degree());
            std::vector<double> old_edge_weight(max_degree());

            int64_t offset = 0;
            for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
                auto new_edge_weight = edge_weight_ + offset;
                auto new_edge_label = edge_label_ + offset;

                for (auto i = 0; i < degree_[new_id]; ++i) {
                    index[i] = i;
                    if (is_edge_labeled_)
                        old_edge_label[i] = new_edge_label[i];

                    if (is_edge_weighted_)
                        old_edge_weight[i] = new_edge_weight[i];
                }

                auto new_neighbors = adj_ + offset;
                std::sort(index.begin(), index.begin() + degree_[new_id], [new_neighbors](intT l, intT r) {
                    return new_neighbors[l] < new_neighbors[r];
                });

                for (auto i = 0; i < degree_[new_id]; ++i) {
                    if (is_edge_labeled_)
                        new_edge_label[i] = old_edge_label[index[i]];

                    if (is_edge_weighted_)
                        new_edge_weight[i] = old_edge_weight[index[i]];
                }

                offset += degree_[new_id];
            }
        }

        // Convert offset and sort neighbors.
        offset_[0] = 0;
        for (auto new_id = 0; new_id < num_vertices(); ++new_id) {
            offset_[new_id + 1] = degree_[new_id] + offset_[new_id];
            std::sort(adj_ + offset_[new_id], adj_ + offset_[new_id + 1]);
        }
    }
}

void Graph::is_sort() {
    for (auto u = 0; u < num_vertices(); ++u) {
        auto nbr = neighbors(u);
        if (nbr.second > 1) {
            for (int i = 1; i < nbr.second; ++i) {
                if (nbr.first[i - 1] > nbr.first[i]) {
                    log_error("The (%d, %d), %d is not sorted.", nbr.first[i - 1], nbr.first[i], u);
                    exit(-1);
                }
            }
        }
    }
}

void Graph::load_csr_edge_weight_max_file(const std::string &max_file_path, double *&max_value) {
    log_info("Load max weight file...");
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream file(max_file_path, std::ios::binary);
    uint64_t size = ((uint64_t)num_vertices_) * sizeof(double);
    max_value = (double*)malloc(size);
    file.read(reinterpret_cast<char *>(max_value), size);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Load vertex label file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}
