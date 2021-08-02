//
// Created by Shixuan Sun on 2020/10/17.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <utility>
#include <vector>
#include <algorithm>
#include "config/type_config.h"
#include "util/utility.h"
#include <bits/stdc++.h>

enum PartitionType {
    Random
};

class Graph {
public:
    /**
     * graph type.
     */
    bool is_directed_;
    bool is_vertex_labeled_;
    bool is_vertex_weighted_;
    bool is_edge_labeled_;
    bool is_edge_weighted_;
    bool is_edge_weight_prefix_summed_;
    bool is_edge_weight_alias_generated_;
    bool is_edge_weight_rejection_generated_;
    bool is_offset_pair_;

    /**
     * The meta information.
     */
    intT num_vertices_;
    int64_t num_edges_;
    intT max_degree_;
    intT node_with_max_degree_;

    /**
     * Graph partition.
     */
    int num_partitions_;
    intT* partition_offset_;

    /**
     * CSR representation.
     */
    int64_t* offset_;
    intT* adj_;
    intT* degree_;

    /**
     * Vertex property: label and weight.
     */
    int* vertex_label_;
    double* vertex_weight_;

    /**
     * Edge property: label and weight.
     */
     int* edge_label_;
     double* edge_weight_;

     /**
      * Prefix sum & offset pair.
      */
     double* edge_weight_prefix_sum_;
     std::pair<int64_t, int64_t>* offset_pair_;

     /**
      * Alias Method Table
      * After sampling a number n, if n smaller than the value_, choose table.first; otherwise choose table.second
      */
      AliasSlot* edge_weight_alias_table_;

     /**
      * Rejection sampling method.
      */
     double* edge_weight_rejection_;
     double* edge_weight_rejection_max_;

public:
    explicit Graph(bool is_vertex_labeled = false, bool is_vertex_weighted = false,
            bool is_edge_labeled = false, bool is_edge_weighted = false, bool is_prefix_summed = false,
            bool is_alias_generated = false, bool is_rejection_generated = false, bool is_offset_pair = false) :
            is_directed_(false),
            is_vertex_labeled_(is_vertex_labeled),
            is_vertex_weighted_(is_vertex_weighted),
            is_edge_labeled_(is_edge_labeled),
            is_edge_weighted_(is_edge_weighted),
            is_edge_weight_prefix_summed_(is_prefix_summed),
            is_edge_weight_alias_generated_(is_alias_generated),
            is_edge_weight_rejection_generated_(is_rejection_generated),
            is_offset_pair_(is_offset_pair),
            num_vertices_(0),
            num_edges_(0),
            node_with_max_degree_(0),
            max_degree_(0),
            num_partitions_(0),
            partition_offset_(nullptr),
            offset_(nullptr),
            adj_(nullptr),
            degree_(nullptr),
            vertex_label_(nullptr),
            vertex_weight_(nullptr),
            edge_label_(nullptr),
            edge_weight_(nullptr),
            edge_weight_prefix_sum_(nullptr),
            offset_pair_(nullptr),
            edge_weight_alias_table_(nullptr),
            edge_weight_rejection_(nullptr),
            edge_weight_rejection_max_(nullptr){
    }

    ~Graph() {
        delete[] partition_offset_;
        partition_offset_ = nullptr;
        delete[] offset_;
        offset_ = nullptr;
        delete[] adj_;
        adj_ = nullptr;
        delete[] degree_;
        degree_ = nullptr;
        delete[] vertex_label_;
        vertex_label_ = nullptr;
        delete[] vertex_weight_;
        vertex_weight_ = nullptr;
        delete[] edge_label_;
        edge_label_ = nullptr;
        delete[] edge_weight_;
        edge_weight_ = nullptr;
        delete[] edge_weight_prefix_sum_;
        edge_weight_prefix_sum_ = nullptr;
        delete[] offset_pair_;
        offset_pair_ = nullptr;
        delete[] edge_weight_alias_table_;
        edge_weight_alias_table_ = nullptr;
        delete[] edge_weight_rejection_;
        edge_weight_rejection_ = nullptr;
        delete[] edge_weight_rejection_max_;
        edge_weight_rejection_max_ = nullptr;
    }

public:
    inline int num_partitions() const {
        return num_partitions_;
    }

    inline intT num_vertices() const {
        return num_vertices_;
    }

    inline int64_t num_edges() const {
        return num_edges_;
    }

    inline double avg_degree() const {
        return num_edges_ / (double)num_vertices_;
    }

    inline intT max_degree() const {
        return max_degree_;
    }

    inline intT node_with_max_degree() const {
        return node_with_max_degree_;
    }

    inline std::pair<intT*, intT> neighbors(intT u) const {
        return {adj_ + offset_[u], degree(u)};
    }

    inline std::pair<int64_t, int64_t> edges(intT u) const {
        return {offset_[u], offset_[u + 1]};
    }

    inline bool is_neighbor(intT u, intT v) const {
        if (degree(u) < degree(v)) {
            std::swap(u, v);
        }
        auto nbrs = neighbors(v);

        int begin = 0;
        int end = nbrs.second - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (nbrs.first[mid] == u)
                return true;
            else if (nbrs.first[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    inline intT degree(intT u) const {
        return degree_[u];
    }

    inline int vertex_label(intT u) const {
        return vertex_label_[u];
    }

    inline double vertex_weight(intT u) const {
        return vertex_weight_[u];
    }

    inline std::pair<int*, intT> edge_labels(intT u) const {
        return {edge_label_ + offset_[u], degree(u)};
    }

    inline std::pair<double*, intT> edge_weights(intT u) const {
        return {edge_weight_ + offset_[u], degree(u)};
    }

    std::pair<intT, intT> get_partition(int pid) const {
        assert(0 <= pid && pid < num_partitions_);
        return {partition_offset_[pid], partition_offset_[pid + 1]};
    }

    int get_pid(intT u) const {
        auto position = std::upper_bound(partition_offset_, partition_offset_ + num_partitions_, u) - partition_offset_;
        return position - 1;
    }

    intT get_partition_size(int pid) const {
        assert(0 <= pid && pid < num_partitions_);
        return partition_offset_[pid + 1] - partition_offset_[pid];
    }

    inline bool is_member(intT u, int pid) const {
        assert(0 <= pid && pid < num_partitions_);
        return partition_offset_[pid] <= u && u < partition_offset_[pid + 1];
    }

    inline int neighbor_position(intT u, intT v) const {
        auto nbrs = neighbors(u);
        return Utility::binary_search(nbrs.first, 0, nbrs.second, v);

    }

    inline double edge_weight(intT u, intT v) const {
        if (degree(v) < degree(u))
            std::swap(u, v);

        int pos = neighbor_position(u, v);
        return edge_weight_[offset_[u] + pos];
    }

    inline double edge_weight(int64_t eid) const {
        return edge_weight_[eid];
    }

    inline int edge_label(intT u, intT v) const {
        if (degree(v) < degree(u))
            std::swap(u, v);

        int pos = neighbor_position(u, v);
        return edge_label_[offset_[u] + pos];
    }

    inline int edge_label(int64_t eid) const {
        return edge_label_[eid];
    }
public:
    /**
     * Graph I/O operations.
     */
    void load_csr(const std::string& graph_dir);

    void store_csr(const std::string& graph_dir);

    void load_partition_csr(const std::string& graph_dir);

    void store_partition_csr(const std::string& graph_dir);

    void load_edge_list(const std::string& graph_dir, char skip='#');

    void store_edge_list(const std::string& graph_dir);

    void store_edge_list_directed(const std::string& graph_dir);

    void store_edge_list_with_weight(const std::string& graph_dir);

    void print_metadata();

    /**
     * Sequential scan the csr, and divide it into a number of partitions such that the degree
     * sum of vertices in a partition is no more than target_num_edges_per_partition.
     * @param target_num_edges_per_partition - the target number of edges per partition
     */
    void sequential_partition(intT target_num_edges_per_partition);

    /**
     * Relabel the vertex id, and reorder the csr/label/weight accordingly. This function
     * can break the partition.
     * @param new_to_old - the mapping relationship from the new vertex id to old vertex id,
     *                      e.g., new_to_old[1] = 5 means that the new id 1 maps to the old id 5.
     */
    void relabel(const intT* new_to_old);
private:
    /**
     * Graph I/O helper functions.
     */
    void load_csr_degree_file(const std::string &file_path, intT *&degree);
    void load_csr_partition_file(const std::string &file_path, intT *&partition);
    void load_csr_adj_file(const std::string &file_path, const intT *degree, int64_t *&offset, intT *&adj);
    void load_csr_vertex_label_file(const std::string &file_path, int *&label);
    void load_csr_vertex_weight_file(const std::string &file_path, double *&weight);
    void load_csr_edge_label_file(const std::string &file_path, int *&label);
    void load_csr_edge_weight_file(const std::string &file_path, double *&weight);
    void load_csr_edge_weight_prefix_sum_file(const std::string &file_path, double *&prefix_sum);
    void load_csr_edge_weight_alias_table_file(const std::string &alias_table_file_path, AliasSlot *&alias_table);
    void load_csr_edge_weight_max_file(const std::string &max_file_path, double *&max_value);

    void load_edge_list(const std::string &file_path, std::vector<std::pair<intT, intT>> &lines,
                        intT &max_vertex_id, char skip);
    void load_edge_list(std::vector<std::pair<intT, intT>> &edge_list, intT max_vertex_id,
                        int64_t *&offset, intT *&degree, intT *&adj);

    void is_sort();

    void collect_metadata();
};


#endif //GRAPH_H
