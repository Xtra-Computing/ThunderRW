//
// Created by Shixuan Sun on 10/31/19.
//

#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <tuple>
#include <random>
#include "util/graph/graph.h"
using namespace std;

const std::string edge_weight_file_name = "b_edge_weight.bin";

int main(int argc, char *argv[]) {
    string input_file_path(argv[1]);
    string input_weight_lower_bound(argv[2]);
    string input_weight_upper_bound(argv[3]);

    string output_file_path_name = input_file_path + "/" + edge_weight_file_name;
    double weight_lower_bound = std::stod(input_weight_lower_bound);
    double weight_upper_bound = std::stod(input_weight_upper_bound);

    Graph graph;
    graph.load_csr(input_file_path);
    graph.print_metadata();

    auto edge_num = graph.num_edges();

    cout << "start assign..." << endl;
    std::vector<double> weights(edge_num);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(weight_lower_bound, weight_upper_bound);

    for (int i = 0; i < edge_num; ++i) {
       weights[i] = dis(gen);
    }

    cout << "start write..." << endl;
    ofstream weights_ofs(output_file_path_name, ios::binary);
    weights_ofs.write(reinterpret_cast<const char *>(&weights.front()), weights.size() * sizeof(double));
    cout << "finish..." << endl;

    return 0;
}
