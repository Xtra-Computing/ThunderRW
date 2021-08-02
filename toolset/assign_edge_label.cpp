//
// Created by sunsx on 12/11/20.
//

#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <tuple>
#include <random>
#include "util/graph/graph.h"
using namespace std;

const std::string edge_label_file_name = "b_edge_label.bin";

int main(int argc, char *argv[]) {
    string input_file_path(argv[1]);
    string input_label_num(argv[2]);

    string output_file_path_name = input_file_path + "/" + edge_label_file_name;
    int label_num = std::stoi(input_label_num);

    Graph graph;
    graph.load_csr(input_file_path);
    graph.print_metadata();

    auto edge_num = graph.num_edges();

    cout << "start assign..." << endl;
    std::vector<int> labels(edge_num);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, label_num - 1);

    for (int i = 0; i < edge_num; ++i) {
        labels[i] = dis(gen);
    }

    cout << "start write..." << endl;
    ofstream label_ofs(output_file_path_name, ios::binary);
    label_ofs.write(reinterpret_cast<const char *>(&labels.front()), labels.size() * sizeof(int));
    cout << "finish..." << endl;

    return 0;
}


