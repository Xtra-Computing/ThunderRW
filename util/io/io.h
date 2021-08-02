//
// Created by Shixuan Sun on 2020/7/1.
//

#ifndef XTRAGRAPHCOMPUTING_IO_H
#define XTRAGRAPHCOMPUTING_IO_H

#include <vector>
#include <string>
#include <fstream>
#include "util/log/log.h"

class IO {
public:
    /**
     * Vector IO.
     */
    template<typename type> static void write(const std::string& file_path, std::vector<type>& vec);
    template<typename type> static void read(const std::string& file_path, std::vector<type>& vec);

    /**
     * Array IO.
     */
    template<typename type> static void write(const std::string& file_path, type* src, uint32_t length);
    template<typename type> static void read(const std::string& file_path, type* &src, uint32_t& length);
};

template<typename type>
void IO::write(const std::string &file_path, std::vector<type> &vec) {
    std::ofstream ofs(file_path, std::ios::binary);
    uint32_t element_size = sizeof(type);
    uint32_t num_element = vec.size();
    uint64_t size = ((uint64_t)num_element) * element_size;

    ofs.write(reinterpret_cast<const char *>(&element_size), 4);
    ofs.write(reinterpret_cast<const char *>(&num_element), 4);
    ofs.write(reinterpret_cast<const char *>(&vec.front()), size);
}

template<typename type>
void IO::read(const std::string &file_path, std::vector<type> &vec) {
    std::ifstream ifs(file_path, std::ios::binary);
    uint32_t element_size;
    ifs.read(reinterpret_cast<char *>(&element_size), 4);

    if (sizeof(type) != element_size) {
        log_error("The type cannot match.");
        exit(-1);
    }

    uint32_t num_element;
    ifs.read(reinterpret_cast<char *>(&num_element), 4);
    vec.resize(num_element);

    uint64_t size = ((uint64_t) num_element) * element_size;
    ifs.read(reinterpret_cast<char *>(&vec.front()), size);
}

template<typename type>
void IO::write(const std::string &file_path, type *src, uint32_t length) {
    std::ofstream ofs(file_path, std::ios::binary);
    uint32_t element_size = sizeof(type);
    uint64_t size = ((uint64_t)length) * element_size;

    ofs.write(reinterpret_cast<const char *>(&element_size), 4);
    ofs.write(reinterpret_cast<const char *>(&length), 4);
    ofs.write(reinterpret_cast<const char *>(src), size);
}

template<typename type>
void IO::read(const std::string &file_path, type *&src, uint32_t &length) {
    std::ifstream ifs(file_path, std::ios::binary);
    uint32_t element_size;
    ifs.read(reinterpret_cast<char *>(&element_size), 4);

    if (sizeof(type) != element_size) {
        log_error("The type cannot match.");
        exit(-1);
    }

    ifs.read(reinterpret_cast<char *>(&length), 4);
    uint64_t size = ((uint64_t) length) * element_size;
    src = (type*)malloc(size);
    ifs.read(reinterpret_cast<char *>(src), size);
}

#endif //XTRAGRAPHCOMPUTING_IO_H
