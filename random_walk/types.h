//
// Created by Shixuan Sun on 11/28/20.
//

#ifndef XTRAGRAPHCOMPUTING_TYPES_H
#define XTRAGRAPHCOMPUTING_TYPES_H

#include "config/type_config.h"

// #define LOG_SEQUENCE
// #define COLLECT_PERF_COUNTER

#define PREFETCH_HINT _MM_HINT_T0

// #define TEST_AMAC

#define ENABLE_INTERLEAVING

#ifdef ENABLE_INTERLEAVING

#define RING_SIZE 64
#define SMALL_RING_SIZE RING_SIZE
#define SEARCH_RING_SIZE 32

#else

#define RING_SIZE 1
#define SMALL_RING_SIZE 1
#define SEARCH_RING_SIZE 1

#endif

/**
 * Hyper parameter.
 */

/**
 * Command line -em: 0 Uniform, 1 Static, 2 Dynamic
 */

enum ExecutionMode {
    Uniform = 0,
    Static = 1,
    Dynamic = 2
};
const char* execution[] = {"Uniform", "Static", "Dynamic"};

/**
 * Command line -sm: 0 Uniform Sampling, 1 Inverse Transformation Sampling, 2 Alias Sampling, 3 Rejection Sampling,
 *                   4 Max Weight Rejection Sampling
 */

enum SamplingMethod {
    UniformSampling = 0,
    InverseTransformationSampling = 1,
    AliasSampling = 2,
    RejectionSampling = 3,
    MaxWeightRejectionSampling = 4
};
const char* sampling[] = {"Uniform", "ITS", "Alias", "RejectionSampling", "MaxRejectionSampling"};

/**
 * Command line -st: 0 All, 1 Single, 2 Multiple
 */

enum SourceType {
    All = 0,
    Single = 1,
    Multiple = 2
};

struct HyperParameter {
    ExecutionMode execution_;
    SamplingMethod sample_;
    int length_;
    std::pair<intT*, int> meta_path_;
};

/**
 * Ring buffer.
 */

struct BufferSlot {
    bool empty_;
    WalkerMeta w_;
    intT local_id_;
    intT prev_;
    int64_t r_;
    double dr_;
    std::pair<int64_t, int64_t> offset_;
    intT* seq_ = nullptr;
    double* weight_ = nullptr;
    AliasSlot* alias_ = nullptr;
};

template<class F>
struct ThreadParameter {
    int id_;
    Graph* graph_;
    double* weight_buffer_;
    intT** seq_buffer_;
    AliasSlot* alias_buffer_;
    WalkerMeta* queries_;
    int num_walkers_;
    uint64_t step_count_;
    F f_;
};

extern HyperParameter g_para;
extern int g_num_threads;

#endif //XTRAGRAPHCOMPUTING_TYPES_H
