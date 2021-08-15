//
// Created by Shixuan Sun on 11/3/20.
//

#ifndef XTRAGRAPHCOMPUTING_CREEPER_H
#define XTRAGRAPHCOMPUTING_CREEPER_H

#include <pthread.h>
#include <emmintrin.h>
#include "config/type_config.h"
#include "util/utility.h"
#include "walker_generator.h"
#include "its_frame.h"
#include "util/log/log.h"
#include "SFMT.h"
#include "util/command_parser.h"
#include "types.h"
#include "uniform_sampling.h"
#include "inverse_transformation_sampling.h"
#include "alias_sampling.h"
#include "rejection_sampling.h"
#include "linux-perf-events.h"

template<class F> void static_gather(Graph* graph, intT v, F &f, double* w) {
    // Loop over the neighbors of current vertex and apply the weight function to each of them.
    auto nbrs = graph->edges(v);
    auto base = graph->adj_ + nbrs.first;
    int64_t loop_count = nbrs.second - nbrs.first;
    WalkerMeta temp_w;

    switch (loop_count & 15) {
        case 15:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 14:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 13:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 12:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 11:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 10:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 9:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 8:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 7:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 6:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 5:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 4:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 3:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 2:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 1:
            w[loop_count - 1] = f.weight(temp_w, v, base[loop_count - 1], nbrs.first + loop_count - 1);
            loop_count -= 1;
        case 0:
            break;
    }

    for (auto x = 0; x < loop_count; x += 16) {
        w[x] = f.weight(temp_w, v, base[x], nbrs.first + x);
        w[x + 1] = f.weight(temp_w, v, base[x + 1], nbrs.first + x + 1);
        w[x + 2] = f.weight(temp_w, v, base[x + 2], nbrs.first + x + 2);
        w[x + 3] = f.weight(temp_w, v, base[x + 3], nbrs.first + x + 3);
        w[x + 4] = f.weight(temp_w, v, base[x + 4], nbrs.first + x + 4);
        w[x + 5] = f.weight(temp_w, v, base[x + 5], nbrs.first + x + 5);
        w[x + 6] = f.weight(temp_w, v, base[x + 6], nbrs.first + x + 6);
        w[x + 7] = f.weight(temp_w, v, base[x + 7], nbrs.first + x + 7);
        w[x + 8] = f.weight(temp_w, v, base[x + 8], nbrs.first + x + 8);
        w[x + 9] = f.weight(temp_w, v, base[x + 9], nbrs.first + x + 9);
        w[x + 10] = f.weight(temp_w, v, base[x + 10], nbrs.first + x + 10);
        w[x + 11] = f.weight(temp_w, v, base[x + 11], nbrs.first + x + 11);
        w[x + 12] = f.weight(temp_w, v, base[x + 12], nbrs.first + x + 12);
        w[x + 13] = f.weight(temp_w, v, base[x + 13], nbrs.first + x + 13);
        w[x + 14] = f.weight(temp_w, v, base[x + 14], nbrs.first + x + 14);
        w[x + 15] = f.weight(temp_w, v, base[x + 15], nbrs.first + x + 15);
    }
}


template<class F> void static_initialization(Graph* graph, F &f) {
    log_info("Initialize...");

    auto start = std::chrono::high_resolution_clock::now();

    if (g_para.sample_ == AliasSampling) {
        graph->edge_weight_alias_table_ = new AliasSlot[graph->num_edges()];
        graph->is_edge_weight_alias_generated_ = true;

#pragma omp parallel
        {
            F lf = f;
            std::queue<intT> small;
            std::queue<intT> large;

            double* local_buffer = new double[graph->max_degree()];

#pragma omp for schedule(dynamic, 1000)
            for (int i = 0; i < graph->num_vertices_; ++i) {
                auto offset = graph->offset_pair_[i];
                static_gather(graph, i, lf, local_buffer);
                alias_initialization(local_buffer, graph->adj_ + offset.first,
                                     graph->edge_weight_alias_table_ + offset.first, offset.second - offset.first,
                                     large, small);
            }

            delete[] local_buffer;
        }
    }
    else if (g_para.sample_ == InverseTransformationSampling) {
        graph->edge_weight_prefix_sum_ = new double[graph->num_edges()];
        graph->is_edge_weight_prefix_summed_ = true;

#pragma omp parallel
        {
            F lf = f;
            double* local_buffer = new double[graph->max_degree()];

#pragma omp for schedule(dynamic, 1000)
            for (int i = 0; i < graph->num_vertices_; ++i) {
                auto offset = graph->offset_pair_[i];
                static_gather(graph, i, lf, local_buffer);
                its_initialization(local_buffer, graph->edge_weight_prefix_sum_ + offset.first,
                        offset.second - offset.first);
            }

            delete[] local_buffer;

        }
    }
    else if (g_para.sample_ == RejectionSampling) {
        graph->edge_weight_rejection_ = new double[graph->num_edges()];
        graph->edge_weight_rejection_max_ = new double[graph->num_vertices()];
        graph->is_edge_weight_rejection_generated_ = true;

#pragma omp parallel
        {
            F lf = f;

#pragma omp for schedule(dynamic, 1000)
            for (int i = 0; i < graph->num_vertices_; ++i) {
                auto offset = graph->offset_pair_[i];
                static_gather(graph, i, lf, graph->edge_weight_rejection_ + offset.first);
                rs_initialization(graph->edge_weight_rejection_ + offset.first, offset.second - offset.first,
                                  graph->edge_weight_rejection_max_[i]);
            }

        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("Preprocessing time: %.6lf seconds",
             std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0);
}


template<class F>
void dynamic_gather(Graph *graph, BufferSlot *ring, F &f, HyperParameter &para, int &num_completed_walkers,
                    int &num_terminated_walkers, std::queue<intT> &large, std::queue<intT> &small) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto w = slot.weight_;
            auto current = slot.w_.current_;

            // Loop over the neighbors of current vertex and apply the weight function to each of them.
            auto nbrs = graph->edges(current);
            auto base = graph->adj_ + nbrs.first;
            int64_t loop_count = nbrs.second - nbrs.first;

            switch (loop_count & 15) {
                case 15:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 14:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 13:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 12:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 11:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 10:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 9:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 8:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 7:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 6:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 5:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 4:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 3:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 2:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 1:
                    w[loop_count - 1] = f.weight(slot.w_, current, base[loop_count - 1], nbrs.first + loop_count - 1);
                    loop_count -= 1;
                case 0:
                    break;
            }

            for (auto x = 0; x < loop_count; x += 16) {
                w[x] = f.weight(slot.w_, current, base[x], nbrs.first + x);
                w[x + 1] = f.weight(slot.w_, current, base[x + 1], nbrs.first + x + 1);
                w[x + 2] = f.weight(slot.w_, current, base[x + 2], nbrs.first + x + 2);
                w[x + 3] = f.weight(slot.w_, current, base[x + 3], nbrs.first + x + 3);
                w[x + 4] = f.weight(slot.w_, current, base[x + 4], nbrs.first + x + 4);
                w[x + 5] = f.weight(slot.w_, current, base[x + 5], nbrs.first + x + 5);
                w[x + 6] = f.weight(slot.w_, current, base[x + 6], nbrs.first + x + 6);
                w[x + 7] = f.weight(slot.w_, current, base[x + 7], nbrs.first + x + 7);
                w[x + 8] = f.weight(slot.w_, current, base[x + 8], nbrs.first + x + 8);
                w[x + 9] = f.weight(slot.w_, current, base[x + 9], nbrs.first + x + 9);
                w[x + 10] = f.weight(slot.w_, current, base[x + 10], nbrs.first + x + 10);
                w[x + 11] = f.weight(slot.w_, current, base[x + 11], nbrs.first + x + 11);
                w[x + 12] = f.weight(slot.w_, current, base[x + 12], nbrs.first + x + 12);
                w[x + 13] = f.weight(slot.w_, current, base[x + 13], nbrs.first + x + 13);
                w[x + 14] = f.weight(slot.w_, current, base[x + 14], nbrs.first + x + 14);
                w[x + 15] = f.weight(slot.w_, current, base[x + 15], nbrs.first + x + 15);
            }

            if (para.sample_ == InverseTransformationSampling) {
                loop_count = nbrs.second - nbrs.first;
                for (auto x = 1; x < loop_count; ++x) {
                    w[x] += w[x - 1];
                }

                // If the sum of the transition probability is 0, then terminate the walk.
                if (w[loop_count - 1] == 0.0) {
                    slot.empty_ = true;
                    num_completed_walkers += 1;
                    num_terminated_walkers += 1;
                }
            }
            else if (para.sample_ == AliasSampling) {
                loop_count = nbrs.second - nbrs.first;
                if (!alias_initialization(w, base, slot.alias_,
                                          loop_count, large, small)) {
                    slot.empty_ = true;
                    num_completed_walkers += 1;
                    num_terminated_walkers += 1;
                }
            }
            else if (para.sample_ == RejectionSampling) {
                loop_count = nbrs.second - nbrs.first;
                slot.dr_ = *std::max_element(w, w + loop_count);

                // If the max value_ is 0, then terminate the walk.
                if (slot.dr_ == 0) {
                    slot.empty_ = true;
                    num_completed_walkers += 1;
                    num_terminated_walkers += 1;
                }
            }
        }
    }
}

template<class F>
void update(Graph *graph, BufferSlot *ring, int &num_completed_walkers, int &current_id, int num_walkers,
            WalkerMeta *walkers, F &f, uint64_t &step_count) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            // Update the status of the walker.
            step_count += 1;
            if (f.update(slot.w_, slot.prev_, slot.w_.current_, 0)) {
                // If the walker completes, then set the slot as empty.
                slot.empty_ = true;
                num_completed_walkers += 1;
#ifdef LOG_SEQUENCE
                if (slot.seq_ != nullptr) {
                    walkers[slot.local_id_].seq_ = slot.seq_;
                    walkers[slot.local_id_].length_ = slot.w_.length_;
                    slot.seq_ = slot.seq_ + slot.w_.length_;

                }
#endif
            }
        }

        // If the slot is empty, then add a new walker to the ring buffer.
        if (slot.empty_) {
            if (current_id < num_walkers) {
                slot.empty_ = false;
                slot.local_id_ = current_id;
                slot.w_ = walkers[current_id++];
                slot.w_.seq_ = slot.seq_;

                if (slot.seq_ != nullptr)
                    slot.seq_[0] = slot.w_.source_;
            }
        }
    }
}

template<class F> void *uniform_compute(void *ptr) {
    auto p = static_cast<ThreadParameter<F>*>(ptr);

    Graph *g = p->graph_;
    intT **s = p->seq_buffer_;
    int num_walkers = p->num_walkers_;
    WalkerMeta *q = p->queries_;
    F f = p->f_;
    int next = 0;
    int num_completed_walkers = 0;
    HyperParameter para = g_para;

    sfmt_t sfmt;
    BufferSlot r[RING_SIZE];

#ifdef ENABLE_INTERLEAVING
#ifdef TEST_AMAC
    AMAC_uniform_frame frames[SEARCH_RING_SIZE];
#endif
#endif
    assert(RING_SIZE < num_walkers);
    sfmt_init_gen_rand(&sfmt, p->id_);

    for (int i = 0; i < RING_SIZE; ++i) {
        r[i].empty_ = false;
        r[i].w_ = q[next];

        if (para.length_ > 0) {
            r[i].local_id_ = next;
            q[next].seq_ = r[i].w_.seq_ = r[i].seq_ = s[i];
            r[i].seq_[0] = r[i].w_.source_;
        }

        next += 1;
    }
    auto start = std::chrono::high_resolution_clock::now();

#ifdef COLLECT_PERF_COUNTER
    std::vector<int> evts;
    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
    std::vector<unsigned long long> creeper_results;
    creeper_results.resize(evts.size());
    unified.start();
#endif

    uint64_t step_count = 0;
    while (num_completed_walkers < num_walkers) {
#ifdef ENABLE_INTERLEAVING
#ifdef TEST_AMAC
        uniform_amac_move(g, r, frames, &sfmt, para.length_);
#else
        uniform_interleaving_move(g, r, &sfmt, para.length_);
#endif
#else
        uniform_move(g, r, &sfmt, para.length_);
#endif

        update(g, r, num_completed_walkers, next, num_walkers, q, f, step_count);
    }

    p->step_count_ = step_count;
    auto end = std::chrono::high_resolution_clock::now();

#ifdef COLLECT_PERF_COUNTER
    unified.end(creeper_results);
    log_info("%.6lf cycles per step, %.6lf instructions per step, %lu LLC access, %lu LLC misses",\
            creeper_results[0] / (double)step_count, creeper_results[1] / (double)step_count,
                    creeper_results[3], creeper_results[4]);
#endif

    double walking_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Thread %d Walking time (seconds): %.6lf", p->id_, walking_time);
    log_info("Thread %d Num of completed walkers: %d", p->id_, num_completed_walkers);
    log_info("Thread %d Num of steps: %zu", p->id_, step_count);
    log_info("Thread %d Throughput (steps per second): %.2lf", p->id_, step_count / walking_time);

    return nullptr;
}

template<class F> void *static_compute(void *ptr) {
    auto p = static_cast<ThreadParameter<F>*>(ptr);
    uint64_t step_count = 0;

     Graph *g = p->graph_;
     intT **s = p->seq_buffer_;
     WalkerMeta *q = p->queries_;
     int num_walkers = p->num_walkers_;
     F f = p->f_;
     int next = 0;
     int num_completed_walkers = 0;
     HyperParameter para = g_para;

     sfmt_t sfmt;
     BufferSlot r[RING_SIZE];

#ifdef ENABLE_INTERLEAVING
     ITSFrame<double> its_frames[SEARCH_RING_SIZE];
     RSFrame<double> rs_frames[SEARCH_RING_SIZE];

#ifdef TEST_AMAC
     AMAC_alias_frame amac_alias_frame[SEARCH_RING_SIZE];
     AMAC_rj_frame amac_rj_frames[SEARCH_RING_SIZE];
     AMAC_its_frame amac_its_frames[SEARCH_RING_SIZE];
#endif

#endif

     assert(RING_SIZE <= num_walkers);
     sfmt_init_gen_rand(&sfmt, p->id_);

     for (int i = 0; i < RING_SIZE; ++i) {
         r[i].empty_ = false;
         r[i].w_ = q[next];

         if (para.length_ > 0) {
             r[i].local_id_ = next;
             r[i].w_.seq_ = r[i].seq_ = s[i];
             r[i].seq_[0] = r[i].w_.source_;
         }

         next += 1;
     }

     /**
      * Start to compute.
      */
     auto start = std::chrono::high_resolution_clock::now();

#ifdef COLLECT_PERF_COUNTER
     std::vector<int> evts;
     evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
     evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
     evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
     evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
     evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
     LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
     std::vector<unsigned long long> creeper_results;
     creeper_results.resize(evts.size());

     unified.start();
#endif



     while (num_completed_walkers < num_walkers) {
         // Gather & move.
         if (para.sample_ == InverseTransformationSampling) {
#ifdef ENABLE_INTERLEAVING
#ifdef TEST_AMAC
             static_its_amac_move(g, r, amac_its_frames, &sfmt, para.length_);
#else
             static_its_interleaving_move(g, r, its_frames, &sfmt, para.length_);
#endif
#else
             static_its_move(g, r, &sfmt, para.length_);
#endif
         } else if (para.sample_ == AliasSampling) {
#ifdef ENABLE_INTERLEAVING
#ifdef TEST_AMAC
             static_alias_AMAC_move(g, r, &sfmt, amac_alias_frame, para.length_);
#else
             static_alias_interleaving_move(g, r, &sfmt, para.length_);
#endif
#else
             static_alias_move(g, r, &sfmt, para.length_);
#endif
         } else if (para.sample_ == RejectionSampling) {
#ifdef ENABLE_INTERLEAVING
#ifdef TEST_AMAC
             static_rs_amac_move(g, r, amac_rj_frames, &sfmt, para.length_);
#else
             static_rs_interleaving_move(g, r, rs_frames, &sfmt, para.length_);
#endif
#else
             static_rs_move(g, r, &sfmt, para.length_);
#endif
         }

         // Update.
         update(g, r, num_completed_walkers, next, num_walkers, q, f, step_count);
     }

     p->step_count_ = step_count;

    auto end = std::chrono::high_resolution_clock::now();

#ifdef COLLECT_PERF_COUNTER
    unified.end(creeper_results);
    log_info("%.6lf cycles per step, %.6lf instructions per step, %lu LLC access, %lu LLC misses",\
            creeper_results[0] / (double)step_count, creeper_results[1] / (double)step_count,
             creeper_results[3], creeper_results[4]);
#endif

    double walking_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Thread %d Walking time (seconds): %.6lf", p->id_, walking_time);
    log_info("Thread %d Num of completed walkers: %d", p->id_, num_completed_walkers);
    log_info("Thread %d Num of steps: %zu", p->id_, step_count);
    log_info("Thread %d Throughput (steps per second): %.2lf", p->id_, step_count / walking_time);

    return nullptr;
}

template<class F> void *dynamic_compute(void *ptr) {
    auto p = static_cast<ThreadParameter<F>*>(ptr);
    Graph* g = p->graph_;
    double* w = p->weight_buffer_;
    intT** s = p->seq_buffer_;
    AliasSlot* a = p->alias_buffer_;
    WalkerMeta* q = p->queries_;
    int num_walkers = p->num_walkers_;
    F f = p->f_;
    int next = 0;
    int num_completed_walkers = 0;
    int num_terminated_walkers = 0;
    HyperParameter para = g_para;

    sfmt_t sfmt;
    BufferSlot r[RING_SIZE];

#ifdef ENABLE_INTERLEAVING
    ITSFrame<double> its_frames[SEARCH_RING_SIZE];
    RSFrame<double> rs_frames[SEARCH_RING_SIZE];
    MaxRSFrame<intT> max_rs_frames[SEARCH_RING_SIZE];
#ifdef TEST_AMAC
    AMAC_max_rj_frame max_rs_amac_frames[SEARCH_RING_SIZE];
#endif
#endif

    assert(RING_SIZE < num_walkers);
    sfmt_init_gen_rand(&sfmt, p->id_);

    for (int i = 0; i < RING_SIZE; ++i) {
        r[i].empty_ = false;
        r[i].w_ = q[next];

        if (w != nullptr)
            r[i].weight_ = w + (int64_t)g->max_degree() * i;

        if (a != nullptr)
            r[i].alias_ = a + (int64_t)g->max_degree() * i;

        if (para.length_ > 0) {
            r[i].local_id_ = next;
            r[i].w_.seq_ = r[i].seq_ = s[i];
            r[i].seq_[0] = r[i].w_.source_;
        }

        next += 1;
    }

    std::queue<intT> large;
    std::queue<intT> small;

    /**
     * Start to compute.
     */
    auto start = std::chrono::high_resolution_clock::now();

#ifdef COLLECT_PERF_COUNTER
    std::vector<int> evts;
    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);

    std::vector<unsigned long long> creeper_results;
    creeper_results.resize(evts.size());
    unified.start();
#endif

    uint64_t step_count = 0;

    while (num_completed_walkers < num_walkers) {
        if (para.sample_ != MaxWeightRejectionSampling) {
            // Gather & move.
            dynamic_gather(g, r, f, para, num_completed_walkers, num_terminated_walkers, large,
                           small);

            if (para.sample_ == InverseTransformationSampling) {
                // dynamic_its_initialization(g, r);

#ifdef ENABLE_INTERLEAVING
                dynamic_its_interleaving_move(g, r, its_frames, &sfmt, para.length_);
#else
                dynamic_its_move(g, r, &sfmt, para.length_);
#endif
            } else if (para.sample_ == AliasSampling) {
                // dynamic_alias_initialization(g, r, large, small);

#ifdef ENABLE_INTERLEAVING
                dynamic_alias_interleaving_move(g, r, &sfmt, para.length_);
#else
                dynamic_alias_move(g, r, &sfmt, para.length_);
#endif
            } else if (para.sample_ == RejectionSampling) {
                // dynamic_rs_initialization(g, r);

#ifdef ENABLE_INTERLEAVING
                dynamic_rs_interleaving_move(g, r, rs_frames, &sfmt, para.length_);
#else
                dynamic_rs_move(g, r, &sfmt, para.length_);
#endif
            }
        }
        else {
            for (int i = 0; i < RING_SIZE; ++i) {
                BufferSlot &slot = r[i];
                if (!slot.empty_) {
                    slot.dr_ = f.max_weight(slot.w_);
                }
            }

#ifdef ENABLE_INTERLEAVING
#ifdef TEST_AMAC
            dynamic_max_weight_rs_amac_move(g, r, max_rs_amac_frames, &sfmt, para.length_, f);
#else
            dynamic_max_weight_rs_interleaving_move(g, r, max_rs_frames, &sfmt, para.length_, f);
#endif
#else
            dynamic_max_weight_rs_move(g, r, &sfmt, g_para.length_, f);
#endif
        }

        // Update.
        update(g, r, num_completed_walkers, next, num_walkers, q, f, step_count);
    }

    p->step_count_ = step_count;
    auto end = std::chrono::high_resolution_clock::now();

#ifdef COLLECT_PERF_COUNTER
    unified.end(creeper_results);
    log_info("%.6lf cycles per step, %.6lf instructions per step, %lu LLC access, %lu LLC misses",\
            creeper_results[0] / (double)step_count, creeper_results[1] / (double)step_count,
             creeper_results[3], creeper_results[4]);
#endif

    double walking_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Thread %d Walking time (seconds): %.6lf", p->id_, walking_time);
    log_info("Thread %d Num of completed walkers: %d", p->id_, num_completed_walkers);
    log_info("Thread %d Num of terminated walkers: %d", p->id_, num_terminated_walkers);
    log_info("Thread %d Num of steps: %zu", p->id_, step_count);
    log_info("Thread %d Throughput (steps per second): %.2lf", p->id_, step_count / walking_time);
    return nullptr;
}

template<class F> void compute(Graph& graph, std::vector<WalkerMeta>& walkers, F f) {
    /**
     * Output the configuration:
     * 1. Execution mode
     * 2. Sampling method
     * 3. Thread number
     * 4. Length
     */
    log_info("Configuration: %s, %s, thread num %d, length %d, cache hint %d",
            execution[g_para.execution_], sampling[g_para.sample_], g_num_threads, g_para.length_, PREFETCH_HINT);

    if (g_para.execution_ == Static) {
        static_initialization(&graph, f);
    }

    /**
     * Initialize the resource...
     */
    int num_threads = g_num_threads;

    double** weight_buffer = nullptr;
    AliasSlot** alias_slot_buffer = nullptr;
    intT*** seq_buffer = nullptr;

    if (g_para.execution_ == Dynamic) {
        weight_buffer = new double*[num_threads];

        if (g_para.sample_ == AliasSampling) {
            alias_slot_buffer = new AliasSlot*[num_threads];
        }

        for (int i = 0; i < num_threads; ++i) {
            weight_buffer[i] = new double[graph.max_degree() * RING_SIZE];

            if (g_para.sample_ == AliasSampling) {
                alias_slot_buffer[i] = new AliasSlot[graph.max_degree() * RING_SIZE];
            }
        }
    }

    std::pair<WalkerMeta*, int> tasks[num_threads];
    int length = walkers.size() / num_threads;
    for (int i = 0; i < num_threads; ++i) {
        tasks[i].first = walkers.data() + i * length;
        tasks[i].second = i == (num_threads - 1) ? walkers.size() - i * length : length;
    }

    /**
    * The buffer maintaining the output. Currently, ThunderRW can only handle the scenario that the memory can
    * maintain all the results. This can be optimized by the double buffer strategy to output the results to disk.
    */
    if (g_para.length_ > 0) {
        seq_buffer = new intT**[num_threads];
        for (int i = 0; i < num_threads; ++i) {

            seq_buffer[i] = new intT*[RING_SIZE];
            for (int j = 0; j < RING_SIZE; ++j) {
#ifdef LOG_SEQUENCE
                int local_walk_num = tasks[i].second;
                // 1.15 as the buffer for the walker with variant length.
                auto per_slot_buffer_size = static_cast<size_t>(local_walk_num * 1.15 / RING_SIZE * g_para.length_);
                seq_buffer[i][j] = new intT[per_slot_buffer_size];
#else
                seq_buffer[i][j] = new intT[g_para.length_];
#endif
            }
        }
    }

    log_info("Compute...");

    ThreadParameter<F> parameters[num_threads];
    pthread_t threads[num_threads];
    pthread_attr_t attr;
    cpu_set_t cpus;
    pthread_attr_init(&attr);

    auto start = std::chrono::high_resolution_clock::now();



    for (int i = 0; i < num_threads; ++i) {
        parameters[i] = {i, &graph, weight_buffer == nullptr ? nullptr : weight_buffer[i],
                         seq_buffer == nullptr ? nullptr : seq_buffer[i], alias_slot_buffer == nullptr ? nullptr : alias_slot_buffer[i],
                         tasks[i].first, tasks[i].second, 0, f};

        // set thread affinity.
        CPU_ZERO(&cpus);
        CPU_SET(i, &cpus);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);

        if (g_para.execution_ == Uniform) {
            pthread_create(&threads[i], &attr, uniform_compute<F>, &parameters[i]);
        }
        else if (g_para.execution_ == Static) {
            pthread_create(&threads[i], &attr, static_compute<F>, &parameters[i]);
        }
        else {
            pthread_create(&threads[i], &attr, dynamic_compute<F>, &parameters[i]);
        }
    }


    uint64_t step_count = 0;
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
        step_count += parameters[i].step_count_;
    }

    auto end = std::chrono::high_resolution_clock::now();

    double walking_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
    log_info("Walking time (seconds): %.6lf", walking_time);
    log_info("Num of steps: %zu", step_count);
    log_info("Throughput (steps per second): %.2lf", step_count / walking_time);

    log_info("Release...");

    if (g_para.execution_ == Dynamic && weight_buffer != nullptr) {
        for (int i = 0; i < num_threads; ++i) {
            delete[] weight_buffer[i];

            if (g_para.sample_ == AliasSampling) {
                delete[] alias_slot_buffer[i];
            }
        }
        delete[] weight_buffer;

        if (g_para.sample_ == AliasSampling) {
            delete[] alias_slot_buffer;
        }
    }

    if (g_para.length_ > 0 && seq_buffer != nullptr) {
        for (int i = 0; i < num_threads; ++i) {
            for (int j = 0; j < RING_SIZE; ++j) {
                delete[] seq_buffer[i][j];
            }
            delete[] seq_buffer[i];
        }
        delete[] seq_buffer;
    }
}

#endif //XTRAGRAPHCOMPUTING_CREEPER_H
