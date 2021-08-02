//
// Created by Shixuan Sun on 11/28/20.
//

#ifndef XTRAGRAPHCOMPUTING_INVERSE_TRANSFORMATION_SAMPLING_H
#define XTRAGRAPHCOMPUTING_INVERSE_TRANSFORMATION_SAMPLING_H

#include "types.h"
#include "its_frame.h"
#include "amac_frame.h"


template<typename T> void its_in_place_initialization(T* src, int length) {
    for (auto i = 1; i < length; ++i) {
        src[i] += src[i - 1];
    }
}

template<typename T> void its_initialization(const T* src, T* dst, int length) {
    dst[0] = src[0];
    for (auto i = 1; i < length; ++i) {
        dst[i] = dst[i - 1] + src[i];
    }
}

void dynamic_its_initialization(Graph *graph, BufferSlot *ring) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& offset = graph->offset_pair_[slot.w_.current_];
            its_in_place_initialization(slot.weight_, offset.second - offset.first);
        }
    }
}

void static_its_interleaving_move(Graph *graph, BufferSlot *ring, ITSFrame<double> *frames, sfmt_t *sfmt, int length) {
    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        // Stage 1: generate random number & prefetch the degree.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *)(graph->offset_pair_ + slot.w_.current_), _MM_HINT_T0);
                slot.prev_ = slot.w_.current_;
                slot.dr_ = sfmt_genrand_real2(sfmt);
            }
        }

        // Stage 2.1: get the neighbor, weight offsets and the degree value_ & prefetch the max weight.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.offset_ = graph->offset_pair_[slot.w_.current_];
                _mm_prefetch((void *)(graph->edge_weight_prefix_sum_ + slot.offset_.second - 1), PREFETCH_HINT);
            }
        }

        // Stage 2.2: generate the target weight & initialize the search.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot& slot = ring[i];
            if (!slot.empty_) {
                slot.dr_ = graph->edge_weight_prefix_sum_[slot.offset_.second - 1] * slot.dr_;
            }
        }
    }

    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == ITSFrame<double>::EMPTY) {
                // Add a task.
                fr->init(graph->edge_weight_prefix_sum_, slot.offset_.first, slot.offset_.second, slot.dr_, i);
                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->run()) {
                        // Find a result.
                        ring[fr->id].r_ = ring[fr->id].offset_.first + fr->result;

                        // Add a new task.
                        fr->init(graph->edge_weight_prefix_sum_, slot.offset_.first, slot.offset_.second, slot.dr_, i);

                        if (search_ring_id == SEARCH_RING_SIZE - 1)
                            search_ring_id = 0;
                        else
                            search_ring_id += 1;

                        break;
                    }

                    if (search_ring_id == SEARCH_RING_SIZE - 1)
                        search_ring_id = 0;
                    else
                        search_ring_id += 1;

                    fr = &frames[search_ring_id];
                }
            }
        }
    }

    // Stage 2.3: start the search.
    bool more_work;
    do {
        more_work = false;
        for (int i = 0; i < SEARCH_RING_SIZE; ++i) {
            ITSFrame<double>& frame = frames[i];
            if (frame.state == ITSFrame<double>::RUNNING) {
                more_work = true;
                if(frame.run()) {
                    // Find the result.
                    ring[frame.id].r_ = ring[frame.id].offset_.first + frame.result;
                    frame.state = ITSFrame<double>::EMPTY;
                }
            }
        }
    } while (more_work);

    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *) (graph->adj_ + slot.r_), PREFETCH_HINT);
            }
        }

        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.w_.current_ = graph->adj_[slot.r_];
                if (slot.w_.length_ < length) {
                    slot.seq_[slot.w_.length_] = slot.w_.current_;
                }
                slot.w_.length_ += 1;
            }
        }
    }
}

void dynamic_its_interleaving_move(Graph *graph, BufferSlot *ring, ITSFrame<double> *frames, sfmt_t *sfmt, int length) {
    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        // Stage 1: generate random number & prefetch the degree.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *)(graph->offset_pair_ + slot.w_.current_), PREFETCH_HINT);
                slot.prev_ = slot.w_.current_;
                slot.dr_ = sfmt_genrand_real2(sfmt);
            }
        }

        // Stage 2.1: get the neighbor, weight offsets and the degree value_ & prefetch the max weight.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.offset_ = graph->offset_pair_[slot.w_.current_];
                slot.offset_.second = slot.offset_.second - slot.offset_.first;
                _mm_prefetch((void *)(slot.weight_ + slot.offset_.second - 1), PREFETCH_HINT);
            }
        }

        // Stage 2.2: generate the target weight & initialize the search.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot& slot = ring[i];
            if (!slot.empty_) {
                slot.dr_ = slot.weight_[slot.offset_.second - 1] * slot.dr_;
            }
        }
    }

    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == ITSFrame<double>::EMPTY) {
                // Add a task.
                fr->init(slot.weight_, 0, slot.offset_.second, slot.dr_, i);
                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->run()) {
                        // Find a result.
                        ring[fr->id].r_ = ring[fr->id].offset_.first + fr->result;

                        // Add a new task.
                        fr->init(slot.weight_, 0, slot.offset_.second, slot.dr_, i);

                        if (search_ring_id == SEARCH_RING_SIZE - 1)
                            search_ring_id = 0;
                        else
                            search_ring_id += 1;

                        break;
                    }

                    if (search_ring_id == SEARCH_RING_SIZE - 1)
                        search_ring_id = 0;
                    else
                        search_ring_id += 1;

                    fr = &frames[search_ring_id];
                }
            }
        }
    }

    // Stage 2.3: start the search.
    bool more_work;
    do {
        more_work = false;
        for (int i = 0; i < SEARCH_RING_SIZE; ++i) {
            ITSFrame<double>& frame = frames[i];
            if (frame.state == ITSFrame<double>::RUNNING) {
                more_work = true;
                if(frame.run()) {
                    // Find the result.
                    ring[frame.id].r_ = ring[frame.id].offset_.first + frame.result;
                    frame.state = ITSFrame<double>::EMPTY;
                }
            }
        }
    } while (more_work);

    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *) (graph->adj_ + slot.r_), PREFETCH_HINT);
            }
        }

        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.w_.current_ = graph->adj_[slot.r_];
                if (slot.w_.length_ < length) {
                    slot.seq_[slot.w_.length_] = slot.w_.current_;
                }
                slot.w_.length_ += 1;
            }
        }
    }
}

void static_its_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto neighbor_base = graph->adj_ + graph->offset_[slot.w_.current_];
            auto weight_base = graph->edge_weight_prefix_sum_ + graph->offset_[slot.w_.current_];
            auto degree = graph->offset_[slot.w_.current_ + 1] - graph->offset_[slot.w_.current_];

            // Generate a random number.
            double max_weight = weight_base[degree - 1];
            auto target_weight = sfmt_genrand_real2(sfmt) * max_weight;

            // Find the position in the weighted array.
            // int selected_position = std::upper_bound(weight_base, weight_base + degree, target_weight) - weight_base;

            int selected_position = Utility::upper_bound_branchless(weight_base, 0, degree, target_weight);

            // Find the neighbor and update.
            slot.w_.current_ = neighbor_base[selected_position];

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

void dynamic_its_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto neighbors = graph->neighbors(slot.w_.current_);
            int max_weight_position = neighbors.second - 1;
            double max_weight = slot.weight_[max_weight_position];

            // Generate a random number.
            auto random_value = sfmt_genrand_real2(sfmt);
            auto target_weight = random_value * max_weight;

            // Find the position in the weighted array.
            // int selected_position = std::upper_bound(slot.weight_, slot.weight_ + neighbors.second, target_weight) - slot.weight_;
            int selected_position = Utility::upper_bound_branchless(slot.weight_, 0, neighbors.second, target_weight);

            // Find the neighbor and update.
            slot.w_.current_ = neighbors.first[selected_position];
            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

void static_its_amac_move(Graph *graph, BufferSlot *ring, AMAC_its_frame *frames, sfmt_t *sfmt, int length) {
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == AMAC_its_frame::Empty) {
                // Init
                fr->init(graph->adj_, graph->offset_pair_, graph->edge_weight_prefix_sum_,
                         slot.w_.current_, i);

                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->state == AMAC_its_frame::S0) {
                        fr->execute_S0();
                    }
                    else if (fr->state == AMAC_its_frame::S1) {
                        fr->execute_S1();
                    }
                    else if (fr->state == AMAC_its_frame::S2) {
                        fr->execute_S2(sfmt);
                    }
                    else if (fr->state == AMAC_its_frame::S3) {
                        fr->execute_S3();
                    }
                    else if (fr->state == AMAC_its_frame::S4) {
                        fr->execute_S4();
                    }
                    else if (fr->state == AMAC_its_frame::S5) {
                        fr->execute_S5();
                    }
                    else if (fr->state == AMAC_its_frame::S6) {
                        fr->execute_S6();
                        // Find a result.
                        auto& t_slot = ring[fr->id];
                        t_slot.prev_ = t_slot.w_.current_;
                        t_slot.w_.current_ = fr->value_;

                        if (t_slot.w_.length_ < length) {
                            t_slot.seq_[t_slot.w_.length_] = t_slot.w_.current_;
                        }
                        t_slot.w_.length_ += 1;

                        // Init
                        fr->init(graph->adj_, graph->offset_pair_, graph->edge_weight_prefix_sum_,
                                 slot.w_.current_, i);
                        break;
                    }


                    if (search_ring_id == SEARCH_RING_SIZE - 1)
                        search_ring_id = 0;
                    else
                        search_ring_id += 1;

                    fr = &frames[search_ring_id];
                }
            }
        }
    }

    // Stage 2.3: start the search.
    bool more_work;
    do {
        more_work = false;
        for (int i = 0; i < SEARCH_RING_SIZE; ++i) {
            AMAC_its_frame& frame = frames[i];
            if (frame.state == AMAC_its_frame::S0) {
                more_work = true;
                frame.execute_S0();
            }
            else if (frame.state == AMAC_its_frame::S1) {
                more_work = true;
                frame.execute_S1();
            }
            else if (frame.state == AMAC_its_frame::S2) {
                more_work = true;
                frame.execute_S2(sfmt);
            }
            else if (frame.state == AMAC_its_frame::S3) {
                more_work = true;
                frame.execute_S3();
            }
            else if (frame.state == AMAC_its_frame::S4) {
                more_work = true;
                frame.execute_S4();
            }
            else if (frame.state == AMAC_its_frame::S5) {
                more_work = true;
                frame.execute_S5();
            }
            else if (frame.state == AMAC_its_frame::S6) {
                frame.execute_S6();
                // Find a result.
                auto& t_slot = ring[frame.id];
                t_slot.prev_ = t_slot.w_.current_;
                t_slot.w_.current_ = frame.value_;

                if (t_slot.w_.length_ < length) {
                    t_slot.seq_[t_slot.w_.length_] = t_slot.w_.current_;
                }
                t_slot.w_.length_ += 1;
            }
        }
    } while (more_work);
}
#endif //XTRAGRAPHCOMPUTING_INVERSE_TRANSFORMATION_SAMPLING_H
