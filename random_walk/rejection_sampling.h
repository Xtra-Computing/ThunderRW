//
// Created by Shixuan Sun on 11/28/20.
//

#ifndef XTRAGRAPHCOMPUTING_REJECTION_SAMPLING_H
#define XTRAGRAPHCOMPUTING_REJECTION_SAMPLING_H

#include "types.h"
#include "rs_frame.h"
#include "amac_frame.h"

template<typename T> void rs_initialization(T* src, int length, T& result) {
    result = *std::max_element(src, src + length);
}

void dynamic_rs_initialization(Graph *graph, BufferSlot *ring) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& offset = graph->offset_pair_[slot.w_.current_];
            double temp;
            rs_initialization(slot.weight_, offset.second - offset.first, temp);
            slot.dr_ = temp;
        }
    }
}

void static_rs_interleaving_move(Graph *graph, BufferSlot *ring, RSFrame<double> *frames,
        sfmt_t *sfmt, int length) {
    // Stage 1: Prefetch the degree.
    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *) (graph->offset_pair_ + slot.w_.current_), _MM_HINT_T0);
            }
        }

        // Stage 2: Set the offset & Prefetch the max weight.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.offset_ = graph->offset_pair_[slot.w_.current_];
                _mm_prefetch((void *) (graph->edge_weight_rejection_max_ + slot.w_.current_), PREFETCH_HINT);
            }
        }

        // Stage 3: Set the max weight.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.dr_ = graph->edge_weight_rejection_max_[slot.w_.current_];
            }
        }
    }


    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == RSFrame<double>::EMPTY) {
                // Add a task.
                fr->init(graph->edge_weight_rejection_, slot.offset_.first, slot.offset_.second, slot.dr_, i, sfmt);
                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->run(sfmt)) {
                        // Find a result.
                        ring[fr->id].r_ = ring[fr->id].offset_.first + fr->result;

                        // Add a new task.
                        fr->init(graph->edge_weight_rejection_, slot.offset_.first, slot.offset_.second, slot.dr_, i, sfmt);

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
            RSFrame<double>& frame = frames[i];
            if (frame.state == RSFrame<double>::RUNNING) {
                more_work = true;
                if(frame.run(sfmt)) {
                    // Find the result.
                    ring[frame.id].r_ = ring[frame.id].offset_.first + frame.result;
                    frame.state = RSFrame<double>::EMPTY;
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

template<typename F> void dynamic_max_weight_rs_interleaving_move(Graph *graph, BufferSlot *ring, MaxRSFrame<intT> *frames,
                                  sfmt_t *sfmt, int length, F& f) {
    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        // Stage 1: Prefetch the degree.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *) (graph->offset_pair_ + slot.w_.current_), PREFETCH_HINT);
            }
        }

        // Stage 2: Set the offset & Prefetch the max weight.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.offset_ = graph->offset_pair_[slot.w_.current_];
            }
        }
    }
    // Stage 3: search
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == MaxRSFrame<intT>::EMPTY) {
                // Add a task.
                fr->init(graph->adj_, slot.offset_.first, slot.offset_.second, slot.dr_, i, sfmt);
                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->run(sfmt, f, ring[fr->id].w_, ring[fr->id].offset_.first)) {
                        // Find a result.
                        ring[fr->id].r_ = fr->result;

                        // Add a new task.
                        fr->init(graph->adj_, slot.offset_.first, slot.offset_.second, slot.dr_, i, sfmt);
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
            auto& frame = frames[i];
            if (frame.state == MaxRSFrame<intT>::RUNNING) {
                more_work = true;
                if(frame.run(sfmt, f, ring[frame.id].w_, ring[frame.id].offset_.first)) {
                    // Find the result.
                    ring[frame.id].r_ = frame.result;
                    frame.state = MaxRSFrame<intT>::EMPTY;
                }
            }
        }
    } while (more_work);



    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot &slot = ring[i];
        if (!slot.empty_) {
            slot.w_.current_ = slot.r_;
            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

void dynamic_rs_interleaving_move(Graph *graph, BufferSlot *ring, RSFrame<double> *frames,
                                 sfmt_t *sfmt, int length) {
    for (int j = 0; j < RING_SIZE; j += SMALL_RING_SIZE) {
        // Stage 1: Prefetch the degree.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                _mm_prefetch((void *) (graph->offset_pair_ + slot.w_.current_), PREFETCH_HINT);
            }
        }

        // Stage 2: Set the offset & Prefetch the max weight.
        for (int i = j; i < j + SMALL_RING_SIZE; ++i) {
            BufferSlot &slot = ring[i];
            if (!slot.empty_) {
                slot.offset_ = graph->offset_pair_[slot.w_.current_];
            }
        }
    }
    // Stage 3: search
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == RSFrame<double>::EMPTY) {
                // Add a task.
                fr->init(slot.weight_, 0, slot.offset_.second - slot.offset_.first, slot.dr_, i, sfmt);
                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->run(sfmt)) {
                        // Find a result.
                        ring[fr->id].r_ = ring[fr->id].offset_.first + fr->result;

                        // Add a new task.
                        fr->init(slot.weight_, 0, slot.offset_.second - slot.offset_.first,
                                 slot.dr_, i, sfmt);

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
            RSFrame<double>& frame = frames[i];
            if (frame.state == RSFrame<double>::RUNNING) {
                more_work = true;
                if(frame.run(sfmt)) {
                    // Find the result.
                    ring[frame.id].r_ = ring[frame.id].offset_.first + frame.result;
                    frame.state = RSFrame<double>::EMPTY;
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

void static_rs_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto degree = graph->offset_[slot.w_.current_ + 1] - graph->offset_[slot.w_.current_];
            auto edge_weight = graph->edge_weight_rejection_ + graph->offset_[slot.w_.current_];
            auto neighbors = graph->adj_ + graph->offset_[slot.w_.current_];
            auto max_weight = graph->edge_weight_rejection_max_[slot.w_.current_];

            // Select a neighbor.

            double r;
            uint32_t p0;
            do {
                // Generate x position.
                auto r0 = sfmt_genrand_uint32(sfmt);
                p0 = r0 % degree;

                // Generate y position.
                r = sfmt_genrand_real2(sfmt) * max_weight;
            } while (r > edge_weight[p0]);

            // Find the neighbor and update.
            slot.w_.current_ = neighbors[p0];

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

template<typename F> void dynamic_max_weight_rs_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length, F& f) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto offset = graph->offset_pair_[slot.w_.current_];
            auto degree = offset.second - offset.first;
            auto neighbors = graph->adj_ + offset.first;
            auto max_weight = slot.dr_;

            // Select a neighbor.
            double r;
            double selected_r;
            uint32_t p0;
            do {
                // Generate x position.
                auto r0 = sfmt_genrand_uint32(sfmt);
                p0 = r0 % degree;
                // Get the value_.
                selected_r = f.weight(slot.w_, slot.w_.current_, neighbors[p0], offset.first + p0);
                // Generate y position.
                r = sfmt_genrand_real2(sfmt) * max_weight;
            } while (r > selected_r);


            // Find the neighbor and update.
            slot.w_.current_ = neighbors[p0];

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

void dynamic_rs_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto degree = graph->offset_[slot.w_.current_ + 1] - graph->offset_[slot.w_.current_];
            auto edge_weight = slot.weight_;
            auto neighbors = graph->adj_ + graph->offset_[slot.w_.current_];
            auto max_weight = slot.dr_;

            // Select a neighbor.
            double r;
            uint32_t p0;
            do {
                // Generate x position.
                auto r0 = sfmt_genrand_uint32(sfmt);
                p0 = r0 % degree;

                // Generate y position.
                r = sfmt_genrand_real2(sfmt) * max_weight;
            } while (r > edge_weight[p0]);


            // Find the neighbor and update.
            slot.w_.current_ = neighbors[p0];

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}


/**
 * Implement the amac rejection sampling for test purpose
 */
void static_rs_amac_move(Graph *graph, BufferSlot *ring, AMAC_rj_frame *frames,
                                 sfmt_t *sfmt, int length) {
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == AMAC_rj_frame::Empty) {
                // Init
                fr->init(graph->adj_, graph->offset_pair_, graph->edge_weight_rejection_max_,
                        graph->edge_weight_, slot.w_.current_, i);

                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->state == AMAC_rj_frame::S0) {
                        fr->execute_S0();
                    }
                    else if (fr->state == AMAC_rj_frame::S1) {
                        fr->execute_S1();
                    }
                    else if (fr->state == AMAC_rj_frame::S2) {
                        fr->execute_S2();
                    }
                    else if (fr->state == AMAC_rj_frame::S3) {
                        fr->execute_S3(sfmt);
                    }
                    else if (fr->state == AMAC_rj_frame::S4) {
                        fr->execute_S4();
                    }
                    else if (fr->state == AMAC_rj_frame::S5) {
                        fr->execute_S5();
                    }
                    else if (fr->state == AMAC_rj_frame::S6) {
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
                        fr->init(graph->adj_, graph->offset_pair_, graph->edge_weight_rejection_max_,
                                 graph->edge_weight_, slot.w_.current_, i);
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
            AMAC_rj_frame& frame = frames[i];
            if (frame.state == AMAC_rj_frame::S0) {
                more_work = true;
                frame.execute_S0();
            }
            else if (frame.state == AMAC_rj_frame::S1) {
                more_work = true;
                frame.execute_S1();
            }
            else if (frame.state == AMAC_rj_frame::S2) {
                more_work = true;
                frame.execute_S2();
            }
            else if (frame.state == AMAC_rj_frame::S3) {
                more_work = true;
                frame.execute_S3(sfmt);
            }
            else if (frame.state == AMAC_rj_frame::S4) {
                more_work = true;
                frame.execute_S4();
            }
            else if (frame.state == AMAC_rj_frame::S5) {
                more_work = true;
                frame.execute_S5();
            }
            else if (frame.state == AMAC_rj_frame::S6) {
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

template<typename F> void dynamic_max_weight_rs_amac_move(Graph *graph, BufferSlot *ring, AMAC_max_rj_frame *frames,
                                                                  sfmt_t *sfmt, int length, F& f) {
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == AMAC_max_rj_frame::Empty) {
                // Init
                fr->init(graph->adj_, graph->offset_pair_, slot.w_.current_, slot.dr_, i);

                if (search_ring_id == SEARCH_RING_SIZE - 1) {
                    search_ring_id = 0;
                }
                else {
                    search_ring_id += 1;
                }
            }
            else {
                for (;;) {
                    if (fr->state == AMAC_max_rj_frame::S0) {
                        fr->execute_S0();
                    }
                    else if (fr->state == AMAC_max_rj_frame::S1) {
                        fr->execute_S1();
                    }
                    else if (fr->state == AMAC_max_rj_frame::S2) {
                        fr->execute_S2(sfmt);
                    }
                    else if (fr->state == AMAC_max_rj_frame::S3) {
                        fr->execute_S3(f, ring[fr->id].w_);
                    }
                    else if (fr->state == AMAC_max_rj_frame::S4) {
                        fr->execute_S4();

                        // Find a result.
                        auto& t_slot = ring[fr->id];
                        t_slot.prev_ = t_slot.w_.current_;
                        t_slot.w_.current_ = fr->value_;

                        if (t_slot.w_.length_ < length) {
                            t_slot.seq_[t_slot.w_.length_] = t_slot.w_.current_;
                        }
                        t_slot.w_.length_ += 1;

                        // Init
                        fr->init(graph->adj_, graph->offset_pair_, slot.w_.current_, slot.dr_, i);
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
            auto& frame = frames[i];
            if (frame.state == AMAC_max_rj_frame::S0) {
                more_work = true;
                frame.execute_S0();
            }
            else if (frame.state == AMAC_max_rj_frame::S1) {
                more_work = true;
                frame.execute_S1();
            }
            else if (frame.state == AMAC_max_rj_frame::S2) {
                more_work = true;
                frame.execute_S2(sfmt);
            }
            else if (frame.state == AMAC_max_rj_frame::S3) {
                more_work = true;
                frame.execute_S3(f, ring[frame.id].w_);
            }
            else if (frame.state == AMAC_max_rj_frame::S4) {
                frame.execute_S4();

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
#endif //XTRAGRAPHCOMPUTING_REJECTION_SAMPLING_H
