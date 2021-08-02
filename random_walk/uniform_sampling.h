//
// Created by Shixuan Sun on 11/28/20.
//

#ifndef XTRAGRAPHCOMPUTING_UNIFORM_SAMPLING_H
#define XTRAGRAPHCOMPUTING_UNIFORM_SAMPLING_H

#include "types.h"
#include "amac_frame.h"

void uniform_interleaving_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    // Stage 1: generate random number & prefetch the degree.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            slot.r_ = sfmt_genrand_uint32(sfmt);
            _mm_prefetch((void*)(graph->offset_pair_ + slot.w_.current_), PREFETCH_HINT);
        }
    }

    // Stage 2: generate the position & prefetch the neighbor.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.offset_ = graph->offset_pair_[slot.w_.current_];
            slot.r_ = slot.offset_.first + (slot.r_ % (slot.offset_.second - slot.offset_.first));
            _mm_prefetch((void*)(graph->adj_ + slot.r_), PREFETCH_HINT);
        }
    }

    // Stage 3: update the walker.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.w_.current_ = graph->adj_[slot.r_];

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }

            slot.w_.length_ += 1;
        }
    }
}

void uniform_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto neighbors = graph->neighbors(slot.w_.current_);
            auto random_value = sfmt_genrand_uint32(sfmt);
            auto selected_position = random_value % neighbors.second;
            slot.w_.current_ = neighbors.first[selected_position];

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

void uniform_amac_move(Graph *graph, BufferSlot *ring, AMAC_uniform_frame* frames, sfmt_t *sfmt, int length) {
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == AMAC_uniform_frame::Empty) {
                // Execute stage S0
                fr->init(graph->offset_pair_, graph->adj_,
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
                    if (fr->state == AMAC_uniform_frame::S0) {
                        fr->execute_S0();
                    }
                    else if (fr->state == AMAC_uniform_frame::S1) {
                        fr->execute_S1(sfmt);
                    }
                    else if (fr->state == AMAC_uniform_frame::S2) {
                        fr->execute_S2();

                        // Find a result.
                        auto& t_slot = ring[fr->id];
                        t_slot.prev_ = t_slot.w_.current_;
                        t_slot.w_.current_ = fr->value_;

                        if (t_slot.w_.length_ < length) {
                            t_slot.seq_[t_slot.w_.length_] = t_slot.w_.current_;
                        }
                        t_slot.w_.length_ += 1;

                        // initialize the slot.
                        fr->init(graph->offset_pair_, graph->adj_,
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
            auto& frame = frames[i];
            if (frame.state == AMAC_uniform_frame::S0) {
                more_work = true;
                frame.execute_S0();
            }
            else if (frame.state == AMAC_uniform_frame::S1) {
                more_work = true;
                frame.execute_S1(sfmt);
            }
            else if (frame.state == AMAC_uniform_frame::S2) {
                frame.execute_S2();

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
#endif //XTRAGRAPHCOMPUTING_UNIFORM_SAMPLING_H
