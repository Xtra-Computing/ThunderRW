//
// Created by Shixuan Sun on 11/28/20.
//

#ifndef XTRAGRAPHCOMPUTING_ALIAS_SAMPLING_H
#define XTRAGRAPHCOMPUTING_ALIAS_SAMPLING_H

#include "types.h"

#include "amac_frame.h"

template<typename T> bool alias_initialization(T* src, intT* nbrs, AliasSlot* alias_table, int length,
                                               std::queue<intT>& large, std::queue<intT>& small) {
    int small_idx, large_idx, index = 0;

    auto weight_avg = std::accumulate(src, src + length, 0.0) / length;

    if (weight_avg == 0.0)
        return false;

    for (int j = 0; j < length; j++) {
        src[j] /= weight_avg;
        if (src[j] < 1)
            small.push(j);
        else
            large.push(j);
    }

    while ((!small.empty()) && (!large.empty())) {
        small_idx = small.front();
        large_idx = large.front();
        alias_table[index] = {src[small_idx], nbrs[small_idx], nbrs[large_idx]};
        index++;
        small.pop();
        src[large_idx] = src[large_idx] + src[small_idx] - 1;
        if (src[large_idx] < 1) {
            large.pop();
            small.push(large_idx);
        }
    }
    while (!small.empty()) {
        small_idx = small.front();
        small.pop();
        alias_table[index] = {1.0, nbrs[small_idx], nbrs[small_idx]};
        index++;
    }
    while (!large.empty()) {
        large_idx = large.front();
        large.pop();
        alias_table[index] = {1.0, nbrs[large_idx], nbrs[large_idx]};
        index++;
    }
    return true;
}

void dynamic_alias_initialization(Graph *graph, BufferSlot *ring, std::queue<intT>& large, std::queue<intT>& small) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& offset = graph->offset_pair_[slot.w_.current_];
            if (!alias_initialization(slot.weight_, graph->adj_ + offset.first, slot.alias_,
                                 offset.second - offset.first, large, small)) {
                slot.empty_ = true;
            }
        }
    }
}

void static_alias_interleaving_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    // Stage 1: generate the first random number & prefetch the degree.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            slot.r_ = sfmt_genrand_uint32(sfmt);
            _mm_prefetch((void*)(graph->offset_pair_ + slot.w_.current_), PREFETCH_HINT);
        }
    }

    // Stage 2: generate the position & prefetch the alias slot.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& offset = graph->offset_pair_[slot.w_.current_];
            slot.r_ = offset.first + (slot.r_ % (offset.second - offset.first));
            slot.dr_ = sfmt_genrand_real2(sfmt);
            _mm_prefetch((void*)(graph->edge_weight_alias_table_ + slot.r_), PREFETCH_HINT);
        }
    }

    // Stage 3: update the walker.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& alias_slot = graph->edge_weight_alias_table_[slot.r_];
            slot.w_.current_ = slot.dr_ <= alias_slot.alias_value_ ?
                    alias_slot.first_ : alias_slot.second_;

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }

            slot.w_.length_ += 1;
        }
    }
}

void dynamic_alias_interleaving_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    // Stage 1: generate the first random number & prefetch the degree.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            slot.r_ = sfmt_genrand_uint32(sfmt);
            _mm_prefetch((void*)(graph->offset_pair_ + slot.w_.current_), PREFETCH_HINT);
        }
    }

    // Stage 2: generate the position & prefetch the alias slot.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& offset = graph->offset_pair_[slot.w_.current_];
            slot.r_ = slot.r_ % (offset.second - offset.first);
            slot.dr_ = sfmt_genrand_real2(sfmt);
            _mm_prefetch((void*)(slot.alias_ + slot.r_), PREFETCH_HINT);
        }
    }

    // Stage 3: update the walker.
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            auto& alias_slot = slot.alias_[slot.r_];
            slot.w_.current_ = slot.dr_ <= alias_slot.alias_value_ ?
                               alias_slot.first_ : alias_slot.second_;

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }

            slot.w_.length_ += 1;
        }
    }
}

void static_alias_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto degree = graph->offset_[slot.w_.current_ + 1] - graph->offset_[slot.w_.current_];
            auto alias_table = graph->edge_weight_alias_table_ + graph->offset_[slot.w_.current_];

            // Select a position.
            auto r0 = sfmt_genrand_uint32(sfmt);
            auto p0 = r0 % degree;

            // Select the vertex.
            auto r1 = sfmt_genrand_real2(sfmt);
            auto vertex = r1 <= alias_table[p0].alias_value_ ? alias_table[p0].first_ : alias_table[p0].second_;

            // Find the neighbor and update.
            slot.w_.current_ = vertex;

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}

void dynamic_alias_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, int length) {
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        if (!slot.empty_) {
            slot.prev_ = slot.w_.current_;
            auto degree = graph->offset_[slot.w_.current_ + 1] - graph->offset_[slot.w_.current_];
            auto alias_table = slot.alias_;

            // Select a position.
            auto r0 = sfmt_genrand_uint32(sfmt);
            auto p0 = r0 % degree;

            // Select the vertex.
            auto r1 = sfmt_genrand_real2(sfmt);
            auto vertex = r1 <= alias_table[p0].alias_value_ ? alias_table[p0].first_ : alias_table[p0].second_;

            // Find the neighbor and update.
            slot.w_.current_ = vertex;

            if (slot.w_.length_ < length) {
                slot.seq_[slot.w_.length_] = slot.w_.current_;
            }
            slot.w_.length_ += 1;
        }
    }
}


/**
 * Implement the ALIAS sampling for test purpose
 */
void static_alias_AMAC_move(Graph *graph, BufferSlot *ring, sfmt_t *sfmt, AMAC_alias_frame* frames,
        int length) {
    int search_ring_id = 0;
    for (int i = 0; i < RING_SIZE; ++i) {
        BufferSlot& slot = ring[i];
        auto* fr = &frames[search_ring_id];
        if (!slot.empty_) {
            if (fr->state == AMAC_alias_frame::Empty) {
                // Execute stage S0
                fr->init(graph->offset_pair_, graph->edge_weight_alias_table_,
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
                    if (fr->state == AMAC_alias_frame::S0) {
                        fr->execute_S0();
                    }
                    else if (fr->state == AMAC_alias_frame::S1) {
                        fr->execute_S1(sfmt);
                    }
                    else if (fr->state == AMAC_alias_frame::S2) {
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
                        fr->init(graph->offset_pair_, graph->edge_weight_alias_table_,
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
            AMAC_alias_frame& frame = frames[i];
            if (frame.state == AMAC_alias_frame::S0) {
                more_work = true;
                frame.execute_S0();
            }
            else if (frame.state == AMAC_alias_frame::S1) {
                more_work = true;
                frame.execute_S1(sfmt);
            }
            else if (frame.state == AMAC_alias_frame::S2) {
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


#endif //XTRAGRAPHCOMPUTING_ALIAS_SAMPLING_H
