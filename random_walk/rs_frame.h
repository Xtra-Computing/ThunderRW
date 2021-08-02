//
// Created by Shixuan Sun on 11/28/20.
//

#ifndef XTRAGRAPHCOMPUTING_RS_FRAME_H
#define XTRAGRAPHCOMPUTING_RS_FRAME_H

#include "SFMT.h"

template<typename T> struct RSFrame {
    enum State {RUNNING, EMPTY};
    const T* base;
    T max_value;
    T cur_value;
    int length;
    int id;
    int result;
    State state = EMPTY;

    inline void init(const T *src, int64_t begin, int64_t end, T max, int slot_id, sfmt_t* sfmt) {
        base = src + begin;
        max_value = max;
        length = end - begin;
        id = slot_id;

        result = sfmt_genrand_uint32(sfmt) % length;
        _mm_prefetch((void*)(base + result), _MM_HINT_T0);
        cur_value = sfmt_genrand_real2(sfmt) * max_value;

        state = RUNNING;
    }

    inline bool run(sfmt_t* sfmt) {
        if (cur_value <= base[result]) {
            return true;
        }
        else {
            result = sfmt_genrand_uint32(sfmt) % length;
            cur_value = sfmt_genrand_real2(sfmt) * max_value;

            _mm_prefetch((void*)(base + result), PREFETCH_HINT);
            return false;
        }
    }
};

template<typename T> struct MaxRSFrame {
    enum State {RUNNING, EMPTY};
    const T* base;
    double max_value;
    double cur_value;
    int length;
    int id;
    int result;
    State state = EMPTY;

    inline void init(const T *src, int64_t begin, int64_t end, double max, int slot_id, sfmt_t* sfmt) {
        base = src + begin;
        max_value = max;
        length = end - begin;
        id = slot_id;

        result = sfmt_genrand_uint32(sfmt) % length;
        _mm_prefetch((void*)(base + result), PREFETCH_HINT);
        cur_value = sfmt_genrand_real2(sfmt) * max_value;

        state = RUNNING;
    }

    template<typename F> bool run(sfmt_t* sfmt, F& f, WalkerMeta& w, int64_t offset) {
        if (cur_value <= f.weight(w, w.current_, base[result], offset + result)) {
           result = base[result]; 
	   return true;
        }
        else {
            result = sfmt_genrand_uint32(sfmt) % length;
            cur_value = sfmt_genrand_real2(sfmt) * max_value;

            _mm_prefetch((void*)(base + result), PREFETCH_HINT);
            return false;
        }
    }
};
#endif //XTRAGRAPHCOMPUTING_RS_FRAME_H
