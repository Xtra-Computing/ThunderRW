//
// Created by Shixuan Sun on 11/22/20.
//

#ifndef XTRAGRAPHCOMPUTING_ITS_FRAME_H
#define XTRAGRAPHCOMPUTING_ITS_FRAME_H

#include "types.h"

#define BIT 1

#define choose(i, a, b) ({			\
	const double *ret;			\
	asm("testl %1,%2 ; cmovne %3,%0"	\
		:"=r" (ret)			\
		:"i" (BIT),			\
		 "g" (i),			\
		 "rm" (a),			\
		 "0" (b));			\
	ret; })

template<typename T> struct ITSFrame {
    enum State {RUNNING, EMPTY};
    const T* original;
    const T* base;
    T value;
    int length;
    int half;
    int result;
    int id;
    State state = EMPTY;

    inline void init(const T *src, int64_t begin, int64_t end, T target, int slot_id) {
        original = base = src + begin;
        value = target;
        length = end - begin;
        half = length >> 1;
        _mm_prefetch((void*)(base + half), PREFETCH_HINT);
        id = slot_id;
        state = RUNNING;
    }

    inline bool run() {
        int c = (base[half] <= value);
        base =  choose(c, &base[half], base);
        length -= half;

        if (length > 1) {
            half = length >> 1;
            _mm_prefetch((void*)(base + half), PREFETCH_HINT);
            return false;
        }
        else {
            result = (*base <= value) + base - original;
            return true;
        }
    }
};


#endif //XTRAGRAPHCOMPUTING_ITS_FRAME_H
