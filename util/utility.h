//
// Created by Shixuan Sun on 2020/10/19.
//

#ifndef UTILITY_H
#define UTILITY_H

#include <xmmintrin.h>
#include <algorithm>
#include <sstream>
#include <random>
#include <iostream>

#include "config/type_config.h"

class Utility {
public:
    template<typename T> static int upper_bound_branchless(T *src, int begin, int end, T target) {
        const T *base = src + begin;
        int n = end - begin;
        while (n > 1) {
            int half = n >> 1;
            base = (base[half] <= target) ? &base[half] : base;
            n -= half;
        }
        return (*base <= target) + base - (src + begin);
    }

    template<typename T> static int test(T *src, int begin, int end, T target) {
        // Initialize.
        const T *base = src + begin;
        int n = end - begin;
        if (n <= 1) {
            return (*base <= target) + base - (src + begin);
        }

        int half = n >> 1;
        // prefetch base[half].

        while (true) {
            // Run.
            base = (base[half] <= target) ? &base[half] : base;
            n -= half;

            if (n <= 1) {
                return (*base <= target) + base - (src + begin);
            }

            half = n >> 1;
            // prefetch base[half].
        }

    }


    template<typename T> static int binary_search(T *src, int begin, int end, T target) {
        int offset_begin = begin;
        int offset_end = end;
        while (offset_end - offset_begin >= 16) {
            auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_begin) + offset_end) / 2);
            _mm_prefetch((char *) &src[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
            _mm_prefetch((char *) &src[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
            if (src[mid] == target) {
                return mid;
            } else if (src[mid] < target) {
                offset_begin = mid + 1;
            } else {
                offset_end = mid;
            }
        }

        // linear search fallback
        for (auto offset = offset_begin; offset < offset_end; ++offset) {
            if (src[offset] >= target) {
                return offset;
            }
        }

        return offset_end;
    }
    //template <typename T> static uint32_t branchfree_search(T*a, uint32_t length, T x);
    static uint32_t branchfree_search(int64_t*a, uint32_t n, int64_t x) {
        using I = uint32_t;
        //const T *base = a;
        const int64_t*base = a;
        while (n > 1) {
            I half = n / 2;
            __builtin_prefetch(base + half / 2, 0, 0);
            __builtin_prefetch(base + half + half / 2, 0, 0);
            base = (base[half] < x) ? base + half : base;
            n -= half;
        }
        return (*base < x) + base - a;
    }

    template<typename T> static void shuffle(T* src, uint32_t n, bool debug_mode = true) {
        std::random_device rd;
        std::mt19937 rng(debug_mode ? 0 : rd());
        std::shuffle(src, src + n, rng);
    }

    template<typename T> static void sequential_prefix_sum(T* input_array, T* output_array, uint32_t length) {
        output_array[0] = input_array[0];
        for(int i = 1; i < length; i++){
            output_array[i] = output_array[i - 1] + input_array[i];
        }
    }

    static void split(const std::string &s, char delimiter, std::vector<std::string>& tokens) {
        std::istringstream iss(s);
        std::string item;
        while (std::getline(iss, item, delimiter)) {
            tokens.emplace_back(item);
        }
    }

};

#endif //UTILITY_H
