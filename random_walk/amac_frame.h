//
// Created by Shixuan Sun on 11/22/20.
//

#ifndef XTRAGRAPHCOMPUTING_AMAC_FRAME_H
#define XTRAGRAPHCOMPUTING_AMAC_FRAME_H

#define AMAC_BIT 1

#define amac_choose(i, a, b) ({			\
	const double *ret;			\
	asm("testl %1,%2 ; cmovne %3,%0"	\
		:"=r" (ret)			\
		:"i" (AMAC_BIT),			\
		 "g" (i),			\
		 "rm" (a),			\
		 "0" (b));			\
	ret; })


struct AMAC_uniform_frame {
    enum State {S0, S1, S2, Empty};
    std::pair<int64_t, int64_t>* offset_pair_;
    intT* adj_;

    intT current_;

    int64_t x_;
    intT value_;
    int id;
    State state = Empty;

    inline void init(std::pair<int64_t, int64_t>* offset_pair, intT* adj,
            intT current, int slot_id) {
        offset_pair_ = offset_pair;
        adj_ = adj;
        current_ = current;
        id = slot_id;
        state = S0;
    }

    inline void execute_S0() {
        _mm_prefetch((void*)(offset_pair_ + current_), PREFETCH_HINT);
        state = S1;
    }

    inline void execute_S1(sfmt_t* sfmt) {
        auto& offset = offset_pair_[current_];
        adj_ = adj_ + offset.first;
        x_ = sfmt_genrand_uint32(sfmt);
        x_ = x_ % (offset.second - offset.first);
        _mm_prefetch((void*)(adj_ + x_), PREFETCH_HINT);
        state = S2;
    }

    inline void execute_S2() {
        value_ = adj_[x_];
        state = Empty;
    }
};

struct AMAC_alias_frame {
    enum State {S0, S1, S2, Empty};
    std::pair<int64_t, int64_t>* offset_pair_;
    AliasSlot* edge_weight_alias_table_;
    intT current_;
    intT value_;
    int64_t x;
    double y;
    int id;
    State state = Empty;

    inline void init(std::pair<int64_t, int64_t>* offset_pair,
                     AliasSlot* edge_weight_alias_table,
                     intT current, int slot_id) {
        offset_pair_ = offset_pair;
        edge_weight_alias_table_ = edge_weight_alias_table;
        current_ = current;
        id = slot_id;
        state = S0;
    }

    inline void execute_S0() {
        _mm_prefetch((void*)(offset_pair_ + current_), PREFETCH_HINT);
        state = S1;
    }

    inline void execute_S1(sfmt_t* sfmt) {
        auto& offset = offset_pair_[current_];
        x = sfmt_genrand_uint32(sfmt);
        x = offset.first + (x % (offset.second - offset.first));
        y = sfmt_genrand_real2(sfmt);
        _mm_prefetch((void*)(edge_weight_alias_table_ + x), PREFETCH_HINT);
        state = S2;
    }

    inline void execute_S2() {
        auto& alias_slot = edge_weight_alias_table_[x];
        value_ = y <= alias_slot.alias_value_ ?
                           alias_slot.first_ : alias_slot.second_;
        state = Empty;
    }
};

struct AMAC_max_rj_frame {
    enum State {S0, S1, S2, S3, S4, Empty};
    intT* adj_;
    std::pair<int64_t, int64_t>* offset_pair_;
    intT current_;
    double max_weight_;

    intT value_;

    std::pair<int64_t, int64_t> offset_;
    int64_t result_;
    int64_t length_;
    double cur_weight_;
    int id;
    State state = Empty;

    inline void init(intT* adj, std::pair<int64_t, int64_t>* offset_pair, intT current,
            double max_weight, int slot_id) {
        adj_ = adj;
        offset_pair_ = offset_pair;
        current_ = current;
        max_weight_ = max_weight;
        id = slot_id;
        state = S0;
    }

    inline void execute_S0() {
        _mm_prefetch((void*)(offset_pair_ + current_), PREFETCH_HINT);
        state = S1;
    }

    inline void execute_S1() {
        offset_ = offset_pair_[current_];
        adj_ = adj_ + offset_.first;
        length_ = offset_.second - offset_.first;
        state = S2;
    }

    inline void execute_S2(sfmt_t* sfmt) {
        result_ = sfmt_genrand_uint32(sfmt) % length_;
        cur_weight_ = sfmt_genrand_real2(sfmt) * max_weight_;
        _mm_prefetch((void*)(adj_ + result_), PREFETCH_HINT);
        state = S3;
    }

    template<typename F> void execute_S3(F& f, WalkerMeta& w) {
        double edge_weight = f.weight(w, current_, adj_[result_], offset_.first + result_);
        state = cur_weight_ <= edge_weight ? S4 : S2;
    }

    inline void execute_S4() {
        value_ = adj_[result_];
        state = Empty;
    }
};

struct AMAC_rj_frame {
    enum State {S0, S1, S2, S3, S4, S5, S6, Empty};
    intT* adj_;
    std::pair<int64_t, int64_t>* offset_pair_;
    double* edge_weight_rejection_max_;
    double* base_;
    intT current_;

    intT value_;

    std::pair<int64_t, int64_t> offset_;
    double max_weight_;
    double cur_weight_;
    int64_t result_;
    int64_t length_;
    int id;
    State state = Empty;

    inline void init(intT* adj, std::pair<int64_t, int64_t>* offset_pair,
                     double* edge_weight_rejection_max, double* base,
                     intT current, int slot_id) {
        adj_ = adj;
        offset_pair_ = offset_pair;
        edge_weight_rejection_max_ = edge_weight_rejection_max;
        base_ = base;
        current_ = current;
        id = slot_id;
        state = S0;
    }

    inline void execute_S0() {
        _mm_prefetch((void*)(offset_pair_ + current_), PREFETCH_HINT);
        state = S1;
    }

    inline void execute_S1() {
        offset_ = offset_pair_[current_];
        base_ = base_ + offset_.first;
        adj_ = adj_ + offset_.first;
        length_ = offset_.second - offset_.first;
        _mm_prefetch((void*)(edge_weight_rejection_max_ + current_), PREFETCH_HINT);
        state = S2;
    }

    inline void execute_S2() {
        max_weight_ = edge_weight_rejection_max_[current_];
        state = S3;
    }

    inline void execute_S3(sfmt_t* sfmt) {
        result_ = sfmt_genrand_uint32(sfmt) % length_;
        cur_weight_ = sfmt_genrand_real2(sfmt) *  max_weight_;
        _mm_prefetch((void*)(base_ + result_), PREFETCH_HINT);
        state = S4;
    }

    inline void execute_S4() {
        state = cur_weight_ <= base_[result_] ? S5 : S3;
    }

    inline void execute_S5() {
        _mm_prefetch((void*)(adj_ + result_), PREFETCH_HINT);
        state = S6;
    }

    inline void execute_S6() {
        value_ = adj_[result_];
        state = Empty;
    }
};

struct AMAC_its_frame {
    enum State {S0, S1, S2, S3, S4, S5, S6, Empty};
    intT* adj_;
    std::pair<int64_t, int64_t>* offset_pair_;
    const double* original_;

    intT current_;

    intT value_;

    std::pair<int64_t, int64_t> offset_;
    double target_;
    int64_t result_;
    int64_t length_;
    int64_t half_;
    const double *base_;
    int id;
    State state = Empty;

    inline void init(intT* adj, std::pair<int64_t, int64_t>* offset_pair,
                     const double* edge_weight_prefix_sum, intT current, int slot_id) {
        adj_ = adj;
        offset_pair_ = offset_pair;
        original_ = edge_weight_prefix_sum;
        current_ = current;
        id = slot_id;
        state = S0;
    }

    inline void execute_S0() {
        _mm_prefetch((void*)(offset_pair_ + current_), PREFETCH_HINT);
        state = S1;
    }

    inline void execute_S1() {
        offset_ = offset_pair_[current_];
        original_ = original_ + offset_.first;
        base_ = original_;
        adj_ = adj_ + offset_.first;
        length_ = offset_.second - offset_.first;
        _mm_prefetch((void*)(original_ + length_ - 1), PREFETCH_HINT);
        state = S2;
    }

    inline void execute_S2(sfmt_t* sfmt) {
        target_ = sfmt_genrand_real2(sfmt);;
        target_ *= original_[length_ - 1];
        state = S3;
    }

    inline void execute_S3() {
        half_ = length_ >> 1;
        _mm_prefetch((void*)(base_ + half_), PREFETCH_HINT);
        state = S4;
    }

    inline void execute_S4() {
        int c = (base_[half_] <= target_);
        base_ = amac_choose(c, &base_[half_], base_);
        length_ -= half_;

        state = length_ > 1 ? S3 : S5;
    }

    inline void execute_S5() {
        result_ = (*base_ <= target_) + base_ - original_;
        _mm_prefetch((void*)(adj_ + result_), PREFETCH_HINT);
        state = S6;
    }

    inline void execute_S6() {
        value_ = adj_[result_];
        state = Empty;
    }
};

#endif //XTRAGRAPHCOMPUTING_AMAC_FRAME_H
