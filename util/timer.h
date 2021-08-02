//
// Created by Shixuan Sun on 2020/6/29.
//

#ifndef XTRAGRAPHCOMPUTING_TIMER_H
#define XTRAGRAPHCOMPUTING_TIMER_H

#include <iostream>
#include <chrono>

class Timer {
public:
    Timer() : beg_(clock_::now()) {}

    void reset() { beg_ = clock_::now(); }

    double elapsed() const {
        return std::chrono::duration_cast<second_>
                (clock_::now() - beg_).count();
    }

    double elapsed_and_reset() {
        double elapsed = std::chrono::duration_cast<second_>
                (clock_::now() - beg_).count();
        beg_ = clock_::now();
        return elapsed;
    }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

#endif //XTRAGRAPHCOMPUTING_TIMER_H
