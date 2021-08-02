//
// Created by Shixuan Sun on 11/3/20.
//

#ifndef XTRAGRAPHCOMPUTING_WALKER_H
#define XTRAGRAPHCOMPUTING_WALKER_H

#include "config/type_config.h"
struct WalkerMeta {
    int id_;
    int source_;
    int current_;
    int length_;
    intT* seq_;
};

#endif //XTRAGRAPHCOMPUTING_WALKER_H
