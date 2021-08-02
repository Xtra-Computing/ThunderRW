//
// Created by Shixuan Sun on 11/5/20.
//

#ifndef XTRAGRAPHCOMPUTING_COMMAND_PARSER_H
#define XTRAGRAPHCOMPUTING_COMMAND_PARSER_H

#include <string>
#include <algorithm>
#include <numeric>

class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i = 1; i < argc; ++i)
            tokens_.emplace_back(argv[i]);
    }

    std::string get_cmd_option(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(tokens_.begin(), tokens_.end(), option);
        if (itr != tokens_.end() && ++itr != tokens_.end()){
            return *itr;
        }
        return "";
    }

    bool check_cmd_option_exists(const std::string &option) const{
        return std::find(tokens_.begin(), tokens_.end(), option)
               != tokens_.end();
    }

    std::string get_cmd() {
        return std::accumulate(tokens_.begin(), tokens_.end(), std::string(" "));
    }

private:
    std::vector<std::string> tokens_;
};

#endif //XTRAGRAPHCOMPUTING_COMMAND_PARSER_H
