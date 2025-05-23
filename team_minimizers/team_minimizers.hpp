#ifndef TEAM_MINIMIZERS_HPP
#define TEAM_MINIMIZERS_HPP

#include <vector>
#include <tuple>
#include <string>

namespace team {

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);
}

#endif