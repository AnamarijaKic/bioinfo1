#ifndef TEAM_ALIGNMENT_HPP
#define TEAM_ALIGNMENT_HPP

#include <string>

// Anamarija Kic
namespace team {
    enum class AlignmentType {
        global, // == NeedlemanWunsch,
        local, // == SmithWaterman,
        semiGlobal // == Gotoh
    };

    int Align(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        AlignmentType type,
        int match,
        int mismatch,
        int gap,
        std::string* cigar = nullptr,
        unsigned int* target_begin = nullptr
    );

    //std::string needleman_wunsch(const std::string& a, const std::string& b);
    //std::string smith_waterman(const std::string& a, const std::string& b);
    //std::string gotoh(const std::string& a, const std::string& b);
}

#endif
