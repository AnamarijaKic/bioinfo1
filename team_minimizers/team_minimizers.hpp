#ifndef TEAM_MINIMIZERS_HPP
#define TEAM_MINIMIZERS_HPP

#include <vector>
#include <tuple>
#include <string>
#include <deque>
#include <unordered_map>
#include <set>

namespace team {

class KMER {
    private:
        bool is_fwd = true;

    public:
        KMER(bool is_fwd_);

        std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
            const char* sequence, unsigned int sequence_len,
            unsigned int kmer_len,
            unsigned int window_len);

        std::string MappKmerBitToString(unsigned int kmer, unsigned int kmer_len);
        unsigned int MappSeqCharPointerToBit(const char* seq, unsigned int kmer_len);
        std::string ReverseComplement(const std::string& kmer);
        std::tuple<unsigned int, unsigned int, bool> GetTupleWithMinFirst(
            const std::deque<std::tuple<unsigned int, unsigned int, bool>>& window);
        std::unordered_map<unsigned int, int> GetMinimizerFrequencies();
        std::set<std::tuple<unsigned int, unsigned int, bool>> GetUniqueMinimizers();
    };
}

#endif