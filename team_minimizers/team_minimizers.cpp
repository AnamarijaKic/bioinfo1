#include "team_minimizers.hpp"
#include <deque>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <iostream>

using namespace std;

namespace team {
    
    // Mapping function for a k-mer
    unsigned int MappKmer(const string& kmer, int index){
        unsigned int mapp = 0;
        // Value mapping for positions (original order)
        unordered_map<char, unsigned int> base_value = {
            {'C', 0},
            {'A', 1},
            {'T', 2},
            {'G', 3}
        };

        for (int i=0; i<kmer.size(); i++){
            unsigned int val;
            val = base_value[kmer[i]];
            mapp = mapp * 10 + val;
        }

        return mapp;
    }

    // Reverse complement of a DNA string
    string ReverseComplement(const string& kmer) { 
        // ako imamo bazu ACG reverse complement nastane kao ACG -> GCA -> CGT!!
        // 1. reversamo sekvencu
        string rc(kmer.rbegin(), kmer.rend()); //https://cplusplus.com/reference/string/string/rbegin/
        // 2. komplementiramo sekvencu
        for (char& c : rc) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'G': c = 'C'; break;
                case 'C': c = 'G'; break;
            }
        }
        return rc;
    }

    // Finding minimizer of window
    tuple<unsigned int, unsigned int, bool> GetTupleWithMinFirst(
        const deque<std::tuple<unsigned int, unsigned int, bool>>& window) {
        
        tuple<unsigned int, unsigned int, bool> min_tuple;
        unsigned int min_value = numeric_limits<unsigned int>::max();

        for (const auto& t : window) {
            if (get<0>(t) < min_value) {
                min_value = get<0>(t);
                min_tuple = t;
            }
        }

        return min_tuple;
    }

    vector<tuple<unsigned int, unsigned int, bool>> Minimize(
        const char* sequence, unsigned int sequence_len,
        unsigned int kmer_len,
        unsigned int window_len) {

        vector<tuple<unsigned int, unsigned int, bool>> minimizers;
        
        
        if (sequence_len < kmer_len || window_len == 0) return minimizers;

        deque<tuple<unsigned int, unsigned int, bool>> window;
        

        // PRVI WINDOW - moramo cijeli prozor napuniti
        for (unsigned int i = 0; i <= window_len-1; i++){
            string kmer(sequence + i, kmer_len); // std::string(const char* s, size_t n);
            string rc_kmer = ReverseComplement(kmer);

            unsigned int mapp_fwd = MappKmer(kmer, i%2);
            unsigned int mapp_rev = MappKmer(rc_kmer, i%2);

            unsigned int min_kmer = min(mapp_fwd, mapp_rev);
            bool is_fwd = (min_kmer == mapp_fwd);

            window.push_back(make_tuple(min_kmer, i+1, is_fwd));

            /*cout<<kmer<<std::endl;
            cout<<mapp_fwd<<std::endl;
            cout<<rc_kmer<<std::endl;
            cout<<mapp_rev<<std::endl;
            cout<<endl;*/
        }

        minimizers.push_back(GetTupleWithMinFirst(window));

        // OSTALI WINDOWI - moramo prvi element popati i zadnji pushati
        for (unsigned int i = window_len; i <= sequence_len-kmer_len; i++){
            window.pop_front();
            string kmer(sequence + i, kmer_len); // std::string(const char* s, size_t n);
            string rc_kmer = ReverseComplement(kmer);

            unsigned int mapp_fwd = MappKmer(kmer, i%2);
            unsigned int mapp_rev = MappKmer(rc_kmer, i%2);

            unsigned int min_kmer = min(mapp_fwd, mapp_rev);
            bool is_fwd = (min_kmer == mapp_fwd);

            window.push_back(make_tuple(min_kmer, i+1, is_fwd));

            minimizers.push_back(GetTupleWithMinFirst(window));

            /*cout<<kmer<<std::endl;
            cout<<mapp_fwd<<std::endl;
            cout<<rc_kmer<<std::endl;
            cout<<mapp_rev<<std::endl;
            cout<<endl;*/
        }

        return minimizers;
    }
}
