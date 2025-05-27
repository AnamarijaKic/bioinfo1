#include "team_minimizers.hpp"
#include <deque>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <bitset>
#include <tuple>

using namespace std;

namespace team {
    // Mapping from bit shifted values to string
    string MappKmerBitToString(unsigned int kmer, unsigned int kmer_len){
        std::string mapp(kmer_len, 'X');
        // Value mapping for positions (original order)
        unordered_map<unsigned int, char> base_value = {
            {0,'C'},
            {1,'A'},
            {2,'T'},
            {3,'G'}
        };

        for (int i = kmer_len - 1; i >= 0; --i){
            unsigned int bits = kmer & 0b11;
            mapp[i] = base_value[bits];
            kmer>>=2;
        }

        return mapp;
    }
    
    // Mapping function for a k-mer
    unsigned int MappKmerStringToBit(const string& kmer){
        unsigned int mapp = 0;
        // Value mapping for positions (original order)
        unordered_map<char, unsigned int> base_value = {
            {'C', 0},
            {'A', 1},
            {'T', 2},
            {'G', 3}
        };

        for (int i=0; i<kmer.size(); i++){
            mapp<<=2;
            mapp |= base_value[kmer[i]];
        }

        return mapp;
    }

    // Mapping function for a k-mer - String
    string MappKmerStringToString(const string& kmer){
        string mapp = "";
        // Value mapping for positions (original order)
        unordered_map<char, unsigned int> base_value = {
            {'C', '0'},
            {'A', '1'},
            {'T', '2'},
            {'G', '3'}
        };

        for (int i=0; i<kmer.size(); i++){
            unsigned int val;
            val = base_value[kmer[i]];
            mapp += val;
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

    // Find kmer inside of number that holds whole sequence, with 2 bits for each sequence
    unsigned int ExtractKmer(unsigned int packed_seq, int i, int k) {
        unsigned int shift = 2 * i; 
        unsigned int mask = (1ULL << (2 * k)) - 1; // k*2 bits set to 1
        return (packed_seq >> shift) & mask;
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

        // Custom Comparator for min-heap because priority_queue is by default max-heap
        auto cmp = [](const tuple<unsigned int, unsigned int, bool>& a,
                        const tuple<unsigned int, unsigned int, bool>& b) {
            return a > b; 
        };


        // Structure that holds all minimizers
        vector<tuple<unsigned int, unsigned int, bool>> minimizers;
        vector<tuple<unsigned int, unsigned int, bool>> minimizers2;

        // Sequence and reverse sequence in string format (but WITH letters for bases)
        string sequence_str(sequence, sequence_len);
        string rc_sequence_str = ReverseComplement(sequence_str);

        // Sequence and reverse sequence in string fomrat (but WITH numbers for bases)
        string mapp_seq_fwd_str = MappKmerStringToString(sequence_str);
        string mapp_seq_rev_str = MappKmerStringToString(rc_sequence_str);

        // Sequence and reverse sequence in number format (but 2 bits for one base)
        unsigned int mapp_seq_fwd_bit = MappKmerStringToBit(sequence_str);
        unsigned int mapp_seq_rev_bit = MappKmerStringToBit(rc_sequence_str);

        const char* rc_sequence = rc_sequence_str.c_str();
        const char* mapp_seq_fwd_str_ptr = mapp_seq_fwd_str.c_str();
        const char* mapp_seq_rev_str_ptr = mapp_seq_rev_str.c_str();

        if (sequence_len < kmer_len || window_len == 0) return minimizers;

       
        deque<tuple<unsigned int, unsigned int, bool>> window_monotonic;
        deque<tuple<unsigned int, unsigned int, bool>> window;
        // It contains already sorted minimizers
        priority_queue< tuple<unsigned int, unsigned int, bool>,
                        vector<tuple<unsigned int, unsigned int, bool>>,
                        decltype(cmp)> 
        window_heap(cmp);

        // BEGINING -> end-minimizer:
        for(unsigned int u=kmer_len; u<window_len;u++){
            deque<tuple<unsigned int, unsigned int, bool>> end_window;
            // It contains already sorted minimizers
            priority_queue< tuple<unsigned int, unsigned int, bool>,
                            vector<tuple<unsigned int, unsigned int, bool>>,
                            decltype(cmp)> 
            end_window_heap(cmp);

            for (unsigned int i = 0; i <= u-kmer_len; i++){
                // string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
                // string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);
                unsigned int mapp_fwd_bit = ExtractKmer(mapp_seq_fwd_bit, sequence_len-kmer_len-i, kmer_len);
                unsigned int mapp_rev_bit = ExtractKmer(mapp_seq_rev_bit, i, kmer_len);

                unsigned int min_kmer = min(mapp_fwd_bit, mapp_rev_bit);
                bool is_fwd = (min_kmer == mapp_fwd_bit);

                end_window.push_back(make_tuple(min_kmer, i+1, is_fwd));
                end_window_heap.push(make_tuple(min_kmer, i+1, is_fwd));
            }
            minimizers.push_back(GetTupleWithMinFirst(end_window));
        }

        // NORMAL WINDOWS
        for (unsigned int i = 0; i <= sequence_len-kmer_len; i++){
            if(i>=window_len){ window.pop_front(); }
            // string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
            // string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);

            // unsigned int min_kmer = min(stoull(mapp_fwd_str), stoull(mapp_rev_str));
            // bool is_fwd = (min_kmer == stoull(mapp_fwd_str));
            // string min_kmer_str = min(mapp_fwd_str, mapp_rev_str);
            // bool is_fwd = (min_kmer_str == mapp_fwd_str);
            unsigned int mapp_fwd_bit = ExtractKmer(mapp_seq_fwd_bit, sequence_len-kmer_len-i, kmer_len);
            unsigned int mapp_rev_bit = ExtractKmer(mapp_seq_rev_bit, i, kmer_len);

            unsigned int min_kmer = min(mapp_fwd_bit, mapp_rev_bit);
            bool is_fwd = (min_kmer == mapp_fwd_bit);

            // unsigned int min_kmer = stoull(min_kmer_str);
            // window_monotonic.emplace_back(min_kmer, i, is_fwd);
            window.push_back(make_tuple(min_kmer, i+1, is_fwd));
            
            // window_heap.push(make_tuple(min_kmer, i+1, is_fwd));
            if(i>=(window_len-1)){
                minimizers.push_back(GetTupleWithMinFirst(window)); 
            }
        }

        // KRAJ -> end-minimizer:
        for(unsigned int u=kmer_len; u<window_len;u++){
            if (sequence_len < u) break;
            deque<tuple<unsigned int, unsigned int, bool>> end_window;
            unsigned int start = sequence_len - u;

            for (unsigned int i = start; i <= sequence_len - kmer_len; ++i){
                unsigned int mapp_fwd_bit = ExtractKmer(mapp_seq_fwd_bit, sequence_len-kmer_len-i, kmer_len);
                unsigned int mapp_rev_bit = ExtractKmer(mapp_seq_rev_bit, i, kmer_len);

                unsigned int min_kmer = min(mapp_fwd_bit, mapp_rev_bit);
                bool is_fwd = (min_kmer == mapp_fwd_bit);

                end_window.push_back(make_tuple(min_kmer, i+1, is_fwd));
            }
            minimizers.push_back(GetTupleWithMinFirst(end_window));
        }
        return minimizers;
    }
}

// KORISNI PRINTOVI I KODOVI

// cout<<"Seq: "<<sequence_str<<endl;
// cout<<"Rc_seq: "<<rc_sequence_str<<endl;
// cout<<"STRING SEKVENCE: "<<endl;
// cout<<mapp_seq_fwd_str<<endl;
// cout<<mapp_seq_rev_str<<endl;
// cout<<"BIT SEKVENCE: "<<endl;
// cout<<"prva: "<<endl;
// cout<<mapp_seq_fwd_bit<<endl;
// cout<<MappKmerBitToString(mapp_seq_fwd_bit, sequence_len)<<endl;
// cout<<"druga: "<<endl;
// cout<<mapp_seq_rev_bit<<endl;
// cout<<MappKmerBitToString(mapp_seq_rev_bit, sequence_len)<<endl;

// cout<<kmer<<std::endl;
// cout<<mapp_fwd<<std::endl;
// cout<<rc_kmer<<std::endl;
// cout<<mapp_rev<<std::endl;
// cout<<endl;

// unsigned int mapp_fwd = MappKmer(kmer);
// unsigned int mapp_rev = MappKmer(rc_seq_kmer);
// cout<<"INTEGERI: "<<endl;
// cout<<mapp_fwd<<endl;
// cout<<mapp_rev<<endl;

// cout<<"STRINGOVI: "<<endl;
// string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
// string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);
// cout<<mapp_fwd_str<<endl;
// cout<<mapp_rev_str<<endl;

// string kmer(sequence + i, kmer_len); // std::string(const char* s, size_t n);
// string rc_kmer = ReverseComplement(kmer);

// cout<<"STRINGOVI: "<<endl;
// string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
// string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);
// cout<<mapp_fwd_str<<endl;
// cout<<mapp_rev_str<<endl;

// auto& za_ispis = GetTupleWithMinFirst(end_window);
// cout<<"MINIMIZERI: "<<endl;
// unsigned int first_value = get<0>(za_ispis); 
// cout<<MappKmerBitToString(first_value,kmer_len)<<endl;
// cout<<MappKmerStringToString(MappKmerBitToString(first_value,kmer_len))<<endl;