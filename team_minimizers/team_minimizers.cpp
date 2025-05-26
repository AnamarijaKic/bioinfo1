#include "team_minimizers.hpp"
#include <deque>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <tuple>

using namespace std;

namespace team {
    
    // Mapping function for a k-mer
    unsigned int MappKmer(const string& kmer){
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

    // Mapping function for a k-mer - String
    string MappKmerString(const string& kmer){
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

        string sequence_str(sequence, sequence_len);
        string rc_sequence_str = ReverseComplement(sequence_str);

        string mapp_seq_fwd_str = MappKmerString(sequence_str);
        string mapp_seq_rev_str = MappKmerString(rc_sequence_str);

        const char* rc_sequence = rc_sequence_str.c_str();
        const char* mapp_seq_fwd_str_ptr = mapp_seq_fwd_str.c_str();
        const char* mapp_seq_rev_str_ptr = mapp_seq_rev_str.c_str();

        if (sequence_len < kmer_len || window_len == 0) return minimizers;

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
                string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
                string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);
                
                unsigned int min_kmer = min(stoull(mapp_fwd_str), stoull(mapp_rev_str));
                bool is_fwd = (min_kmer == stoull(mapp_fwd_str));

                end_window.push_back(make_tuple(min_kmer, i+1, is_fwd));
                end_window_heap.push(make_tuple(min_kmer, i+1, is_fwd));
            }
            // cout<<"MINIMIZERI: "<<endl;
            // const auto& top = end_window_heap.top();
            // unsigned int first_value = get<0>(top);  // First element of the tuple
            // cout << "First value from top tuple: " << first_value << endl;

            // const auto& top2 = GetTupleWithMinFirst(end_window);
            // unsigned int first_value2 = get<0>(top2);  // First element of the tuple
            // cout << "Value from my function: "<< first_value2<<endl;

            minimizers.push_back(GetTupleWithMinFirst(end_window));
        }

        // UNUTARNJI WINDOWI 
        for (unsigned int i = 0; i <= sequence_len-kmer_len; i++){
            if(i>=window_len){ window.pop_front(); }
            string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
            string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);

            unsigned int min_kmer = min(stoull(mapp_fwd_str), stoull(mapp_rev_str));
            bool is_fwd = (min_kmer == stoull(mapp_fwd_str));

            window.push_back(make_tuple(min_kmer, i+1, is_fwd));
            window_heap.push(make_tuple(min_kmer, i+1, is_fwd));

            if(i>=(window_len-1)){ minimizers.push_back(GetTupleWithMinFirst(window)); }
        }

        // KRAJ -> end-minimizer:
        for(unsigned int u=kmer_len; u<window_len;u++){
            if (sequence_len < u) break;
            deque<tuple<unsigned int, unsigned int, bool>> end_window;
            unsigned int start = sequence_len - u;

            for (unsigned int i = start; i <= sequence_len - kmer_len; ++i){
                string mapp_fwd_str(mapp_seq_fwd_str_ptr+i, kmer_len);
                string mapp_rev_str(mapp_seq_rev_str_ptr+sequence_len-kmer_len-i, kmer_len);

                unsigned int min_kmer = min(stoull(mapp_fwd_str), stoull(mapp_rev_str));
                bool is_fwd = (min_kmer == stoull(mapp_fwd_str));

                end_window.push_back(make_tuple(min_kmer, i+1, is_fwd));
            }
            minimizers.push_back(GetTupleWithMinFirst(end_window));
        }

        return minimizers;
    }
}

// KORISNI PRINTOVI I KODOVI

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

// cout<<"Seq: "<<sequence_str<<endl;
// cout<<"Rc_seq: "<<rc_sequence_str<<endl;
// cout<<"STRING SEKVENCE: "<<endl;
// cout<<mapp_seq_fwd_str<<endl;
// cout<<mapp_seq_rev_str<<endl;