#include "team_minimizers.hpp"
#include <deque>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <bitset>
#include <tuple>
#include <set>

using namespace std;

// Anamarija Kic
namespace team {

    // class KMER{
        // private:
            bool is_fwd = true;
            bool frequencies = true;
            unordered_map<unsigned int, int> minimizer_frequencies;
            set<tuple<unsigned int, unsigned int, bool>> unique_minmizers;

        // public:
            // Constructor for class
            KMER::KMER(bool is_fwd_) : is_fwd(is_fwd_){}

            // Return unique minimizers
            set<tuple<unsigned int, unsigned int, bool>> KMER::GetUniqueMinimizers(){
                return unique_minmizers;
            }

            // Return frequencies of minimizers
            unordered_map<unsigned int, int> KMER::GetMinimizerFrequencies(){
                return minimizer_frequencies;
            }

            // Enable calculating frequencies of minimizers
            void KMER::SetFrequenciesCount(bool set){
                frequencies = set;
            }

            // Mapping from bit shifted values to string
            string KMER::MappKmerBitToString(unsigned int kmer, unsigned int kmer_len){
                std::string mapp(kmer_len, 'X');
                // Value mapping for positions (original order)
                // unordered_map<unsigned int, char> base_value = {
                //     {0,'C'},
                //     {1,'A'},
                //     {2,'T'},
                //     {3,'G'}
                // };
                unordered_map<unsigned int, char> base_value = {
                    {0,'0'},
                    {1,'1'},
                    {2,'2'},
                    {3,'3'}
                };

                for (int i = kmer_len - 1; i >= 0; --i){
                    unsigned int bits = kmer & 0b11;
                    mapp[i] = base_value[bits];
                    kmer>>=2;
                }

                return mapp;
            }

            // Mapping function for a k-mer
            unsigned int KMER::MappSeqCharPointerToBit(const char* seq, unsigned int kmer_len){
                unsigned int mapp = 0;
                // Value mapping for positions (original order)
                unordered_map<char, unsigned int> base_value = {
                    {'C', 0},
                    {'A', 1},
                    {'T', 2},
                    {'G', 3}
                };

                for (int i=0; i<kmer_len; i++){
                    mapp<<=2;
                    mapp |= base_value[seq[i]];
                }

                return mapp;
            }

            // Reverse complement of a DNA string
            string KMER::ReverseComplement(const string& kmer) { 
                // if we have base ACG reverse complement is made in this way => ACG -> GCA -> CGT!!
                // 1. reversing the sequence
                string rc(kmer.rbegin(), kmer.rend());
                // 2. complement sequence
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
            tuple<unsigned int, unsigned int, bool> KMER::GetTupleWithMinFirst(
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

            vector<tuple<unsigned int, unsigned int, bool>> KMER::Minimize(
                        const char* sequence, unsigned int sequence_len,
                        unsigned int kmer_len,
                        unsigned int window_len) {

                // Clear this two private variables
                minimizer_frequencies.clear();
                unique_minmizers.clear();

                
                // Structure that holds all minimizers
                vector<tuple<unsigned int, unsigned int, bool>> minimizers;

                // Sequence and reverse sequence in string format (but WITH letters for bases)
                // string sequence_str(sequence, sequence_len);
                // string rc_sequence_str = ReverseComplement(sequence_str);
                // const char* rc_sequence = rc_sequence_str.c_str();

                if (sequence_len < kmer_len || window_len == 0) return minimizers;
            
                deque<tuple<unsigned int, unsigned int, bool>> window;
                // deque<tuple<unsigned int, unsigned int, bool>> window_rev;    

                // BEGINING -> end-minimizer:
                for(unsigned int u=kmer_len; u<(window_len+kmer_len-1);u++){
                    deque<tuple<unsigned int, unsigned int, bool>> end_window;
                    // deque<tuple<unsigned int, unsigned int, bool>> end_window_rev;

                    for (unsigned int i = 0; i <= u-kmer_len; i++){
                            
                            unsigned int mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);

                            // unsigned int mapp_bit;
                            // if(is_fwd){ mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                            // } else{     mapp_bit = MappSeqCharPointerToBit(sequence+sequence_len-kmer_len-i, kmer_len);
                            // }
                            end_window.push_back(make_tuple(mapp_bit, i+1, is_fwd));
                        
                        }

                        auto& best = GetTupleWithMinFirst(end_window);
                        minimizers.push_back(best);
                        unique_minmizers.insert(best);

                        if(frequencies){
                            unsigned int first_value = get<0>(best); 
                            minimizer_frequencies[first_value]++;
                        }
                }

                // NORMAL WINDOWS
                for (unsigned int i = 0; i <= sequence_len-kmer_len; i++){
                    if(i>=window_len){ window.pop_front();}

                    unsigned int mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                    // unsigned int mapp_bit;
                    // if(is_fwd){ mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                    // } else{     mapp_bit = MappSeqCharPointerToBit(sequence+sequence_len-kmer_len-i, kmer_len);
                    // }
                    
                    window.push_back(make_tuple(mapp_bit, i+1, is_fwd));

                    if(i>=(window_len-1)){
                        auto& best = GetTupleWithMinFirst(window);
                        minimizers.push_back(best);
                        unique_minmizers.insert(best);

                        if(frequencies){
                            unsigned int first_value = get<0>(best); 
                            minimizer_frequencies[first_value]++;
                        }
                    }
                }

                // END -> end-minimizer:
                for(unsigned int u=kmer_len; u<(window_len+kmer_len-1);u++){
                    if (sequence_len < u) break;
                    deque<tuple<unsigned int, unsigned int, bool>> end_window;
                    // deque<tuple<unsigned int, unsigned int, bool>> end_window_rev;
                    unsigned int start = sequence_len - u;

                    for (unsigned int i = start; i <= sequence_len - kmer_len; ++i){

                        unsigned int mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);

                        // unsigned int mapp_bit;
                        // if(is_fwd){ mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                        // } else{     mapp_bit = MappSeqCharPointerToBit(sequence+sequence_len-kmer_len-i, kmer_len);
                        // }

                        end_window.push_back(make_tuple(mapp_bit, i+1, is_fwd)); 
                    }
                    auto& best = GetTupleWithMinFirst(end_window);
                    minimizers.push_back(best);
                    unique_minmizers.insert(best);

                    if(frequencies){
                        unsigned int first_value = get<0>(best); 
                        minimizer_frequencies[first_value]++;
                    }
                }
                return minimizers;

            }

    // };
}