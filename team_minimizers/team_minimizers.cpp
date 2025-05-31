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

    // class KMER{
        // private:
            bool is_fwd = true;
        // public:
            // Constructor for class
            KMER::KMER(bool is_fwd_) : is_fwd(is_fwd_){}

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
                // ako imamo bazu ACG reverse complement nastane kao ACG -> GCA -> CGT!!
                // 1. reversamo sekvencu
                string rc(kmer.rbegin(), kmer.rend());
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
                            unsigned int mapp_bit;
                            if(is_fwd){ mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                            } else{     mapp_bit = MappSeqCharPointerToBit(sequence+sequence_len-kmer_len-i, kmer_len);
                            }
                            end_window.push_back(make_tuple(mapp_bit, i+1, is_fwd));
                        
                        }
                        minimizers.push_back(GetTupleWithMinFirst(end_window));
                   
                }

                // NORMAL WINDOWS - najbolje za sad
                for (unsigned int i = 0; i <= sequence_len-kmer_len; i++){
                    if(i>=window_len){ window.pop_front();}

                    unsigned int mapp_bit;
                    if(is_fwd){ mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                    } else{     mapp_bit = MappSeqCharPointerToBit(sequence+sequence_len-kmer_len-i, kmer_len);
                    }
                    
                    window.push_back(make_tuple(mapp_bit, i+1, is_fwd));

                    if(i>=(window_len-1)){
                        minimizers.push_back(GetTupleWithMinFirst(window));
                    }
                }

                // KRAJ -> end-minimizer:
                for(unsigned int u=kmer_len; u<(window_len+kmer_len-1);u++){
                    if (sequence_len < u) break;
                    deque<tuple<unsigned int, unsigned int, bool>> end_window;
                    // deque<tuple<unsigned int, unsigned int, bool>> end_window_rev;
                    unsigned int start = sequence_len - u;

                    for (unsigned int i = start; i <= sequence_len - kmer_len; ++i){
                        unsigned int mapp_bit;
                        if(is_fwd){ mapp_bit = MappSeqCharPointerToBit(sequence+i, kmer_len);
                        } else{     mapp_bit = MappSeqCharPointerToBit(sequence+sequence_len-kmer_len-i, kmer_len);
                        }

                            end_window.push_back(make_tuple(mapp_bit, i+1, is_fwd)); 
                    }
                    minimizers.push_back(GetTupleWithMinFirst(end_window));
                }
                return minimizers;
            }

    // };
}

// FUCTIONS THAT ARE NO LONGER NEEDED:

// unsigned int MappKmerStringToBit(const string& kmer){ // Mapping function for a k-mer
//         unsigned int mapp = 0;
//         // Value mapping for positions (original order)
//         unordered_map<char, unsigned int> base_value = {
//             {'C', 0},
//             {'A', 1},
//             {'T', 2},
//             {'G', 3}
//         };
//         for (int i=0; i<kmer.size(); i++){
//             mapp<<=2;
//             mapp |= base_value[kmer[i]];
//         }
//         return mapp;
//     }

// unsigned int ReverseComplementSeqCharPointerToBit(const char* seq, unsigned int seq_len) {  // Reverse complement of a DNA string
//         unsigned int rc = 0;
//         // Value mapping for positions (original order)
//         unordered_map<char, unsigned int> base_value = {
//             {'C', 3},
//             {'A', 2},
//             {'T', 1},
//             {'G', 0}
//         };
//         for (unsigned int i=seq_len-1; i>0; i--){
//             rc<<=2;
//             rc |= base_value[seq[i]];
//             if(i==1){
//                 rc<<=2;
//                 rc |= base_value[seq[i-1]];
//             }
//         }
//         return rc;
//     }

// string MappKmerStringToString(const string& kmer){ // Mapping function for a k-mer - String
//         string mapp = "";
//         // Value mapping for positions (original order)
//         unordered_map<char, unsigned int> base_value = {
//             {'C', '0'},
//             {'A', '1'},
//             {'T', '2'},
//             {'G', '3'}
//         };
//         for (int i=0; i<kmer.size(); i++){
//             unsigned int val;
//             val = base_value[kmer[i]];
//             mapp += val;
//         }
//         return mapp;
//     }

// unsigned int ExtractKmer(unsigned int packed_seq, int i, int k) { // Find kmer inside of number that holds whole sequence, with 2 bits for each sequence
//         unsigned int shift = 2 * i; 
//         unsigned int mask = (1ULL << (2 * k)) - 1; // k*2 bits set to 1
//         return (packed_seq >> shift) & mask;
//     }
// unsigned int ExtractKmerWithMask(unsigned int packed_seq, int i, unsigned int mask) {
//         unsigned int shift = 2 * i; 
//         return (packed_seq >> shift) & mask;
//     }





// KORISNI PRINTOVI I KODOVI

// cout<<"unsigned int value rc: "<<mapp_seq_rev_bit_2<<endl;
// cout<<"tocna unsigned int value rc: "<<mapp_seq_rev_bit<<endl;
// cout<<"STRING SEKVENCE: "<<endl;
// cout<<"Seq: "<<sequence_str<<endl;
// cout<<"Rc_seq: "<<rc_sequence_str<<endl;
// cout<<"BIT SEKVENCE: "<<endl;
// cout<<"prva: "<<endl;
// cout<<mapp_seq_fwd_bit<<endl;
// cout<<MappKmerBitToString(mapp_seq_fwd_bit, sequence_len)<<endl;
// cout<<"prva NOVI POSTUPAK: "<<endl;
// cout<<mapp_seq_fwd_bit_2<<endl;
// cout<<MappKmerBitToString(mapp_seq_fwd_bit_2, sequence_len)<<endl;
// cout<<"druga: "<<endl;
// cout<<mapp_seq_rev_bit<<endl;
// cout<<MappKmerBitToString(mapp_seq_rev_bit, sequence_len)<<endl;
// cout<<"druga NOVI POSTUPAK"<<endl;
// cout<<mapp_seq_rev_bit_2<<endl;
// cout<<MappKmerBitToString(mapp_seq_rev_bit_2, sequence_len)<<endl;


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