#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <random>
#include <ctime>  
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <set>
#include <queue>
#include <unordered_map>
#include "bioparser/include/bioparser/fasta_parser.hpp"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "team_alignment/team_alignment.hpp" 
#include "team_minimizers/team_minimizers.hpp"
#include "team_minimizers/team_minimizers.hpp"

#define VERSION "0.1.0"
#define PROGRAM_NAME "toolForGenomeAllignment"

using namespace std;
using namespace team;

//Marta Kekić
std::vector<std::tuple<unsigned int, unsigned int, bool>> remove_duplicates(
    const std::vector<std::tuple<unsigned int, unsigned int, bool>>& input) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result;
    std::unordered_set<std::string> seen;

    for (const auto& tup : input) {
        //serialize the tuple into a string to use as a unique key
        auto [hash, pos, strand] = tup;
        std::string key = std::to_string(hash) + "_" + std::to_string(pos) + "_" + std::to_string(strand);
        
        //if this key hasn't been seen before, add the tuple to the result
        if (seen.insert(key).second) {
            result.push_back(tup); 
        }
    }
    return result;
}

void PrintReferenceIndex(const std::unordered_map<unsigned int, std::set<std::pair<unsigned int, bool>>>& reference_index) {
    for (const auto& [minimizer_hash, positions_set] : reference_index) {
        std::cout << "Minimizer hash: " << minimizer_hash << "\n";
        for (const auto& [position, is_original_strand] : positions_set) {
            std::cout << "  Position: " << position
                      << ", Strand: " << (is_original_strand ? "+" : "-") << "\n";
        }
    }
}

//Anamarija Kic
// Reverse complement of a DNA string
    string ReverseComplement(const string& kmer) { 
        //ACG -> GCA -> CGT!!
        // 1. reverse seq
        string rc(kmer.rbegin(), kmer.rend()); //https://cplusplus.com/reference/string/string/rbegin/
        // 2. compl seq
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

// Mapping from bit shifted values to string
string MappKmerBitToStringFWD(unsigned int kmer, unsigned int kmer_len){
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

void PrintMinimizersVector(const vector<tuple<unsigned int, unsigned int, bool>>& minimizers, unsigned int k) {
    for (const auto& t : minimizers) {
        unsigned int hash = get<0>(t);
        bool is_fwd = get<2>(t);
        string mapp = MappKmerBitToStringFWD(hash, k);
        // if(is_fwd){
        //     mapp = MappKmerBitToStringFWD(hash, k);
        // } else{
        //     mapp = MappKmerBitToStringREV(hash, k);
        // }
        
        unsigned int position = get<1>(t);
        
        cout << "Hash: " << mapp
             << ", Position: " << position
             << ", Strand: " << (is_fwd ? "FWD" : "REV") << endl;
    }
}

void PrintMinimizersSet(const std::set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
    for (const auto& t : minimizers) {
        unsigned int hash = std::get<0>(t);
        unsigned int position = std::get<1>(t);
        bool is_fwd = std::get<2>(t);

        std::cout << "Hash: " << hash
                  << ", Position: " << position
                  << ", Strand: " << (is_fwd ? "FWD" : "REV") << std::endl;
    }
}

struct SequenceFASTA {  // or any other name
    public:
        SequenceFASTA(  // required arguments
            const char* name, std::uint32_t name_len,
            const char* data, std::uint32_t data_len) {
                name_ = string(name,name_len);
                data_ = string(data, data_len);
        }
        const string& name() const {return name_;}
        const string& data() const {return data_;}
    private:
        string name_;
        string data_;
};

struct SequenceFASTQ {  // or any other name
    public:
        SequenceFASTQ(  // required arguments
            const char* name, uint32_t name_len,
            const char* sequence, uint32_t sequence_len,
            const char* quality, uint32_t quality_len) {
                name_ = string(name,name_len);
                sequence_ = string(sequence, sequence_len);
                quality_ = string(quality, quality_len);
        }
        const string& name() const {return name_;}
        const string& sequence() const {return sequence_;}
    private:
        string name_;
        string sequence_;
        string quality_;
};

void printHelp(){
    cout<<endl;
    cout<<"Usage: "<<PROGRAM_NAME<<"[options] <file1> <file2>"<<endl;
    cout<<"NOTE: file1 needs to be in FASTA format, while the second file will contain a set of fragments in either FASTA or FASTQ format."<<endl;
    cout<<"Options: \n"
        << "\t  -a, --alignment TYPE     Alignment type: global, local, semiGlobal\n"
        << "\t  -m MATCH                 Match score (default: 1)\n"
        << "\t  -n MISMATCH              Mismatch penalty (default: -1)\n"
        << "\t  -g GAP                   Gap penalty (default: -1)\n"
        << "\t  -k KMER                  k-mer length for minimizers (default: 15)\n"
        << "\t  -w WINDOW                window size for minimizers (default: 5)\n"
        << "\t  -f FREQUENCY_THRESHOLD   Frequency threshold factor (default: 0.001)\n"
        << "\t  -c                       Output CIGAR string\n"
        << "\t  -h, --help               Show this help message\n"
        << "\t  --version                Show version information\n";
}

void printBasicStatisticFASTA(string file){
    auto p = bioparser::Parser<SequenceFASTA>::Create<bioparser::FastaParser>(file);
    // while(true){
    //     // parse whole file
    auto s = p->Parse(-1);
    // if (s.empty()) break;

    size_t totalNumberOfBases = 0;
    size_t totalNumberOfSequences = 0;
    size_t maxLength = 0;
    size_t minLength = SIZE_MAX;
    vector<size_t> allLengths;
    for (const auto& seq : s){
        size_t size = seq->data().size();

        allLengths.push_back(size);
        totalNumberOfBases += size;
        totalNumberOfSequences +=1;

        cout<<"SequenceFASTA name: "<<seq->name()<<endl;
        cout<<"Length of sequence: "<<size<<endl;

        if (maxLength<size){ maxLength = size; }
        if (minLength>size){ minLength = size; }            
    }
        
    cout << "Total number of sequences: " << totalNumberOfSequences << endl;
    cout << "Average length of sequences: " << totalNumberOfBases/totalNumberOfSequences << endl;
    cout << "Maximal length of sequence: " << maxLength <<endl;
    cout << "Minimal length of sequence: " << minLength <<endl;
    // std::cout << "Total bases: " << totalNumberOfBases << std::endl;

    // N50 value
    size_t cumulativeSum = 0;
    sort(allLengths.begin(), allLengths.end());
    for (size_t currentLength : allLengths){
        cumulativeSum += currentLength;
        if (cumulativeSum >= (totalNumberOfBases/2)){
            cout << "N50 length: " << currentLength <<endl;
            break;
        }
    }

        
    // }
}

//Marta Kekić
void printBasicStatisticFASTQ(string file){
    auto p = bioparser::Parser<SequenceFASTQ>::Create<bioparser::FastqParser>(file);
    
    // parse in chunks
    std::vector<std::unique_ptr<SequenceFASTQ>> s;
    while (true) {
    auto c = p->Parse(1ULL << 30);  // 1 GB
    if (c.empty()) {
        break;
    }
    s.insert(
        s.end(),
        std::make_move_iterator(c.begin()),
        std::make_move_iterator(c.end()));
    }

    size_t totalNumberOfBases = 0;
    size_t totalNumberOfSequences = 0;
    size_t maxLength = 0;
    size_t minLength = SIZE_MAX;
    vector<size_t> allLengths;
    for (const auto& seq : s){
        size_t size = seq->sequence().size();

        allLengths.push_back(size);
        totalNumberOfBases += size;
        totalNumberOfSequences +=1;

        cout<<"SequenceFASTQ name: "<<seq->name()<<endl;
        cout<<"Length of sequence: "<<size<<endl;

        if (maxLength<size){ maxLength = size; }
        if (minLength>size){ minLength = size; }            
    }
        
    cout << "Total number of sequences: " << totalNumberOfSequences << endl;
    cout << "Average length of sequences: " << totalNumberOfBases/totalNumberOfSequences << endl;
    cout << "Maximal length of sequence: " << maxLength <<endl;
    cout << "Minimal length of sequence: " << minLength <<endl;
    // std::cout << "Total bases: " << totalNumberOfBases << std::endl;

    // N50 value
    size_t cumulativeSum = 0;
    sort(allLengths.begin(), allLengths.end());
    for (size_t currentLength : allLengths){
        cumulativeSum += currentLength;
        if (cumulativeSum >= (totalNumberOfBases/2)){
            cout << "N50 length: " << currentLength <<endl;
            break;
        }
    }

        
    // }
}

//Longest Increasing Subsequence
vector<pair<unsigned int, unsigned int>> FindLIS(const vector<pair<unsigned int, unsigned int>>& matches) {
    size_t n = matches.size();
    if (n == 0) return {};

    vector<int> lis(n, 1);    //lis[i] stores the length of the LIS ending at position i
    vector<int> prev(n, -1); //prev[i] stores the index of the previous element in the LIS ending at i

    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            //check strictly increasing order on the reference positions (second)
            //and disallow multiple fragment positions (first) from matching to the same one
            if (matches[i].second > matches[j].second && lis[i] < lis[j] + 1 && matches[i].first != matches[j].first
            && matches[i].first - matches[j].first < 5000 && matches[i].second - matches[j].second < 5000) {
                lis[i] = lis[j] + 1;
                prev[i] = j;
            }
        }
    }

    //find the index of the maximum LIS length    
    int max_index = distance(lis.begin(), max_element(lis.begin(), lis.end()));

    //reconstruct the LIS by backtracking through the prev array
    vector<pair<unsigned int, unsigned int>> result;
    for (int i = max_index; i >= 0; i = prev[i]) {
        result.push_back(matches[i]);
        if (prev[i] == -1) break;
    }

    reverse(result.begin(), result.end()); //because we reconstructed it backwards

    return result;
}

void PrintLIS(const std::vector<std::pair<unsigned int, unsigned int>>& lis_matches) {
    std::cout << "Longest Increasing Subsequence (LIS) of matches:\n";
    for (const auto& match : lis_matches) {
        std::cout << "Fragment position: " << match.first
                  << ", Reference position: " << match.second << "\n";
    }
}


int main(int argc, char* argv[]) {
    team::AlignmentType align_type_ = team::AlignmentType::global;  // default;
    int match = 1, mismatch = -1, gap = -1;
    string file1, file2;
    unsigned int k = 15, w = 5;
    double f = 0.001;
    bool output_cigar = false;

    //Anamarija Kic
    if (argc < 2) {
        std::cerr << "Error: Not enough arguments"<<endl;
        printHelp();
        return 1;
    }

    string arg1 = argv[1];
    if (arg1 == "-h" || arg1 == "--help") {
        printHelp();
        return 0;
    } else if (arg1 == "--version") {
        std::cout << PROGRAM_NAME <<" v" << VERSION << endl;
        return 0;
    }

    if (argc < 3) {
        std::cerr << "Error: Expected two input files\n";
        return 1;
    }

    for (int i=1; i<argc; ++i){
        if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            int number = ++i;
            if (strcmp(argv[number], "global") == 0) { align_type_ = team::AlignmentType::global;
            } else if (strcmp(argv[number], "local") == 0) { align_type_ = team::AlignmentType::local;
            } else if (strcmp(argv[number], "semiGlobal") == 0) { align_type_ = team::AlignmentType::semiGlobal;
            } else { std::cerr << "Error: Expected Alignment type: global, local, semiGlobal\n";
            }
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            match = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            mismatch = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            gap = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k = static_cast<unsigned int>(std::atoi(argv[++i]));
        } else if (strcmp(argv[i], "-w") == 0 && i + 1 < argc) {
            w = static_cast<unsigned int>(std::atoi(argv[++i]));
        } else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            f = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-c") == 0) {
            output_cigar = true;
        }else if (file1.empty()) {
            file1 = argv[i];
        } else if (file2.empty()) {
            file2 = argv[i];
        } else {
            std::cerr << "Unknown or extra argument: " << argv[i] << "\n";
            printHelp();
            return 1;
        }

    }
    if (file1.empty() || file2.empty()) {
        std::cerr << "Error: Two input files are required.\n";
        printHelp();
        return 1;
    }

    // Output for debug
    std::cout << "Files: " << file1 << ", " << file2 << "\n";
    std::cout << "Align type: **************" 
              << ", Match: " << match
              << ", Mismatch: " << mismatch
              << ", Gap: " << gap << "\n";
    // std::cout << "Align type: " << align_type_;

    cout << "Processing files: " << file1 << " and " << file2 << endl;

    // The first file will contain a reference genome in FASTA format
    // --------------------------------------------------------------
    cout << "Basic statistic for reference genome"<<endl;
    cout << "------------------------------------"<<endl;
    auto referenceParser = bioparser::Parser<SequenceFASTA>::Create<bioparser::FastaParser>(file1);
    // parse whole file
    auto referenceSequence = referenceParser->Parse(-1);
    printBasicStatisticFASTA(file1);

    //Marta Kekić
    //map for minimizer index in the reference genome (key-minimizer hash, (position, orientation))
    unordered_map<unsigned int, set<pair<unsigned int, bool>>> reference_index_fwd;
    unordered_map<unsigned int, int> minimizer_frequencies_fwd;
    unordered_map<unsigned int, set<pair<unsigned int, bool>>> reference_index_rev;
    unordered_map<unsigned int, int> minimizer_frequencies_rev;

    auto& reference = referenceSequence.front()->data();
    //minimizers in the reference
    KMER ref(true);
    auto minimizers_fwd =ref.Minimize(reference.c_str(), reference.length(), k, w);
    //PrintMinimizersVector(minimizers_fwd, k);

    //complement of reference
    auto& reference_rev = ReverseComplement(reference);
    //minimizers in the complement reference
    KMER ref_rev(false);
    auto minimizers_rev = ref_rev.Minimize(reference_rev.c_str(), reference_rev.length(), k, w);
    PrintMinimizersVector(minimizers_rev, k);


    cerr << "DEBUG: Number of minimizers in ref fwd before f = " << minimizers_fwd.size() << endl;
    cerr << "DEBUG: Number of minimizers in ref rev before f = " << minimizers_rev.size() << endl;

    //how frequent is the minimizer in the reference
    for (const auto& [hash, pos, strand] : minimizers_fwd) {
        minimizer_frequencies_fwd[hash]++;
    }
    for (const auto& [hash, pos, strand] : minimizers_rev) {
        minimizer_frequencies_rev[hash]++;
    }

    int frequency_threshold = static_cast<int>(f * reference.length());
    //filtering minimizers that show up too often
    for (const auto& [hash, pos, strand] : minimizers_fwd) {
        if (minimizer_frequencies_fwd[hash] <= frequency_threshold) {
            //reference_index[hash].emplace_back(pos, strand);
            reference_index_fwd[hash].insert({pos, strand});
            //cerr << "hash: " << hash << " pos: " << pos << " strand: " << strand << endl;
        }
    }

    //filtering minimizers that show up too often
    for (const auto& [hash, pos, strand] : minimizers_rev) {
        if (minimizer_frequencies_rev[hash] <= frequency_threshold) {
            //reference_index[hash].emplace_back(pos, strand);
            reference_index_rev[hash].insert({pos, strand});
            //cerr << "hash: " << hash << " pos: " << pos << " strand: " << strand << endl;
        }
    }

    

    //Anamarija Kic
    // while the second file will contain a set of fragments in either FASTA or FASTQ format.
    // --------------------------------------------------------------------------------------
    cout << endl << "Basic statistic for fragments of genome"<<endl;
    cout << "------------------------------------"<<endl;
    bool isFastq = true;
    bool isFasta = true;
    vector<unique_ptr<SequenceFASTA>> fragmentSequencesFASTA;
    vector<std::unique_ptr<SequenceFASTQ>> fragmentSequencesFASTQ;
    try{
        auto p = bioparser::Parser<SequenceFASTQ>::Create<bioparser::FastqParser>(file2);
        // parse in chunks
        while (true) {
            auto c = p->Parse(1ULL << 30);  // 1 GB
            if (c.empty()) {
                break;
            }
            fragmentSequencesFASTQ.insert(
                fragmentSequencesFASTQ.end(),
                std::make_move_iterator(c.begin()),
                std::make_move_iterator(c.end()));
        }
    } catch (const std::exception& e) {
        cerr << "Failed to parse as FASTQ: " << e.what() << endl;
        isFastq = false;
        try{
            auto fragmentParser = bioparser::Parser<SequenceFASTA>::Create<bioparser::FastaParser>(file2);
            // parse whole file
            fragmentSequencesFASTA = fragmentParser->Parse(-1);
        } catch (const std::exception& e) {
                isFasta = false;
                cerr << "Failed to parse as FASTA: " << e.what() << endl;
        }
    }
    if (isFasta == false && isFastq == false){ return 0; } // end of program

    if(isFasta){
        cout<< "ODKOMENTIRAJ BASIC STATISTIC PO POTREBI ***********************************"<<endl;
        // printBasicStatisticFASTA(file2);

        // MINIMIZER ZA MALI PRIMJER ********************************************
        // vector<tuple<unsigned int, unsigned int, bool>> minimizers = team::Minimize(fragmentSequencesFASTA[0]->data().c_str(), fragmentSequencesFASTA[0]->data().size(), k, w);
        // PrintMinimizersVector(minimizers, k);

        if (fragmentSequencesFASTA.size()<2){  std::cerr << "Need at least 2 sequences to pick two different ones.\n"; return 0;}

        /*mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
        uniform_int_distribution<> dist(0, fragmentSequencesFASTA.size() - 1);

        int index1;
        do {
            index1 = rand() % fragmentSequencesFASTA.size();
        } while (fragmentSequencesFASTA[index1]->data().size()>5000);

        int index2;
        do {
            index2 = rand() % fragmentSequencesFASTA.size();
        } while ((fragmentSequencesFASTA[index1]->data().size()>5000)
                    || (index1==index2));

        //cout << "ZA DEBUGIRANJE indexi 19438 i 6323 *************************************************"<<endl;
        //index1 = 19438;
        //index2 = 6323; //27


        auto& seq1 = fragmentSequencesFASTA[index1];
        auto& seq2 = fragmentSequencesFASTA[index2];

        cout << "First random index: " << index1 << " -> " << seq1->name() << "\n";
        cout << "Second random index: " << index2 << " -> " << seq2->name() << "\n";

       
        // std::string cigar;
        // unsigned int target_begin;
        // // trebam svaku od ovih seq1 i seq2 pripasati k pravoj referenci***********************************************************************************************************+
        // int score = team::Align(
        //     seq1->data().c_str(), seq1->data().length(),
        //     seq2->data().c_str(), seq2->data().length(),
        //     align_type_,  // Note: enum must also be qualified
        //     match, mismatch, gap,
        //     &cigar, &target_begin
        // );

        cout<<seq1->data()<<endl;
        // for(int i=0; i<target_begin;i++){cout<<" ";}
        cout<<seq2->data()<<endl;
        cout<<cigar<<endl;
        cout<<score<<endl;*/

        //Marta Kekić
        cerr << "DEBUG: Number of fragment sequences = " << fragmentSequencesFASTA.size() << endl;
        cerr << "DEBUG: Reference length = " << reference.length() << endl;
        cerr << "DEBUG: Reference index fwd size = " << reference_index_fwd.size() << endl;
        cerr << "DEBUG: Reference index rev size = " << reference_index_rev.size() << endl;

        //PrintReferenceIndex(reference_index_fwd);



        for (const auto& seq : fragmentSequencesFASTA) {
            cout << "pocinje" << endl;
            //minimizers in the seq fragment - returns vector tuple(hash, position, strand)
            KMER frag(true);
            auto raw_frag_min = frag.Minimize(seq->data().c_str(), seq->data().length(), k, w);
            auto frag_min = remove_duplicates(raw_frag_min);

            PrintMinimizersVector(frag_min, k);

            cerr << "DEBUG: Number of minimizers in fragment = " << frag_min.size() << endl;
            cerr << "DEBUG: Number of minimizers in reference fwd = " << reference_index_fwd.size() << endl;
            cerr << "DEBUG: Number of minimizers in reference rev = " << reference_index_rev.size() << endl;


            //finding matches in the reference genome for the minimizers in frag_min

            //*****************************************************************MATCHES KAD SE SALJU SAMO MANJI MINIMIZERI */
            /*for (const auto& [hash, f_pos, f_strand] : frag_min) {
                if (reference_index.count(hash)) {
                    for (const auto& [r_pos, r_strand] : reference_index[hash]) {
                        //saving (position in the fragnment, position in the reference)
                        matches.emplace_back(f_pos, r_pos);
                    }
                }
            }*/

            vector<pair<unsigned int, unsigned int>> matches_fwd;
            vector<pair<unsigned int, unsigned int>> matches_rev;
            for (const auto& [hash, f_pos, f_strand] : frag_min) {
                if (reference_index_fwd.count(hash)) {
                    for (const auto& [r_pos, r_strand] : reference_index_fwd[hash]) {
                            matches_fwd.emplace_back(f_pos, r_pos);
                        }
                    for (const auto& [r_pos, r_strand] : reference_index_rev[hash]) {
                            matches_rev.emplace_back(f_pos, r_pos);   
                    }
                }
            }

            //matches print
            std::cout << "Matches (fragment_pos, reference_pos):\n";
            for (const auto& [f_pos, r_pos] : matches_fwd) {
                std::cout << "(" << f_pos << ", " << r_pos << ")\n";
            }

            cerr << "DEBUG: Number of matches fwd found = " << matches_fwd.size() << endl;
            cerr << "DEBUG: Number of matches rev found = " << matches_rev.size() << endl;


            //Longest Increasing Subsequence

            /*auto chain = FindLIS(matches);
            PrintLIS(chain);
            cerr << "DEBUG: LIS size = " << chain.size() << endl;
            if (chain.empty()) continue;*/

            auto chain_fwd = FindLIS(matches_fwd);
            PrintLIS(chain_fwd);
            auto chain_rev = FindLIS(matches_rev);
            PrintLIS(chain_rev);
            vector<pair<unsigned int, unsigned int>> chain;
            if (chain_fwd.size() >= chain_rev.size()){
                chain = chain_fwd;
            }
            else{
                chain = chain_rev;
            }
            if (chain.empty()) continue;

            //beginning and end positions in the fragment and the reference (target)
            unsigned int q_begin = chain.front().first-1;
            unsigned int q_end = chain.back().first + k-2;
            unsigned int t_begin = chain.front().second-1;
            unsigned int t_end = chain.back().second + k-2;

            if (q_end > seq->data().size()) q_end = seq->data().size();
            if (t_end > reference.size()) t_end = reference.size();

            if (q_end <= q_begin || t_end <= t_begin) {
                cerr << "ERROR: Invalid range for alignment" << endl;
                continue;
            }

            cerr << "DEBUG: q_begin = " << q_begin << ", q_end = " << q_end
                << ", t_begin = " << t_begin << ", t_end = " << t_end << endl;


            string cigar;
            unsigned int ref_offset;
            int score;
            string strand;

            try{ 
                if(chain == chain_fwd){
                    strand = "+";
                    score = team::Align(
                        seq->data().c_str() + q_begin, q_end - q_begin + 1,
                        reference.c_str() + t_begin, t_end - t_begin + 1,
                        align_type_, match, mismatch, gap,
                        output_cigar ? &cigar : nullptr, &ref_offset);
                }
                else{
                    strand = "-";
                    score = team::Align(
                        seq->data().c_str() + q_begin, q_end - q_begin + 1,
                        reference_rev.c_str() + t_begin, t_end - t_begin + 1,
                        align_type_, match, mismatch, gap,
                        output_cigar ? &cigar : nullptr, &ref_offset);
                }
                
                
                } catch (const std::exception& e) {
                    cerr << "ERROR: Exception during Align: " << e.what() << endl;
                    continue;
                }

            cerr << "DEBUG: Alignment score = " << score << ", ref_offset = " << ref_offset << endl;


            cerr << "DEBUG: Reached output stage" << endl;
            cout << seq->name() << "\t" << seq->data().size() << "\t" << q_begin << "\t" << (q_end + 1)
                << "\t"<< strand <<"\t" << referenceSequence.front()->name() << "\t" << reference.length()
                << "\t" << (chain==chain_fwd ? t_begin : reference_rev.length() - t_end - 1) 
                << "\t" <<(chain==chain_fwd ?  (t_end + 1) : reference_rev.length() - t_begin) 
                << "\t" << score << "\t" << (q_end - q_begin + 1)
                << "\t60";
            if (output_cigar) {
                cout << "\tcg:Z:" << cigar;
            }
            cout << endl;
        }

        /* Anamarija Kic
        // Get unique minimizers and Number of occurance of unique minimizers:
        // std::set<std::tuple<unsigned int, unsigned int, bool>> unique_minimizers;
        unordered_map<unsigned int, size_t> minimizer_counts;
        cout<<"Iteration through second file!"<<endl;
        cout<<"FASTA SIZE: "<<fragmentSequencesFASTA.size()<<endl;
        std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> all_minimizers;
        for (size_t i = 0; i < fragmentSequencesFASTA.size(); ++i) {
            // cout<<"i: "<<i<<endl;
            const char* sequence = fragmentSequencesFASTA[i]->data().c_str();
            size_t length = fragmentSequencesFASTA[i]->data().size();
            
            auto minimizers = team::Minimize(sequence, length, k, w);
            all_minimizers.push_back(minimizers);

            // Count and append unique minimizers
            for (const auto& minimizer : minimizers) {
                unsigned int mapp = get<0>(minimizer);
                minimizer_counts[mapp]++;
            }
            PrintMinimizersVector(minimizers, k);
            cout<<endl;
            cout<<endl;
            
        }

        size_t num_distinct = minimizer_counts.size();
        cout<<"Number of distinct minimizers: "<<num_distinct<<endl;

        // Fraction of singletons
        size_t num_singletons = count_if(
            minimizer_counts.begin(), minimizer_counts.end(),
            [](const pair<const unsigned int, size_t>& entry) {
                return entry.second == 1;
            }
        );
        double singleton_fraction = 1.0 * num_singletons / num_distinct;
        cout<<"Fraction of singletons: "<<singleton_fraction<<endl;

        // Most frequent minimizer (excluding top f)
        priority_queue<pair<size_t, unsigned int>> frequencies;
        for (const auto& kv : minimizer_counts) {
            frequencies.emplace(kv.second, kv.first);  
        }

        // Poping the top `f` frequent minimizers
        for (int i = 0; i < f*frequencies.size() && !frequencies.empty(); ++i) {
            frequencies.pop();
        }

        if (!frequencies.empty()) {
            size_t freq = frequencies.top().first;
            cout << "Most frequent minimizer (excluding top " << f << "): count = " << freq << endl;
        } else {
            cout << "Not enough minimizers to exclude top " << f << "!" << endl;
        }*/
    }

    //Marta Kekić
    if(isFastq){
        printBasicStatisticFASTQ(file2);
        for (const auto& seq : fragmentSequencesFASTQ) {
            KMER frag(true);
            auto raw_frag_min = frag.Minimize(seq->sequence().c_str(), seq->sequence().length(), k, w);
            auto frag_min = remove_duplicates(raw_frag_min);

            vector<pair<unsigned int, unsigned int>> matches_fwd;
            vector<pair<unsigned int, unsigned int>> matches_rev;
            for (const auto& [hash, f_pos, f_strand] : frag_min) {
                if (reference_index_fwd.count(hash)) {
                    for (const auto& [r_pos, r_strand] : reference_index_fwd[hash]) {
                        matches_fwd.emplace_back(f_pos, r_pos);
                    }
                }
                if (reference_index_rev.count(hash)){
                    for (const auto& [r_pos, r_strand] : reference_index_rev[hash]) {
                        matches_rev.emplace_back(f_pos, r_pos);   
                    }
                }
            }

            auto chain_fwd = FindLIS(matches_fwd);
            auto chain_rev = FindLIS(matches_rev);
            vector<pair<unsigned int, unsigned int>> chain;
            if (chain_fwd.size() >= chain_rev.size()){
                chain = chain_fwd;
            }
            else{
                chain = chain_rev;
            }
            if (chain.empty()) continue;

            unsigned int q_begin = chain.front().first-1;
            unsigned int q_end = chain.back().first + k-2;
            unsigned int t_begin = chain.front().second-1;
            unsigned int t_end = chain.back().second + k-2;

            string cigar;
            unsigned int ref_offset;
            int score;
            string strand;

            try{ 
                if(chain == chain_fwd){
                    strand = "+";
                    score = team::Align(
                        seq->sequence().c_str() + q_begin, q_end - q_begin + 1,
                        reference.c_str() + t_begin, t_end - t_begin + 1,
                        align_type_, match, mismatch, gap,
                        output_cigar ? &cigar : nullptr, &ref_offset);
                }
                else{
                    strand = "-";
                    score = team::Align(
                        seq->sequence().c_str() + q_begin, q_end - q_begin + 1,
                        reference_rev.c_str() + t_begin, t_end - t_begin + 1,
                        align_type_, match, mismatch, gap,
                        output_cigar ? &cigar : nullptr, &ref_offset);
                }
                
                
                } catch (const std::exception& e) {
                    cerr << "ERROR: Exception during Align: " << e.what() << endl;
                    continue;
                }

            cerr << "DEBUG: Alignment score = " << score << ", ref_offset = " << ref_offset << endl;


            cerr << "DEBUG: Reached output stage" << endl;
            cout << seq->name() << "\t" << seq->sequence().size() << "\t" << q_begin << "\t" << (q_end + 1)
                << "\t"<< strand <<"\t" << referenceSequence.front()->name() << "\t" << reference.length()
                << "\t" << (chain==chain_fwd ? t_begin : reference_rev.length() - t_end - 1) 
                << "\t" <<(chain==chain_fwd ?  (t_end + 1) : reference_rev.length() - t_begin) 
                << "\t" << score << "\t" << (q_end - q_begin + 1)
                << "\t60";
            if (output_cigar) {
                cout << "\tcg:Z:" << cigar;
            }
            cout << endl;

        }

    }

    return 0;
}

