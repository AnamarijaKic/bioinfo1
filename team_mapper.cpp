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
#include <omp.h>

#define VERSION "3.1.0"
#define PROGRAM_NAME "toolForGenomeAllignment"

using namespace std;
using namespace team;

// Marta Kekić
// Removes duplicates from vector of minimizers
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

// Anamarija Kic
// Reverse complement of a DNA string
    string ReverseComplement(const string& kmer) { 
        //ACG -> GCA -> CGT!!
        // 1. reverse sequence
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

// Mapping from bit shifted values to string for printing out values
string MappKmerBitToStringFWD(unsigned int kmer, unsigned int kmer_len){
    std::string mapp(kmer_len, 'X');
    // // Value mapping for positions (original order)
    // unordered_map<unsigned int, char> base_value = {
    //     {0,'C'},
    //     {1,'A'},
    //     {2,'T'},
    //     {3,'G'}
    // };
    // Value mapping for positions (original order)
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

// Marta Kekić
// Printing out minimizers from reference genome
void PrintReferenceIndex(const std::unordered_map<unsigned int, std::set<std::pair<unsigned int, bool>>>& reference_index, unsigned int kmer_len) {
    for (const auto& [minimizer_hash, positions_set] : reference_index) {
        std::cout << "Minimizer: " << MappKmerBitToStringFWD(minimizer_hash, kmer_len) << "\n";
        for (const auto& [position, is_original_strand] : positions_set) {
            std::cout << "  Position: " << position
                      << ", Strand: " << (is_original_strand ? "+" : "-") << "\n";
        }
    }
}

// Anamarija Kic
// Printing minimizers from KMER.Minimize function
void PrintMinimizersVector(const vector<tuple<unsigned int, unsigned int, bool>>& minimizers, unsigned int k) {
    for (const auto& t : minimizers) {
        unsigned int hash = get<0>(t);
        bool is_fwd = get<2>(t);
        string mapp = MappKmerBitToStringFWD(hash, k);
        unsigned int position = get<1>(t);
        
        cout << "Minimizer: " << mapp
             << ", Position: " << position
             << ", Strand: " << (is_fwd ? "FWD" : "REV") << endl;
    }
}

// Printing set of minimizers
void PrintMinimizersSet(const std::set<std::tuple<unsigned int, unsigned int, bool>>& minimizers, unsigned int kmer_len) {
    for (const auto& t : minimizers) {
        unsigned int hash = std::get<0>(t);
        unsigned int position = std::get<1>(t);
        bool is_fwd = std::get<2>(t);

        std::cout << "Minimizer: " << MappKmerBitToStringFWD(hash, kmer_len)
                  << ", Position: " << position
                  << ", Strand: " << (is_fwd ? "FWD" : "REV") << std::endl;
    }
}

// Structure needed for bioparsers reading from FASTA files
struct SequenceFASTA {  
    public:
        SequenceFASTA(  
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

// Structure needed for bioparsers reading from FASTQ files
struct SequenceFASTQ {  
    public:
        SequenceFASTQ(  
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

// Printing help for input
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
        << "\t  --version                Show version information\n"
        << "\t  -s                       Basic statistic for first and second file\n"
        ;
}

// Printing basic statistic for FASTA file
void printBasicStatisticFASTA(string file){
    auto p = bioparser::Parser<SequenceFASTA>::Create<bioparser::FastaParser>(file);
    auto s = p->Parse(-1);
    
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
        
    cout << "Total number of sequences: " << totalNumberOfSequences <<endl;
    cout << "Average length of sequences: " << totalNumberOfBases/totalNumberOfSequences <<endl;
    cout << "Maximal length of sequence: " << maxLength <<endl;
    cout << "Minimal length of sequence: " << minLength <<endl;
    // cout << "Total number of bases: " << totalNumberOfBases <<endl;

    // N50 value
    size_t cumulativeSum = 0;
    sort(allLengths.begin(), allLengths.end(), greater<size_t>());
    for (size_t currentLength : allLengths){
        cumulativeSum += currentLength;
        if (cumulativeSum >= (totalNumberOfBases/2)){
            cout << "N50 length: " << currentLength <<endl;
            break;
        }
    }
}

// Marta Kekić
// Printing basic statistic for FASTQ file
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
        
    cout << "Total number of sequences: " << totalNumberOfSequences <<endl;
    cout << "Average length of sequences: " << totalNumberOfBases/totalNumberOfSequences <<endl;
    cout << "Maximal length of sequence: " << maxLength <<endl;
    cout << "Minimal length of sequence: " << minLength <<endl;
    // cout << "Total number of bases: " << totalNumberOfBases <<endl;

    // N50 value
    size_t cumulativeSum = 0;
    sort(allLengths.begin(), allLengths.end(), greater<size_t>());
    for (size_t currentLength : allLengths){
        cumulativeSum += currentLength;
        if (cumulativeSum >= (totalNumberOfBases/2)){
            cout << "N50 length: " << currentLength <<endl;
            break;
        }
    }
}

//Longest Increasing Subsequence Problem Algorithm
vector<pair<unsigned int, unsigned int>> FindLIS(const vector<pair<unsigned int, unsigned int>>& matches) {
    size_t n = matches.size(); 
    if (n == 0) return {};

    vector<int> lis(n, 1);   // lis[i] stores the length of the LIS ending at position i
    vector<int> prev(n, -1); // prev[i] stores the index of the previous element in the LIS ending at i

    // matches.first - index from fragment genome, sorted ascending
    // matches.second - index from reference genome, checking strictly increasing order on the reference positions
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            // Disallow multiple fragment positions (first) from matching to the same one
            if (matches[i].second > matches[j].second && lis[i] < lis[j] + 1 && matches[i].first != matches[j].first
                && matches[i].first-matches[j].first<5000 && matches[i].second-matches[j].second<5000) {
                lis[i] = lis[j] + 1;
                prev[i] = j;
            }
        }
    }

    // Find the index of the maximum LIS length    
    int max_index = distance(lis.begin(), max_element(lis.begin(), lis.end()));

    // Reconstruct the LIS by backtracking through the prev array
    vector<pair<unsigned int, unsigned int>> result;
    for (int i = max_index; i >= 0; i = prev[i]) {
        result.push_back(matches[i]);
        if (prev[i] == -1) break;
    }

    reverse(result.begin(), result.end()); // Because we reconstructed it backwards

    return result;
}

// Printing found LIS as (fragment_position, reference_position)
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
    bool output_cigar = false, statistic = false;

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
            } else { std::cerr << "Error: Expected Alignment type: global, local, semiGlobal\n"; printHelp(); return 1;
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
        } else if (strcmp(argv[i], "-s") == 0) {
            statistic = true;
        } else if (file1.empty()) {
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

    // The first file will contain a reference genome in FASTA format
    // --------------------------------------------------------------
    auto referenceParser = bioparser::Parser<SequenceFASTA>::Create<bioparser::FastaParser>(file1);
    auto referenceSequence = referenceParser->Parse(-1); // parse whole file
    
    if(statistic){
        cout << "Basic statistic for reference genome"<<endl;
        cout << "------------------------------------"<<endl;
        printBasicStatisticFASTA(file1);
    }

    // Marta Kekić
    // Map for minimizer index in the reference genome (key-minimizer hash, (position, orientation))
    unordered_map<unsigned int, set<pair<unsigned int, bool>>> reference_index_fwd;
    unordered_map<unsigned int, set<pair<unsigned int, bool>>> reference_index_rev;

    auto& reference = referenceSequence.front()->data();
    //minimizers in the reference
    KMER ref(true);
    auto minimizers_fwd =ref.Minimize(reference.c_str(), reference.length(), k, w);

    unordered_map<unsigned int, int> minimizer_frequencies_fwd = ref.GetMinimizerFrequencies();
    // PrintMinimizersVector(minimizers_fwd, k);

    //complement of reference
    auto reference_rev = ReverseComplement(reference);
    //minimizers in the complement reference
    KMER ref_rev(false);
    auto minimizers_rev = ref_rev.Minimize(reference_rev.c_str(), reference_rev.length(), k, w);

    unordered_map<unsigned int, int> minimizer_frequencies_rev = ref_rev.GetMinimizerFrequencies();
    // PrintMinimizersVector(minimizers_rev, k);

    //filtering minimizers that show up too often
    int frequency_threshold_fwd = static_cast<int>(f * ref.GetUniqueMinimizers().size());
    int frequency_threshold_rev = static_cast<int>(f * ref_rev.GetUniqueMinimizers().size());

    // creating a vector intead of a map sto that we can sort it
    std::vector<std::pair<unsigned int, int>> freq_vec_fwd(minimizer_frequencies_fwd.begin(), minimizer_frequencies_fwd.end());
    std::vector<std::pair<unsigned int, int>> freq_vec_rev(minimizer_frequencies_rev.begin(), minimizer_frequencies_rev.end());

    // sorting minimizers by frequency
    std::sort(freq_vec_fwd.begin(), freq_vec_fwd.end(),
        [](const auto& a, const auto& b) {
            return a.second > b.second; // descending
        });

    // banned_hashes_fwd => frequency_threshold_fwd most common minimizers
    unordered_set<unsigned int> banned_hashes_fwd;
    for (int i = 0; i < std::min(frequency_threshold_fwd, static_cast<int>(freq_vec_fwd.size())); ++i) {
        banned_hashes_fwd.insert(freq_vec_fwd[i].first);
    }

    // reference_index_fwd => only minimizers not in banned_hashes_fwd
    for (const auto& [hash_fwd, pos_fwd, strand_fwd] : minimizers_fwd) {
        if (!banned_hashes_fwd.count(hash_fwd)) {
            reference_index_fwd[hash_fwd].insert({pos_fwd, strand_fwd});
        }
    }

    // for reverse complement
    // sorting minimizers by frequency
    std::sort(freq_vec_rev.begin(), freq_vec_rev.end(),
        [](const auto& a, const auto& b) {
            return a.second > b.second; // descending
        });

    // banned_hashes_rev => frequency_threshold_rev most common minimizers
    unordered_set<unsigned int> banned_hashes_rev;
    for (int i = 0; i < std::min(frequency_threshold_rev, static_cast<int>(freq_vec_rev.size())); ++i) {
        banned_hashes_rev.insert(freq_vec_fwd[i].first);
    }

    // reference_index_rev => only minimizers not in banned_hashes_rev
    for (const auto& [hash_rev, pos_rev, strand_rev] : minimizers_rev) {
        if (!banned_hashes_rev.count(hash_rev)) {
            reference_index_rev[hash_rev].insert({pos_rev, strand_rev});
        }
    }

    cout<< "DEBUG: GOTOV SAM S OVIM!"<<endl;

    // Anamarija Kic
    if(statistic){
        // Get unique minimizers and number of occurance of unique minimizers:
        cout<<"Number of distinct minimizers for forward strand: "<<minimizer_frequencies_fwd.size()<<endl;
        cout<<"Number of distinct minimizers for reverse complement: "<<minimizer_frequencies_rev.size()<<endl;

        // Fraction of singletons
        int count_fwd = 0;
        for (const auto& pair : minimizer_frequencies_fwd) {
            if (pair.second == 1) {
                ++count_fwd;
            }
        }
        int count_rev = 0;
        for (const auto& pair : minimizer_frequencies_rev) {
            if (pair.second == 1) {
                ++count_rev;
            }
        }
        double singleton_fraction_fwd = 1.0 * count_fwd / minimizer_frequencies_fwd.size();
        cout<<"Fraction of singletons on forward strand: "<<singleton_fraction_fwd<<endl;
        double singleton_fraction_rev = 1.0 * count_rev / minimizer_frequencies_rev.size();
        cout<<"Fraction of singletons on reverse complement: "<<singleton_fraction_rev<<endl;

        //  Number of occurrences of the most frequent minimizer when the top f frequent minimizers are not taken in account 
        unsigned int best_key_fwd = 0;
        int max_value_fwd = numeric_limits<int>::min();
        for (const auto& pair : minimizer_frequencies_fwd) {
            if (reference_index_fwd.find(pair.first) != reference_index_fwd.end() && pair.second > max_value_fwd) {
                max_value_fwd = pair.second;
                best_key_fwd = pair.first;
            }
        }
        unsigned int best_key_rev = 0;
        int max_value_rev = numeric_limits<int>::min();
        for (const auto& pair : minimizer_frequencies_rev) {
            if (reference_index_rev.find(pair.first) != reference_index_rev.end() && pair.second > max_value_rev) {
                max_value_rev = pair.second;
                best_key_rev = pair.first;
            }
        }
        if(max_value_fwd!=numeric_limits<int>::min()){ cout << "Minimizer on forward strand with max value: " << MappKmerBitToStringFWD(best_key_fwd,k) << ", Value: " << max_value_fwd <<endl;
        } else { cout << "There are no minimizeres on forward strand after removing (1-f) percent of most frequent minimizers."<<endl; }
        if(max_value_rev!=numeric_limits<int>::min()){ cout << "Minimizer on reverse complement with max value: " << MappKmerBitToStringFWD(best_key_rev,k) << ", Value: " << max_value_rev <<endl;
        } else { cout << "There are no minimizeres on forward strand after removing (1-f) percent of most frequent minimizers."<<endl; }
    }
    
    // while the second file will contain a set of fragments in either FASTA or FASTQ format -> so we need to determine format
    // ------------------------------------------------------------------------------------------------------------------------------
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
        // cerr << "Failed to parse as FASTQ: " << e.what() << endl;
        isFastq = false;
        try{
            auto fragmentParser = bioparser::Parser<SequenceFASTA>::Create<bioparser::FastaParser>(file2);
            fragmentSequencesFASTA = fragmentParser->Parse(-1); // parse whole file
        } catch (const std::exception& e) {
                isFasta = false;
                // cerr << "Failed to parse as FASTA: " << e.what() << endl;
        }
    }
    if (isFasta == false && isFastq == false){ cerr << "Given file is not in FASTA or FASTQ format! "<<endl; return 1; } // end of program

    if(isFasta){
        if(statistic){
            cout << endl << "Basic statistic for fragments of genome"<<endl;
            cout << "------------------------------------"<<endl;
            printBasicStatisticFASTA(file2);
        }
        
        if (fragmentSequencesFASTA.size()<2){  std::cerr << "Need at least 2 sequences to pick two different ones.\n"; return 0;}

        // Two random sequences from the second input file, with lengths not exceeding 5000 base pairs, need to be aligned and the result reported.
        // mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
        // uniform_int_distribution<> dist(0, fragmentSequencesFASTA.size() - 1);
        // int index1;
        // do { index1 = rand() % fragmentSequencesFASTA.size();
        // } while (fragmentSequencesFASTA[index1]->data().size()>5000);
        // int index2;
        // do { index2 = rand() % fragmentSequencesFASTA.size();
        // } while ((fragmentSequencesFASTA[index1]->data().size()>5000) || (index1==index2));
        // auto& seq1 = fragmentSequencesFASTA[index1];
        // auto& seq2 = fragmentSequencesFASTA[index2];
        // cout << "First random index: " << index1 << " -> " << seq1->name() << "\n";
        // cout << "Second random index: " << index2 << " -> " << seq2->name() << "\n";
        // std::string cigar;
        // unsigned int target_begin;
        // int score = team::Align(
        //     seq1->data().c_str(), seq1->data().length(),
        //     seq2->data().c_str(), seq2->data().length(),
        //     align_type_,  // Note: enum must also be qualified
        //     match, mismatch, gap,
        //     &cigar, &target_begin
        // );
        // cout<<seq1->name()<<endl;
        // cout<<seq2->name()<<endl;
        // cout<<cigar<<endl;
        // cout<<score<<endl;

        // Marta Kekić
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < fragmentSequencesFASTA.size(); ++i) {
            auto& seq = fragmentSequencesFASTA[i];
            KMER frag(true);
            //Anamarija Kic
            if(!statistic){ frag.SetFrequenciesCount(false);}
            else{ frag.SetFrequenciesCount(true);}

            // Marta Kekić
            //minimizers in the seq fragment - returns vector tuple(hash, position, strand)
            auto raw_frag_min = frag.Minimize(seq->data().c_str(), seq->data().length(), k, w);
            auto frag_min = remove_duplicates(raw_frag_min);

            // Anamarija Kic
            if(statistic){
                // Get unique minimizers and number of occurance of unique minimizers:
                auto minimizer_frequencies_fwd = frag.GetMinimizerFrequencies();
                cout<<"Number of distinct minimizers for forward strand: "<<minimizer_frequencies_fwd.size()<<endl;
            
                // Fraction of singletons
                int count_fwd = 0;
                for (const auto& pair : minimizer_frequencies_fwd) {
                    if (pair.second == 1) {
                        ++count_fwd;
                    }
                }
                double singleton_fraction_fwd = 1.0 * count_fwd / minimizer_frequencies_fwd.size();
                cout<<"Fraction of singletons on forward strand: "<<singleton_fraction_fwd<<endl;
            }

            // Marta Kekić
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
            auto chain_fwd = FindLIS(matches_fwd);
            //PrintLIS(chain_fwd);
            auto chain_rev = FindLIS(matches_rev);
            //PrintLIS(chain_rev);
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

            #pragma omp critical
            {
                std::cout << seq->name() << "\t" << seq->data().size() << "\t" << q_begin << "\t" << (q_end + 1)
                << "\t"<< strand <<"\t" << referenceSequence.front()->name() << "\t" << reference.length()
                << "\t" << (chain==chain_fwd ? t_begin : reference_rev.length() - t_end - 1) 
                << "\t" <<(chain==chain_fwd ?  (t_end + 1) : reference_rev.length() - t_begin) 
                << "\t" << score << "\t" << (q_end - q_begin + 1)
                << "\t60";
            
                if (output_cigar) {
                    std::cout << "\tcg:Z:" << cigar;
                }
                std::cout<<std::endl;
            }
            
        }        
    }

    if(isFastq){
        if(statistic){
            std::cout << std::endl << "Basic statistic for fragments of genome"<<std::endl;
            std::cout << "------------------------------------"<<std::endl;
            printBasicStatisticFASTQ(file2);
        }
        
        #pragma omp parallel for schedule(dynamic)
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


            #pragma omp critical
            {
                std::cout << seq->name() << "\t" << seq->sequence().size() << "\t" << q_begin << "\t" << (q_end + 1)
                    << "\t" << strand << "\t" << referenceSequence.front()->name() << "\t" << reference.length()
                    << "\t" << (chain==chain_fwd ? t_begin : reference_rev.length() - t_end - 1)  
                    << "\t" << (chain==chain_fwd ?  (t_end + 1) : reference_rev.length() - t_begin) 
                    << "\t" << score << "\t" << (q_end - q_begin + 1)
                    << "\t60";
                if (output_cigar) {
                    std::cout << "\tcg:Z:" << cigar;
                }
                std::cout << std::endl;
            }
        }

    }
    
    return 0;
}

