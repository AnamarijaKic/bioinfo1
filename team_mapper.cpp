#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <random>
#include <ctime>  
#include <unordered_map>
#include "bioparser/include/bioparser/fasta_parser.hpp"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "team_alignment/team_alignment.hpp" 
#include "team_minimizers/team_minimizers.hpp"

#define VERSION "0.1.0"
#define PROGRAM_NAME "toolForGenomeAllignment"

using namespace std;

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
//vector<pair<unsigned int, unsigned int>> FindLIS(const vector<pair<unsigned int, unsigned int>>& matches) {
    /*vector<unsigned int> dp, prev(matches.size());
    vector<size_t> lis;*/

    //looking for Increasing Subsequence using binary search O(nlog(n))
    //for (size_t i = 0; i < matches.size(); ++i) {
        //std::cerr << "dp size: " << dp.size() << ", i: " << i << "\n";

        /*auto it = lower_bound(dp.begin(), dp.end(), matches[i].second,
            [&](unsigned int a, unsigned int b) { return matches[dp[a]].second < b; });
        */
       /*  auto it = std::lower_bound(dp.begin(), dp.end(), i, 
        [&](unsigned int a, unsigned int b) { return matches[a].second < matches[b].second; });

        if (it == dp.end()) {
            prev[i] = dp.empty() ? -1 : dp.back();
            dp.push_back(i);
        } else {
            *it = i;
            prev[i] = (it == dp.begin()) ? -1 : dp[it - dp.begin() - 1];
        }
    }

    //chain reconstruction
    vector<pair<unsigned int, unsigned int>> result;
    for (int i = dp.empty() ? -1 : dp.back(); i >= 0; i = prev[i]) {
        result.push_back(matches[i]);
    }
    reverse(result.begin(), result.end());
    return result;
}*/

vector<pair<unsigned int, unsigned int>> FindLIS(const vector<pair<unsigned int, unsigned int>>& matches) {
    size_t n = matches.size();
    if (n == 0) return {};

    vector<int> lis(n, 1);   // LIS dužine, inicijalno svi 1
    vector<int> prev(n, -1); // prethodni element u LIS-u

    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (matches[i].second > matches[j].second && lis[i] < lis[j] + 1) {
                lis[i] = lis[j] + 1;
                prev[i] = j;
            }
        }
    }

    // Nađi indeks najveće vrijednosti u lis[]
    int max_index = distance(lis.begin(), max_element(lis.begin(), lis.end()));

    // Rekonstruiraj LIS slijedeći prev[]
    vector<pair<unsigned int, unsigned int>> result;
    for (int i = max_index; i >= 0; i = prev[i]) {
        result.push_back(matches[i]);
        if (prev[i] == -1) break;
    }

    reverse(result.begin(), result.end()); // jer je išao unazad

    return result;
}



int main(int argc, char* argv[]) {
    team::AlignmentType align_type_ = team::AlignmentType::global;  // default;
    int match = 1, mismatch = -1, gap = -1;
    string file1, file2;
    unsigned int k = 15, w = 5;
    double f = 0.001;
    bool output_cigar = false;


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

    //map for minimizer index in the reference genome (key-minimizer hash, (position, orientation))
    unordered_map<unsigned int, vector<pair<unsigned int, bool>>> reference_index;
    unordered_map<unsigned int, int> minimizer_frequencies;

    auto& reference = referenceSequence.front()->data();
    //minimizers in the reference
    auto minimizers = team::Minimize(reference.c_str(), reference.length(), k, w);

    cout << "Minimizers for reference: " << referenceSequence.front()->name() << "\n";
            for (const auto& [hash, pos, strand] : minimizers) {
                cout << "  hash: " << hash
                    << ", position: " << pos
                    << ", strand: " << (strand ? "+" : "-") << "\n";
            }


    //how frequent is the minimizer in the reference
    for (const auto& [hash, pos, strand] : minimizers) {
        minimizer_frequencies[hash]++;
    }

    int frequency_threshold = static_cast<int>(f * reference.length());
    //filtering minimizers that show up too often
    for (const auto& [hash, pos, strand] : minimizers) {
        if (minimizer_frequencies[hash] <= frequency_threshold) {
            reference_index[hash].emplace_back(pos, strand);
        }
    }

    

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

        if (fragmentSequencesFASTA.size()<2){  std::cerr << "Need at least 2 sequences to pick two different ones.\n";}

        /*mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
        uniform_int_distribution<> dist(0, fragmentSequencesFASTA.size() - 1);

        int index1;
        do {
            index1 = dist(rng);
        } while (fragmentSequencesFASTA[index1]->data().size()>5000);

        int index2;
        do {
            index2 = dist(rng);
        } while ((fragmentSequencesFASTA[index1]->data().size()>5000)
                    || (index1==index2));

        //cout << "ZA DEBUGIRANJE indexi 19438 i 6323 *************************************************"<<endl;
        //index1 = 19438;
        //index2 = 6323; //27


        auto& seq1 = fragmentSequencesFASTA[index1];
        auto& seq2 = fragmentSequencesFASTA[index2];

        cout << "First random index: " << index1 << " -> " << seq1->name() << "\n";
        cout << "Second random index: " << index2 << " -> " << seq2->name() << "\n";

       
        std::string cigar;
        unsigned int target_begin;
        // trebam svaku od ovih seq1 i seq2 pripasati k pravoj referenci***********************************************************************************************************+
        int score = team::Align(
            seq1->data().c_str(), seq1->data().length(),
            seq2->data().c_str(), seq2->data().length(),
            align_type_,  // Note: enum must also be qualified
            match, mismatch, gap,
            &cigar, &target_begin
        );

        cout<<seq1->data()<<endl;
        // for(int i=0; i<target_begin;i++){cout<<" ";}
        cout<<seq2->data()<<endl;
        cout<<cigar<<endl;
        cout<<score<<endl;*/

        cerr << "DEBUG: Number of fragment sequences = " << fragmentSequencesFASTA.size() << endl;
        cerr << "DEBUG: Reference length = " << reference.length() << endl;
        cerr << "DEBUG: Reference index size = " << reference_index.size() << endl;




        for (const auto& seq : fragmentSequencesFASTA) {
            cout << "pocinje" << endl;
            //minimizers in the seq fragment - returns vector tuple(hash, position, strand)
            auto frag_min = team::Minimize(seq->data().c_str(), seq->data().length(), k, w);

            cout << "Minimizers for fragment: " << seq->name() << "\n";
            for (const auto& [hash, pos, strand] : frag_min) {
                cout << "  hash: " << hash
                    << ", position: " << pos
                    << ", strand: " << (strand ? "+" : "-") << "\n";
            }

            cerr << "DEBUG: Number of minimizers in fragment = " << frag_min.size() << endl;

            vector<pair<unsigned int, unsigned int>> matches;
            //finding matches in the reference genom for the minimizers in frag_min 
            for (const auto& [hash, f_pos, f_strand] : frag_min) {
                if (reference_index.count(hash)) {
                    for (const auto& [r_pos, r_strand] : reference_index[hash]) {
                        //saving (position in the fragnment, position in the reference)
                        matches.emplace_back(f_pos, r_pos);
                    }
                }
            }

            cerr << "DEBUG: Number of matches found = " << matches.size() << endl;


            //Longest Increasing Subsequence
            auto chain = FindLIS(matches);
            cerr << "DEBUG: LIS size = " << chain.size() << endl;
            if (chain.empty()) continue;

            //beginning and end positions in the fragment and the reference (target)
            unsigned int q_begin = chain.front().first;
            unsigned int q_end = chain.back().first + k;
            unsigned int t_begin = chain.front().second;
            unsigned int t_end = chain.back().second + k;

            //cerr << "DEBUG: " << q_begin << q_end << t_begin << t_end << endl;

            // Ograniči q_end i t_end ako prelaze duljinu sekvence
            if (q_end > seq->data().size()) {
                cerr << "DEBUG: q_end prevelik (" << q_end << "), postavljam na " << seq->data().size() << endl;
                q_end = seq->data().size();
            }
            if (t_end > reference.size()) {
                cerr << "DEBUG: t_end prevelik (" << t_end << "), postavljam na " << reference.size() << endl;
                t_end = reference.size();
            }

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

            try{ 
                score = team::Align(
                    seq->data().c_str() + q_begin, q_end - q_begin,
                    reference.c_str() + t_begin, t_end - t_begin,
                    align_type_, match, mismatch, gap,
                    output_cigar ? &cigar : nullptr, &ref_offset);
                
                } catch (const std::exception& e) {
                    cerr << "ERROR: Exception during Align: " << e.what() << endl;
                    continue;
                }

            cerr << "DEBUG: Alignment score = " << score << ", ref_offset = " << ref_offset << endl;


            cerr << "DEBUG: Reached output stage" << endl;
            cout << seq->name() << "\t" << seq->data().size() << "\t" << q_begin << "\t" << q_end
                << "\t+\t" << referenceSequence.front()->name() << "\t" << reference.length()
                << "\t" << t_begin << "\t" << t_end << "\t" /*<< score << "\t"*/ << (q_end - q_begin)
                << "\t60";
            if (output_cigar) {
                cout << "\tcg:Z:" << cigar;
            }
            cerr << "DEBUG: Output completed" << endl;
            cout << endl;
        }


    }

    if(isFastq){
        printBasicStatisticFASTQ(file2);
        for (const auto& seq : fragmentSequencesFASTQ) {
            auto frag_min = team::Minimize(seq->sequence().c_str(), seq->sequence().length(), k, w);

            vector<pair<unsigned int, unsigned int>> matches;
            for (const auto& [hash, f_pos, f_strand] : frag_min) {
                if (reference_index.count(hash)) {
                    for (const auto& [r_pos, r_strand] : reference_index[hash]) {
                        matches.emplace_back(f_pos, r_pos);
                    }
                }
            }

            auto chain = FindLIS(matches);
            if (chain.empty()) continue;

            unsigned int q_begin = chain.front().first;
            unsigned int q_end = chain.back().first + k;
            unsigned int t_begin = chain.front().second;
            unsigned int t_end = chain.back().second + k;

            string cigar;
            unsigned int ref_offset;

            int score = team::Align(
                seq->sequence().c_str() + q_begin, q_end - q_begin,
                reference.c_str() + t_begin, t_end - t_begin,
                align_type_, match, mismatch, gap,
                output_cigar ? &cigar : nullptr, &ref_offset);

            cout << seq->name() << "\t" << seq->sequence().size() << "\t" << q_begin << "\t" << q_end
                << "\t+\t" << referenceSequence.front()->name() << "\t" << reference.length()
                << "\t" << t_begin << "\t" << t_end << "\t" << score << "\t" << (q_end - q_begin)
                << "\t60";
            if (output_cigar) {
                cout << "\tcg:Z:" << cigar;
            }
            cout << endl;
        }

    }
    


    return 0;
}

