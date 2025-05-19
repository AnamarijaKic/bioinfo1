#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <random>
#include <ctime>  
#include "bioparser/include/bioparser/fasta_parser.hpp"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "team_alignment/team_alignment.hpp"  

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



int main(int argc, char* argv[]) {
    team::AlignmentType align_type_;
    int match = 1, mismatch = -1, gap = -1;
    string file1, file2;

    if (argc < 2) {
        std::cerr << "Error: Not enough argumen"<<endl;
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
        isFasta = false;
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

        mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
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

        // cout << "ZA DEBUGIRANJE indexi 19438 i 6323 *************************************************"<<endl;
        // index1 = 19438;
        // index2 = 6323; //27


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
        cout<<score<<endl;


    }
    
    
    
    
    
    


    return 0;
}

