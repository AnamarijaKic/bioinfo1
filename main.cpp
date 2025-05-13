#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include "bioparser/include/bioparser/fasta_parser.hpp"
#include "bioparser/include/bioparser/fastq_parser.hpp"

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
    cout<<"Options: "<<endl;
    cout<<"\t-h or --help\t\t display the help"<<endl;
    cout<<"\t--version \t\t display the version of program"<<endl;
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

    string file1 = argv[1];
    string file2 = argv[2];

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
        printBasicStatisticFASTA(file2);
    }
    
    
    
    


    return 0;
}

