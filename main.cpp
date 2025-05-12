#include <iostream>
#include <string>
#include "bioparser/include/bioparser/fasta_parser.hpp"

#define VERSION "0.1.0"
#define PROGRAM_NAME "toolForGenomeAllignment"

using namespace std;

void printHelp(){
    cout<<"Usage: "<<PROGRAM_NAME<<"[options] <file1> <file2>"<<endl;
    cout<<"NOTE: file1 needs to be in FASTA format, while the second file will contain a set of fragments in either FASTA or FASTQ format."<<endl;
    cout<<"Options: "<<endl;
    cout<<"\t-h or --help\t\t display the help"<<endl;
    cout<<"\t--version \t\t display the version of program"<<endl;

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


    return 0;
}

