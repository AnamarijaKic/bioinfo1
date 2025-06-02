#include "team_alignment.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <algorithm>

#define MATCH 0  // symbol for match
#define INSERT 1 // symbol for insert 
#define DELETE 2 // symbol for delete 


// Anamarija Kic
namespace team {
    typedef struct {
        int cost; // cost of reaching this cell 
        int parent; // parent cell 
    } cell;

    int match_func(char a, char b, int match, int mismatch) {
        if (a == b) return match;   // No penalty for a match
        else return mismatch;       // Penalty for a mismatch
    };

    int indel(char c, int gap) {
        if (c == '-') return 0;     // gaps have no cost
        else return gap;            
    }

    void printAlignmentMatrix(const char* query, unsigned int query_len,
                                const char* target, unsigned int target_len, std::vector<std::vector<team::cell>>& matrix){
        std::cout<<"\n\t \t ";
        for (int j=0; j<target_len; j++){
            std::cout<<target[j]<<"\t ";
        }
        std::cout<<std::endl;

        for(int i=0; i<query_len+1; i++){
            if(i>0){std::cout<<query[i-1]<<"\t ";}
            else if(i==0){std::cout<<"\t ";}
            for (int j=0; j<target_len+1; j++){
                std::cout<<matrix[i][j].cost<<"\t ";
            }
            std::cout<<std::endl<<std::endl;
        }
    }


    int Align( const char* query, unsigned int query_len,
                const char* target, unsigned int target_len,
                team::AlignmentType type,
                int match,
                int mismatch,
                int gap,
                std::string* cigar,
                unsigned int* target_begin){

        int initialization;
        std::string query_str(query, query_len);
        std::string target_str(target, target_len);

        switch (type) {
            case AlignmentType::global:
                initialization = gap;
                break;
            case AlignmentType::local:
                initialization = 0;
                break;
            case AlignmentType::semiGlobal:
                initialization = 0;
                break;
            default:
                throw std::invalid_argument("Unknown AlignmentType provided.");
        }

        // Matrix initialization
        std::vector<std::vector<cell>> m(query_len + 1, std::vector<cell>(target_len + 1)); 
        int i,j,k; // counters 
        int gi=0, gj=0; // goal cell coordinates 
        int opt[3]; // cost of the three options 

        // Initialize first column (i.e., m[i][0])
        for (unsigned int i = 0; i <= query_len; ++i) {
            m[i][0].cost = i * initialization;            
            m[i][0].parent = DELETE;           
        }

        // Initialize first row (i.e., m[0][j])
        for (unsigned int j = 0; j <= target_len; ++j) {
            m[0][j].cost = j * initialization;          
            m[0][j].parent = INSERT;          
        }

        std::string result;
        std::string compressed;
        int max_cost=std::numeric_limits<int>::min();
        int global_i = query_len;
        int global_j = target_len;
        int local_i=0, local_j=0; 
        switch (type) {
            case AlignmentType::global:
                for (i=1; i<=query_len; i++){
                    for (j=1; j<=target_len; j++) {
                        opt[MATCH] = m[i-1][j-1].cost + match_func(query_str[i-1],target_str[j-1],match,mismatch);
                        opt[INSERT] = m[i][j-1].cost + indel(target_str[j-1], gap);
                        opt[DELETE] = m[i-1][j].cost + indel(query_str[i-1], gap);
                        m[i][j].cost = opt[0];
                        m[i][j].parent = 0;
                        for (k=1; k<=2; k++) {
                            if (opt[k] > m[i][j].cost) {
                                m[i][j].cost = opt[k];
                                m[i][j].parent = k;
                            }
                        }
                    }   
                }
                gi = query_len;
                gj = target_len;
                if (target_begin != nullptr) {
                    *target_begin = 0;
                }
                if (cigar != nullptr) {
                    while (global_i > 0 || global_j > 0) {
                        if (global_i > 0 && global_j > 0 && m[global_i][global_j].parent == MATCH) {
                            result.push_back('M');
                            global_i--;
                            global_j--;
                        } else if (global_j > 0 && m[global_i][global_j].parent == INSERT) {
                            result.push_back('I');
                            global_j--;
                        } else if (global_i > 0 && m[global_i][global_j].parent == DELETE) {
                            result.push_back('D');
                            global_i--;
                        } else {
                            throw std::invalid_argument("Unknown error in determining cigar string.");
                            break;
                        }
                    }

                    
                    std::reverse(result.begin(), result.end());
                    
                    //*cigar = result;

                    char prev = result[0];
                    int count = 1;
                    // Compressing cigar for better output
                    for (size_t i = 1; i < result.size(); ++i) {
                        if (result[i] == prev) {
                            count++;
                        } else {
                            compressed += std::to_string(count) + prev;
                            prev = result[i];
                            count = 1;
                        }
                    }

                    // Add the last group
                    compressed += std::to_string(count) + prev;
                    *cigar = compressed;
                }
                

                // printAlignmentMatrix(query, query_len, target, target_len, m);
                // if(cigar!= nullptr){ std::cout<<"CIGAR: "<<result<<std::endl<<std::endl; }
                // if(cigar!= nullptr){ std::cout<<"CIGAR: "<<compressed<<std::endl<<std::endl; }
                // std::cout<<"Score: "<<m[gi][gj].cost<<std::endl<<std::endl;
            
                return( m[gi][gj].cost ); /* Steven Skiena, http://www.algorithm.cs.sunysb.edu/computationalbiology/*/
                break;
            case AlignmentType::local:  
                for (i=1; i<=query_len; i++){
                    for (j=1; j<=target_len; j++) {
                        opt[MATCH] = m[i-1][j-1].cost + match_func(query_str[i-1],target_str[j-1],match,mismatch);
                        opt[INSERT] = m[i][j-1].cost + indel(target_str[j-1], gap);
                        opt[DELETE] = m[i-1][j].cost + indel(query_str[i-1], gap);
                        m[i][j].cost = opt[0]; 
                        m[i][j].parent = 0; 
                        for (k=1; k<=2; k++) {
                            if (opt[k] > m[i][j].cost) {
                                m[i][j].cost = opt[k];
                                m[i][j].parent = k;
                            }
                        }
                        if(m[i][j].cost<0){ m[i][j].cost=0;}
                        if (m[i][j].cost>max_cost){
                            max_cost = m[i][j].cost;
                            gi = i;
                            gj = j;
                            local_i = gi;  // start = i for bactracking in determening CIGAR string
                            local_j = gj;  // start = j for bactracking in determening CIGAR string
                        }
                    }   
                } 
                
                // Now (si+1, sj+1) is the starting position of the local alignment
                if (target_begin != nullptr) {
                    *target_begin = local_j + 1;
                }

                if(cigar!=nullptr){
                    while (m[local_i][local_j].cost > 0) {
                        int dir = m[local_i][local_j].parent;
                        if (dir == MATCH) {
                            result.push_back('M');
                            local_i--; local_j--;
                        } else if (dir == INSERT) {
                                result.push_back('I');
                                local_j--;
                        } else if (dir == DELETE) {
                                result.push_back('D');
                                local_i--;
                        } else {
                                throw std::invalid_argument("Unknown error in determining cigar string.");
                                break; 
                        }
                    }
                    std::reverse(result.begin(), result.end());
                    
                    //*cigar = result;

                    char prev = result[0];
                    int count = 1;
                    // Compressing cigar for better output
                    for (size_t i = 1; i < result.size(); ++i) {
                        if (result[i] == prev) {
                            count++;
                        } else {
                            compressed += std::to_string(count) + prev;
                            prev = result[i];
                            count = 1;
                        }
                    }

                    // Add the last group
                    compressed += std::to_string(count) + prev;
                    *cigar = compressed;
                }

                // printAlignmentMatrix(query, query_len, target, target_len, m);
                // if(cigar!= nullptr){ std::cout<<"CIGAR: "<<compressed<<std::endl<<std::endl; }
                // if(cigar!= nullptr){ std::cout<<"CIGAR: "<<result<<std::endl<<std::endl; }
                // std::cout<<"Score: "<<m[gi][gj].cost<<std::endl<<std::endl;

                return( m[gi][gj].cost ); /* Steven Skiena, http://www.algorithm.cs.sunysb.edu/computationalbiology/*/
                            
                break;
            
            case AlignmentType::semiGlobal:
                for (i=1; i<=query_len; i++){
                    for (j=1; j<=target_len; j++) {
                        opt[MATCH] = m[i-1][j-1].cost + match_func(query_str[i-1],target_str[j-1],match,mismatch);
                        opt[INSERT] = m[i][j-1].cost + indel(target_str[j-1], gap);
                        opt[DELETE] = m[i-1][j].cost + indel(query_str[i-1], gap);
                        m[i][j].cost = opt[0];
                        m[i][j].parent = 0;
                        for (k=1; k<=2; k++) {
                            if (opt[k] > m[i][j].cost) {
                                m[i][j].cost = opt[k];
                                m[i][j].parent = k;
                            }
                        }
                    }   
                }
                for (i=0; i<=query_len; i++){
                    if(m[i][target_len].cost>max_cost){
                        max_cost = m[i][target_len].cost;
                        gi = i;
                        gj = target_len;
                    }
                }
                for (j=0; j<=target_len; j++) {
                    if(m[query_len][j].cost>max_cost){
                        max_cost = m[query_len][j].cost;
                        gi = query_len;
                        gj = j;
                    }
                }

                global_i = gi;
                global_j = gj;

                if (target_begin != nullptr) { 
                    *target_begin = 0;
                }
                if (cigar != nullptr) {
                    while (global_i > 0 || global_j > 0) {
                        if (global_i > 0 && global_j > 0 && m[global_i][global_j].parent == MATCH) {
                            result.push_back('M');
                            global_i--;
                            global_j--;
                        } else if (j > 0 && m[global_i][global_j].parent == INSERT) {
                            result.push_back('I');
                            global_j--;
                        } else if (global_i > 0 && m[global_i][global_j].parent == DELETE) {
                            result.push_back('D');
                            global_i--;
                        } else {
                            throw std::invalid_argument("Unknown error in determining cigar string.");
                            break;
                        }
                    }
                    std::reverse(result.begin(), result.end());
                    
                    // Go to the end of column/row
                    if(gj!=target_len || gi!=query_len){ 
                        // Result is in last row
                        if(gi==query_len){
                            result.append(target_len - gj, 'I');
                        } 
                        // Result is in last column
                        else if (gj==target_len) {
                            result.append(query_len - gi, 'D');
                        }
                    }

                    //*cigar = result;

                    char prev = result[0];
                    int count = 1;
                    // Compressing cigar for better output
                    for (size_t i = 1; i < result.size(); ++i) {
                        if (result[i] == prev) {
                            count++;
                        } else {
                            compressed += std::to_string(count) + prev;
                            prev = result[i];
                            count = 1;
                        }
                    }

                    // Add the last group
                    compressed += std::to_string(count) + prev;
                    *cigar = compressed;
                }
                
                // printAlignmentMatrix(query, query_len, target, target_len, m);
                // if(cigar!= nullptr){ std::cout<<"CIGAR: "<<result<<std::endl<<std::endl; }
                // if(cigar!= nullptr){ std::cout<<"CIGAR: "<<compressed<<std::endl<<std::endl; }
                // std::cout<<"Score: "<<m[gi][gj].cost<<std::endl<<std::endl;

                return( m[gi][gj].cost ); /* Steven Skiena, http://www.algorithm.cs.sunysb.edu/computationalbiology/*/

                break;
            
            default:
                throw std::invalid_argument("Unknown AlignmentType provided.");
                return -1; 
        }       
    }
}
