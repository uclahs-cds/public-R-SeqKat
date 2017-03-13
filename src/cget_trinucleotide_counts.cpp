#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//[[Rcpp::export]]
Rcpp::NumericVector cget_trinucleotide_counts(std::vector<std::string> key, std::string chr, const int length = 3, long int start = 1, long int end = -1) {

    std::map<std::string, long int> count;

    //initialize
    for (int i = 0; i < key.size(); i++) 
        count.insert( std::pair<std::string, long int>(key[i], 0) );

    ifstream f(chr.c_str());

    char tri_nucl_c[length];
    std::string tri_nucl;

    if(f.is_open()) {
        
        //find the length of genome, chr_length=length_of_chr + 1(end_of_file);
        f.seekg(0, f.end);
        long int chr_length = f.tellg();

        f.seekg(start - 1, f.beg);

        long int back_move_step = length - (length / 2);

        if (end < 0 || end > chr_length) end = chr_length;
        if (start < 1) start = 1;

        for (long int i = start - 1; i <= end - length; i++) {
            f.get(tri_nucl_c, 4); //extract the next 3 charactors, file pointer moves 3 chars
            tri_nucl = tri_nucl_c;
            transform(tri_nucl.begin(), tri_nucl.end(), tri_nucl.begin(), ::toupper);
            f.seekg(- back_move_step, ios_base::cur); //move pointer back to 'next of beginning of this loop'
            //printf("%s\n", tri_nucl.c_str());
            ++count.at(tri_nucl);
        }

    } else {
        printf("Fail to open file!");
        exit(EXIT_FAILURE);
    }

    f.close();

    Rcpp::NumericVector count_vec(key.size(),0);

    for (int i = 0; i < key.size(); ++i)
        count_vec[i] = count.at(key[i]);

    return count_vec;
}