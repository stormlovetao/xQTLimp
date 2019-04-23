#ifndef IMPZ_H
#define IMPZ_H
#include <cmath>
#include <stdio.h>
#include <string>
#include "impg.h"
void impz(int chrom , string out_dir , string last_gene_name , vector<int>* p_snps_flag , double** weight ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> impute_flags,vector<long long int>* p_useful_typed_snps);

#endif
