#ifndef IMPUTE_H
#define IMPUTE_H
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string>


#include "utils.h"

typedef struct answer
{
	int row;
	double** weight;
	double** last_sigma_it;
}ans;

using namespace std;

ans zgenbt(vector<typed_snp> *p_maf_snps ,double maf , double lam ,vector<int>* p_snps_flag ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> convert_flags,vector<char> impute_flags ,
vector<string> haps,vector<long long int>* p_useful_typed_snps,
map<long long int,int> *p_m_all_snps,map<long long int , int> *p_m_typed_snps,
double** last_sigma_it , int row);


void impz(vector<typed_snp> *p_maf_snps , vector<snps> *p_ignore_snps , int chrom , string out_dir , string last_name , vector<int>* p_snps_flag , double** weight ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> impute_flags,vector<long long int>* p_useful_typed_snps);

#endif
