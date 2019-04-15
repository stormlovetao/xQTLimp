#ifndef ZGENBT_H
#define ZGENBT_H
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>


#include "impg.h"
#include "linsubs.h"


typedef struct answer
{
	double** weight;
	double** last_sigma_it;
}ans;

using namespace std;

ans zgenbt(vector<int>* p_snps_flag ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> convert_flags,vector<char> impute_flags ,
vector<string> haps,vector<long long int>* p_useful_typed_snps,
map<long long int,int> *p_m_all_snps,map<long long int , int> *p_m_typed_snps,
double** last_sigma_it);
#endif
