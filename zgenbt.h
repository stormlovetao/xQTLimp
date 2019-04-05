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
}ans;

using namespace std;

ans zgenbt( vector<int>* p_snps_flag ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> convert_flags,vector<char> impute_flags ,
vector<string> haps,vector<long long int>* p_useful_typed_snps);
#endif
