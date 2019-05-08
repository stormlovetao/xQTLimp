/*
###############################################################################################
## This file provide basic functions for imputation process ,including matrix operations and operations on LD
## This file is adapted from work of Bogdan Pasaniuc et. al [1]. 
## Reference: [1]. Pasaniuc, Bogdan, et al. "Fast and accurate imputation of summary statistics enhances evidence of functional enrichment." Bioinformatics 30.20 (2014): 2906-2914.
################################################################################################
*/

#ifndef  UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string>
#include <string.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <map>

using namespace std;

#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define ROTATE(a,i,j,k,l) g=a[i*n+j];h=a[k*(n)+l];a[i*(n)+j]=g-s*(h+g*tau);\
	a[k*(n)+l]=h+s*(g-h*tau);

const int verbose = 1;

typedef long long int pos_t;

// provided by user
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	string ref_allele;
	string alt_allele;
	double zscore;
} typed_snp;

typedef struct{
	string snp_name;
	pos_t snp_pos;
	string ref_allele;
	string alt_allele;
	double zscore;
}snps; 

// from internal file
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	string ref_allele;
	string alt_allele;
} ref_snp;

// from internal file
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	double zscore;
} zscore_typed_snp;


// to store a genotype
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	string zero;
	string one;
	string two;
} genotype;


FILE* safe_fopen(const char* filename, const char *op);

void* safe_calloc(size_t nelements, size_t sizeof_element);

// search for the index of the snp whose pos is specified by pos
// assume the snps in typed_snps are sorted by snp positions
bool search_by_pos(const vector<typed_snp> &typed_snps,
			pos_t pos, size_t &idx);

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<ref_snp> &all_snps, pos_t pos, size_t &idx);

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<genotype> &all_genos, pos_t pos, size_t &idx);

// mark snps that are typed in snp_flags to 1
// mark snps that need to convert z-scores in convert_flags to 1
size_t mark_snps(vector<typed_snp> &typed_snps, 
			const vector<ref_snp> &snp_ref_alt,vector<char>& convert_flags,
			vector<char>& impute_flags);

// load haplotypes from the 1000 genome project reference panel file
size_t load_haplotypes(const char *filename,
					vector<string> &haps, size_t num_total_snps);



/* linear algebra */
void pdinv(double *cinv, double *coeff, int n) ;

// compute pseudo inverse
void pinv(double *A_inv, double *A, size_t n);

// compute pseudo inverse using jacobi
void pinv_jacobi(double *A_inv, double *A, size_t n);

/* numer recipes p 97 */
int choldc (double *a, int n, double p[]);
void cholsl (double *a, int n, double p[], double b[], double x[]);
void cholesky(double *cf, double *a, int n) ;
void pmat(double *mat, int n)   ;


// load z-score typed snps file
size_t load_zscore_typed_snps(const vector<typed_snp>& typed_snps,
	 vector<zscore_typed_snp> &ztyped_snps,
	vector<long long int>* p_useful_typed_snps);

// compute allele frequencies for all the snps
void get_all_freqs(const vector<string>& haps, vector<double> &freqs);

// compute the allele frequency of a snp, specified by snp_idx
double get_freq(const vector<string>& haps, size_t snp_idx);

// compute the h frequency between two snps, specified by snp_idx1, snp_idx2
double get_h_freq(const vector<string>& haps,size_t snp_idx1,size_t snp_idx2);

// get variance
double get_var(const vector<string>& haps,
		const vector<double>& freqs, int idx);


#endif



