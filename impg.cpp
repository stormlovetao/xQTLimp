#include <string.h>
#include <cmath>

#include "impg.h"

using namespace std;

// load z-score typed snps file
size_t load_zscore_typed_snps(const vector<typed_snp>& typed_snps,
 vector<zscore_typed_snp> &ztyped_snps
,vector<long long int>* p_useful_typed_snps) {
	size_t idx = 0;
	int len = (*p_useful_typed_snps).size();
	int typed_len = typed_snps.size();
	for(int i = 0;i < len;i++)
	{
		zscore_typed_snp tmp;
		tmp.idx = idx;
		for(int j = 0;j < typed_len;j++)
		{
			if((*p_useful_typed_snps)[i] == typed_snps[j].snp_pos)
			{
				tmp.snp_name = typed_snps[j].snp_name;
				tmp.snp_pos = typed_snps[j].snp_pos;
				tmp.zscore = typed_snps[j].zscore;
				ztyped_snps.push_back(tmp);
				idx++;
				break; 
			}
		}
		
	}
	
	size_t num_ztyped_snps = ztyped_snps.size();
	
	return num_ztyped_snps;
}

// load genotypes
size_t load_genotypes(const char *filename, 
	vector<genotype> &genotypes) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	// read in the rest
	size_t snp_idx = 0;
	char snp_name[1024];
	char zero[8];
	char one[8];
	char two[8];
	pos_t snp_pos;
	while(fscanf(fin, "%s %lld %s %s %s", snp_name, &snp_pos,
			zero, one, two)>0) {
		genotype tmp;
		tmp.idx = snp_idx;
		tmp.snp_name = string(snp_name);
		tmp.snp_pos = snp_pos;
		tmp.zero = string(zero);
		tmp.one = string(one);
		tmp.two = string(two);
		genotypes.push_back(tmp);
		snp_idx++;
	}
	fclose(fin);
	
	size_t num_genotypes = genotypes.size();
	if(verbose) {
		printf("Info: Loaded %u genotypes from %s...\n",
			(unsigned int)num_genotypes, filename);
	}
	
	return num_genotypes;
}

// load the ld matrix from file to memory
void load_ld_mat(const char *filename, size_t num_typed_snps, double *ld_mat) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = i+1; j < num_typed_snps; j++) {
			double corr;
			if(fscanf(fin, "%lf", &corr)) {
				ld_mat[i*num_typed_snps+j] = corr;
				ld_mat[j*num_typed_snps+i] = corr;
			}
			else {
				fprintf(stderr, "Error: Wrong LD matrix...\n");
				exit(1);
			}
		}
		ld_mat[i*num_typed_snps+i] = 1.0;
	}
	
	if(verbose) {
		printf("Info: Loaded LD matrix from %s...\n", filename);
	}
}

// compute allele frequencies for all the snps
void get_all_freqs(const vector<string>& haps, vector<double> &freqs) {
	size_t num_total_snps = freqs.size();
	for(size_t i = 0; i < num_total_snps; i++) {
		freqs[i] = get_freq(haps, i);
	}
}

// compute the allele frequency of a snp, specified by snp_idx
double get_freq(const vector<string>& haps, size_t snp_idx) {
	double freq = 0.0;
	size_t nhaps = haps[0].size();
	size_t i;
	for(i = 0; i < nhaps; i++) {
		if(haps[snp_idx][i] == '1') {
			freq += 1.0;
		}
	}
	freq = freq/((double)nhaps);
	return freq;
}

// compute the h frequency of two snps, specified by snp_idx1, snp_idx2
double get_h_freq(const vector<string>& haps, 
			size_t snp_idx1, size_t snp_idx2) {
	double freq = 0.0;
	size_t nhaps = haps[0].size();
	size_t i;
	for(i = 0 ; i < nhaps; i++) {
		if(haps[snp_idx1][i] == '1' && haps[snp_idx2][i] == '1') {
			freq += 1.0;
		}
	}
	freq = freq/((double)nhaps);
	return freq;
}

double get_var(const vector<string>& haps,
		const vector<double>& freqs, int idx) {
	double mn = freqs[idx];
	double vr = 0.0;
	size_t nhaps = haps.size();
	for(size_t i = 0; i < nhaps; i++) {
		if(haps[i][idx] == '1') {
			vr += (1.0-mn)*(1.0-mn);
		}
		else {
			vr += (0.0-mn)*(0.0-mn);
		}
	}
	return (vr/((double)nhaps-1.0));
}


// get the command line input for gen_beta program
int get_gen_beta_cmd_line(int argc, char **argv, char **IN_HAP_FILE,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE,
		char **OUT_FILE_PREFIX, char **MAF_TH, char **LAMBDA) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"h:m:t:p:f:l:")) != -1) {
		switch(opt) {
			case 'h': // haplotype file
				*IN_HAP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'm': // all snps file
				*IN_ALL_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 't': // typed snps file
				*IN_TYPED_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'p': // prefix for output files
				*OUT_FILE_PREFIX = (char *) strdup(optarg);
				nflags++;
				break;
            case 'f': // MAF threshold
                *MAF_TH = (char*) strdup(optarg);
                nflags++;
                break;
            case 'l': // LAMBDA
                *LAMBDA = (char*) strdup(optarg);
                nflags++;
                break;
		}
	}
	if(!(*IN_HAP_FILE) || !(*IN_ALL_SNP_FILE) || !(*IN_TYPED_SNP_FILE) ||
       !(*OUT_FILE_PREFIX)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-h (required) specify haplotype file\n");
		fprintf(stderr, "\t-m (required) specify SNP mapping file\n");
		fprintf(stderr, "\t-t (required) specify typed SNP file\n");
		fprintf(stderr, "\t-p (required) specify output file prefix\n");
        fprintf(stderr, "\t-f (optional) specify minimum MAF (0.01 by default)\n");
		fprintf(stderr, "\t-l (optional) specify lambda (0.1 by default)\n");
        exit(1);
	}
	
	return nflags;
}

// get command line input for impute z-scores
int get_imp_cmd_line(int argc, char **argv, char **PREFIX,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE, char **OUT_FILE) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"p:m:t:o:")) != -1) {
		switch(opt) {
			case 'p': // prefix used in gen_beta
				*PREFIX = (char *) strdup(optarg);
				nflags++;
				break;
			case 'm': // snp mapping file
				*IN_ALL_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 't': // typed snp file
				*IN_TYPED_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'o': // output file name
				*OUT_FILE = (char *) strdup(optarg);
				nflags++;
				break;
		}
	}
	if(!(*PREFIX) || !(*IN_ALL_SNP_FILE) || !(*IN_TYPED_SNP_FILE) ||
       !(*OUT_FILE)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-p (required) specify input file prefix\n");
		fprintf(stderr, "\t-m (required) specify SNP mapping file\n");
		fprintf(stderr, "\t-t (required) specify typed SNP file\n");
		fprintf(stderr, "\t-o (required) specify output file name\n");
		exit(1);
	}
	
	return nflags;
}

// get the command line input for gen_beta with ld program
int get_gen_beta_ld_cmd_line(int argc, char **argv, char **IN_HAP_FILE,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE, char **IN_LD_FILE,
		char **OUT_FILE_PREFIX, char **MAF_TH) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"h:m:t:p:l:f:")) != -1) {
		switch(opt) {
			case 'h': // haplotype file
				*IN_HAP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'm': // snp mapping file
				*IN_ALL_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 't': // typed snps file
				*IN_TYPED_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'l': // name for ld
				*IN_LD_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'p': // prefix for output files
				*OUT_FILE_PREFIX = (char *) strdup(optarg);
				nflags++;
				break;
            case 'f': // MAF threshold
                *MAF_TH = (char*) strdup(optarg);
                nflags++;
                break;

		}
	}
	if(!(*IN_HAP_FILE) || !(*IN_ALL_SNP_FILE) || !(*IN_TYPED_SNP_FILE) ||
       !(*IN_LD_FILE) || !(*OUT_FILE_PREFIX)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-h (required) specify haplotype file\n");
		fprintf(stderr, "\t-m (required) specify SNP mapping file\n");
		fprintf(stderr, "\t-t (required) specify typed SNP file\n");
		fprintf(stderr, "\t-l (required) specify LD statistics file\n");
		fprintf(stderr, "\t-p (required) specify output file prefix -p\n");
		fprintf(stderr, "\t-f (optional) specify minimum MAF (0.01 by default)\n");
        exit(1);
	}
	
	return nflags;
}
