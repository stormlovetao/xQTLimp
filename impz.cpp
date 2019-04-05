#include "impz.h"
#include <iostream>
using namespace std;

void impz(string out_dir , string last_gene_name , vector<int>* p_snps_flag , double** weight ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> impute_flags,vector<long long int>* p_useful_typed_snps) {


	size_t num_total_snps = all_snps.size();
	string OUT_FILE = out_dir + last_gene_name;
	
	
	// load z-score typed snps
	vector<zscore_typed_snp> ztyped_snps;
	size_t num_ztyped_snps = load_zscore_typed_snps(typed_snps, ztyped_snps , p_useful_typed_snps);
	
	// impute for all snps
	double exp_zscore_mean = 0.0;
	FILE* fout = safe_fopen(OUT_FILE.c_str(), "w");
	fprintf(fout, "SNP_name SNP_pos Ref_Allele Alt_Allele Z-Score r2pred\n");
	
//	size_t typed_snp_idx = 0;
	int index = 0;
	for(size_t idx = 0; idx < num_total_snps; idx++) {
		
		// get rid of snp name and snp pos in the first two columns
		// get beta and impute z-score
		if((*p_snps_flag)[idx] == 0)
		{
			continue;
		}
		
		double imp_zscore = 0.0;
		double beta;
		for(size_t i = 0; i < num_ztyped_snps; i++){
			beta = weight[index][i];
			imp_zscore += beta*(ztyped_snps[i].zscore - exp_zscore_mean);
		}
		
		// get variance
		double var;
		var = weight[index][num_ztyped_snps];
		
		// use original z-score for snps that don't need imputation
		if(impute_flags[index] == 0) {
//			imp_zscore = ztyped_snps[typed_snp_idx].zscore;
			var = 1.0;
//			typed_snp_idx++;
		}
		index++;
		// print to file
		fprintf(fout, "%s %lld %c %c %.6lf %.6lf\n",
			all_snps[idx].snp_name.c_str(), all_snps[idx].snp_pos, all_snps[idx].ref_allele,
			all_snps[idx].alt_allele, imp_zscore, var);
	}
	
	fclose(fout);
	for(int i = 0;i < index;i++)
	{
		free(weight[i]);
	}
	free(weight);
	
	

}
