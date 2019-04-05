#include "zgenbt.h"

ans zgenbt(vector<int>* p_snps_flag ,vector<typed_snp> typed_snps,vector<ref_snp> all_snps,
vector<char> convert_flags,vector<char> impute_flags , 
	vector<string> haps,vector<long long int>* p_useful_typed_snps)
{
	size_t num_typed_snps = typed_snps.size();
	size_t num_total_snps = all_snps.size();
	double maf_th = 0.01;
	double lambda = 0.1;
	


	
	// compute allele frequencies for all snps
	vector<double> freqs;
	freqs.resize(num_total_snps, 0.0);
	get_all_freqs(haps, freqs);
	
	int counter = 0;
	for(int i = 0;i < num_total_snps;i++)
  	{
		if(freqs[i] < maf_th || freqs[i] > 1 - maf_th)
		{
			(*p_snps_flag).push_back(0);
			counter++;
		}
		else
		{
			(*p_snps_flag).push_back(1);
		}
	}
	
	

	// t,t estimate sigma_t matrix using maf filtered typed snps
	size_t sigma_t_tmp_size = num_typed_snps*num_typed_snps;
	double *sigma_t_tmp = (double*)safe_calloc(sigma_t_tmp_size,
		sizeof(double));
	
	for(size_t i = 0; i < num_typed_snps; i++) {
		size_t idxi = typed_snps[i].idx;
		double pi = freqs[idxi];
		for(size_t j = i+1; j < num_typed_snps; j++) {
			size_t idxj = typed_snps[j].idx;
			double pj = freqs[idxj];
			double pij = get_h_freq(haps, idxi, idxj);
			double r = (pij-pi*pj)/sqrt(pi*(1.0-pi)*pj*(1.0-pj));  
			sigma_t_tmp[i*num_typed_snps+j] = r;
			sigma_t_tmp[j*num_typed_snps+i] = sigma_t_tmp[i*num_typed_snps+j];
		}
		sigma_t_tmp[i*num_typed_snps+i] = 1.0;
	}
	// output the betas and vars

		
	// for filtering
	vector<char> filter_flags;
	filter_flags.resize(num_typed_snps, 0);

	// filter perfectly correlated snps
	/*
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = i+1; j < num_typed_snps; j++) {
			if(sigma_t_tmp[i*num_typed_snps+j] > 0.99) {
				filter_flags[j] = 1;
			}
		}
	}
	*/
	// filter maf snps
	for(size_t i = 0; i < num_typed_snps; i++) {
		double freq = freqs[typed_snps[i].idx];
		if(freq < maf_th || freq > 1-maf_th) {
			filter_flags[i] = 1;
			
		}
		else{
			(*p_useful_typed_snps).push_back(typed_snps[i].snp_pos);
		}
	}
	
	// count the number of unfiltered snps
	size_t num_snps = 0;
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			num_snps++;
		}
	}
	
	// for debugging
	
	double** weight = (double**)malloc((num_total_snps - counter) * sizeof(double*));
	for(int i = 0 ; i < (num_total_snps - counter) ; i++)
	{
		weight[i] = (double*)malloc((num_snps + 1)* sizeof(double));
	}
	size_t sigma_t_size = num_snps * num_snps;
	double *sigma_t = (double*)safe_calloc(sigma_t_size, sizeof(double));
	double *sigma_t_inv = (double*)safe_calloc(sigma_t_size, sizeof(double));
	double *sigma_it = (double*)safe_calloc(num_snps, sizeof(double));
	double *beta = (double*)safe_calloc(num_snps, sizeof(double));

	// construct sigma for typed snps
	size_t ii = 0;
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = 0; j < num_typed_snps; j++) {
			if(filter_flags[i] != 1 && filter_flags[j] != 1) {
				sigma_t[ii] = sigma_t_tmp[i*num_typed_snps+j];
				if(i == j) {
					sigma_t[ii] += lambda;
				}
				ii++;
			}
		}
	}

	// compute the inverse
	pdinv(sigma_t_inv, sigma_t, num_snps);
		
	//compute the weight
	int mem = 0;
	int index1 = 0; 
	for(int idx = 0; idx < num_total_snps; idx++) {
		// for SNPs that need imputation
		if((*p_snps_flag)[idx] == 0)
		{
			continue;
		}
		
		if(impute_flags[idx] == 1) {
			// get sigma_it
			ii = 0;
			double pi = freqs[idx];
			for(size_t i = 0; i < num_typed_snps; i++) {
				if(filter_flags[i] != 1) {
						size_t sidx = typed_snps[i].idx;
						double pj = freqs[sidx];
						double r = 0.0;
						if(pi < maf_th || pi > 1-maf_th) 
						{
						// do nothing
						}
						else 
						{
							double pij = get_h_freq(haps, idx, sidx);
							r = (pij-pi*pj)/sqrt(pi*(1.0-pi)*pj*(1.0-pj));
						}
						sigma_it[ii] = r;
						ii++;
				}
			}
			//got sigma_it 
			//
			for(size_t i = 0; i < num_snps; i++) {
				beta[i] = 0.0;
					for(size_t j = 0; j < num_snps; j++) 
					{					
					// valid because sigma_t_inv is symmetric
						beta[i] += sigma_it[j]*sigma_t_inv[num_snps*i+j];
					}
				weight[index1][i] = beta[i];
			}
		}
		// for SNPs that don't need imputation
		else {
			int qw = 0;
			for(size_t i = 0; i < num_typed_snps; i++) {
				if(filter_flags[i] == 1)
				{
					continue;
				}
				if(typed_snps[i].idx == idx) {
					weight[index1][qw] = 1;
					qw++;
				}
				else {
					weight[index1][qw] = 0;
					qw++;
				}
			}
		}
		
// ----------------------------------------------------------------		
		// print variance
		// for SNPs that need imputation
		if(impute_flags[idx] == 1) {
			double var = 0.0;
			for(size_t i = 0; i < num_snps; i++) {
				for(size_t j = 0; j < num_snps; j++) {
					// valid because sigma_t is symmetric			
					var += beta[i]*beta[j]*sigma_t[num_snps*i+j];
				}
			}
			weight[index1][num_snps] = var;
		}
		// for SNPs that don't need imputation
		else {
			weight[index1][num_snps] = 1.0;
		}
			
		index1++;
	} //for
	
		
	
	free(sigma_t_tmp);	
	free(sigma_t);
	free(sigma_t_inv);
	free(sigma_it);
	free(beta);
		
	ans values;
	values.weight = weight;
	
	return values;	

}

