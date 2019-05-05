#ifndef COMP_H
#define COMP_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <dirent.h>
#include <fstream>
#include <vector>
#include <map> 
#include <stdlib.h>
#include <math.h>
#include <zlib.h> 
#include <list>
#include <time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include "utils.h" 

using namespace std;
#define window_size 500000 
string read_pos(string line);
void split_line(string line_list[],string line);
void init_line_list(string line_list[]);
void extract(long long int* pos,string* gene_name,string line);
bool read_pos_name_alleles(string line , long long int* pos , string* name,string allele[2]);
void calculate_window(long long int window[2],string last_gene_name ,
				map<string,long long int*> *gene_pos_map);
void clean_all_vector(	vector<typed_snp> *p_maf_snps,      
							vector<snps> *p_origin_typed_snps , 
							vector<typed_snp> *p_typed_snps ,
							vector<snps>*p_ignore_snps , 
							vector<ref_snp>*p_snp_map,
							vector<string>* p_hap,
							vector<long long int>* p_useful_typed_snps,
							vector<int>* p_snps_flag);
int get_chrom(string line);
void split_chrom(string eqtl_path ,long long int chrom[]);
void make_output_dir(char *Out);
void organize_files(string out);

#endif
