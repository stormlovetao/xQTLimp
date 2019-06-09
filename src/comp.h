#ifndef COMP_H
#define COMP_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <dirent.h>
#include <fstream>
#include <vector>
#include <map> 
#include <stdlib.h>
#include <zlib.h> 
#include <list>
#include <time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include "utils.h" 

typedef struct {
	string name;
	int chrom;
	long long int start;
	long long int end;
	long long int l_win;
	long long int r_win;
}file_record;

using namespace std;
string read_pos(string line);
bool my_cmp(file_record a,file_record b);
void split_line(string line_list[],string line);
void init_line_list(string line_list[]);
void extract(long long int* pos,string* name,string line);
bool read_pos_name_alleles(string line , long long int* pos , string* name,string allele[2]);
void calculate_window(int window_size , long long int window[2],string last_name ,
				map<string,long long int*> *pos_map);
void clean_all_vector(	vector<typed_snp> *p_maf_snps,      
							vector<snps> *p_origin_typed_snps , 
							vector<typed_snp> *p_typed_snps ,
							vector<snps>*p_ignore_snps , 
							vector<ref_snp>*p_snp_map,
							vector<string>* p_hap,
							vector<long long int>* p_useful_typed_snps,
							vector<int>* p_snps_flag);
int get_chrom(string line);
int split_chrom(string Xqtl_path ,long long int chrom[]);
void make_output_dir(long long int chrom[] , int chrom_num  ,char *Out ,int chr);
void organize_files(vector<string>* p_files , int chrom , string out , map<string,long long int*> pos_map);
void print_usage (FILE * stream, int exit_code);

#endif
