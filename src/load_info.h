#ifndef LOAD_INFO_H
#define LOAD_INFO_H
#include <iostream>
#include <dirent.h>
#include <stdio.h>
#include <string>
#include <string.h>
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

#define GZ_BUF_SIZE 20000


bool gzLoadVcfFile(char* tem , const char* gzfn , map<string,long long int>* p_VcfIndex ,vector<string>* p_VcfFile );
void load_gene_pos_map(string gene_annotation , map<string,long long int*> *gene_pos_map);
void record_all(vector<snps> *p_origin_typed_snps , string  pos , string ref, string alt  ,   string z);
void filter_snps(string file_name , string last_gene_name , vector<snps> *p_origin_typed_snps ,
					vector<typed_snp> *p_typed_snps , 
						vector<snps> *p_ignore_snps,
							 long long int window[2],
						map<string,long long int>* p_VcfIndex,
  						 vector<string>* p_VcfFile);

void read_hap(string line,vector<string>* p_hap);
void gen_map_hap(string ref_file , string last_gene_name ,  vector<ref_snp> *p_snp_map , 
			vector<string>*p_hap , long long int window[2],
			map<string,long long int>* p_VcfIndex,
  						 vector<string>* p_VcfFile);
int travel_eqtl(long long int start , long long int end ,
			string eqtl_path ,vector<long long int>* p_batch_bonder ,int batch);

#endif

