#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <map> 
#include <malloc.h>
#include <stdlib.h>
#include <pthread.h> 
#include <semaphore.h>
#include <math.h>
#include <zlib.h> 
#include <list>
#include <time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include "zgenbt.h"
#include "impz.h"
#include "util.h" 

#define window_size 500000 

#define GZ_BUF_SIZE 20000

using namespace std;
int batch =  35;
int real_batch = 0;
double rate = 0.1;

struct thread_data{
   map<string,long long int>* p_VcfIndex;
   vector<string>* p_VcfFile;
   map<string,long long int*> *gene_pos_map;
   int chrom; 
   string eqtl_path; 
   string ref_file;
   string out_dir;
   long long int start; 
   long long int end;
   sem_t * bin_sem;
};


string read_pos(string line)
{
	string pos = "";
	int i = 0;
	while(line[i] != '\t')
	{
		i++;
	}
	while(line[i] == '\t')
	{
		i++;
	}
	while(line[i] != '\t')
	{
		pos += line[i];
		i++;
	}
	return pos;
	
} 

bool gzLoadVcfFile(const char* gzfn , map<string,long long int>* p_VcfIndex ,vector<string>* p_VcfFile )
{
	//open .gz file
	gzFile gzfp = gzopen(gzfn,"rb");
	if(!gzfp)
	{
		cout << "can't open\n";
		return false;
	}
	char buf[GZ_BUF_SIZE];

	//erase head
	gzgets(gzfp , buf , GZ_BUF_SIZE);
    while(buf[0] == '#')
    {
        gzgets(gzfp , buf , GZ_BUF_SIZE);
    }
        
    string line = string(buf);
    string pos = read_pos(line);
    (*p_VcfIndex)[pos] = 0;
    (*p_VcfFile).push_back(line);
    long long int index = 1;


	while( gzgets(gzfp , buf , GZ_BUF_SIZE))
	{
		string line = string(buf);
		string pos = read_pos(line);
		if((*p_VcfIndex).count(pos) == 1)
		{
		//	cout << "conflict information in vcf file!\n";
		}
		else
		{
			(*p_VcfIndex)[pos] = index;
			index++;
			(*p_VcfFile).push_back(line);	
		}
	}
	
	//close .gz file
	gzclose(gzfp);
	return true;
}

void split_line(string line_list[],string line)
{
	int flag = 0;
	int index = 0;
	
	for(int i = 0;i < line.length();i++)
	{
		flag = 0;
		if(line[i] == '\t')
		{
			flag = 1;
			i++;
		}
		
		if(flag == 1)
		{ 
			index++;
		}
		line_list[index] += line[i];
	}
	
}

void init_line_list(string line_list[])
{
	for(int i = 0;i < 9;i++)
	{
		line_list[i] = "";
	}
}

//extract info from given line
void extract(long long int* pos,string* gene_name,string line)
{
	//read start pos
	int index = 0;
	while(line[index] != '\t')
	{
		index++;
	}
	while(line[index] == '\t')
	{
		index++;
	}
	string pos1 = "";
	while(line[index] != '\t')
	{
		pos1 += line[index];
		index++;
	}
	while(line[index] == '\t')
	{
		index++;
	} 
	pos[0] = atoi(pos1.c_str());
	//read end pos
	string pos2 = "";
	while(line[index] != '\t')
	{
		pos2 += line[index];
		index++; 
	}
	while(line[index] == '\t' || line[index] == '+' || line[index] == '-')
	{
		index++;
	}
	pos[1] = atoi(pos2.c_str());
	//read gene name
	while(line[index] != '\t')
	{
		(*gene_name) += line[index];
		index++;
	}
}

//record postions of all gene
void load_gene_pos_map(string gene_annotation , map<string,long long int*> *gene_pos_map)
{
	ifstream fin(gene_annotation.c_str());
	string line;
	getline(fin,line);
	while(line != "")
	{
		//pos0:start pos  pos1:end pos
		long long int* pos = (long long int*)malloc(2 * sizeof(long long int));
		string gene_name = "";
		string* cur_gene_name = &gene_name;
		//extract info from the line
		extract(pos , cur_gene_name ,line);
		//record it into map
		(*gene_pos_map)[gene_name] = pos;
		getline(fin,line);
	}
	fin.close();
}


void record_all(vector<snps> *p_origin_typed_snps , string info , string slope, string slope_se)
{
	string line_list[4];
	int index = 0;
	int i = 1; 
	while(info[index] != '_')
	{
		line_list[0] += info[index];
		index++;
	}

	while(i < 4)
	{
		index++;
		while(info[index] != '_')
		{
			line_list[i] += info[index];
			index++;
		}
		i++;
		
	}

	snps temp;
	temp.snp_pos = atoi(line_list[1].c_str());
	temp.ref_allele = line_list[2];
	temp.alt_allele = line_list[3];
	temp.zscore = double(atof(slope.c_str())) / double(atof(slope_se.c_str()));
	(*p_origin_typed_snps).push_back(temp);	
	
	
}

bool read_pos_name_alleles(string line , long long int* pos , string* name,string allele[2])
{
	int index = 0;
	while(line[index] != '\t')
	{
		index++;
	}
	while(line[index] == '\t')
	{
		index++;
	}
	
	
	string str_pos = "";
	while(line[index] != '\t')
	{
		str_pos += line[index];
		index++;
	}
	while(line[index] == '\t')
	{
		index++;
	}
	(*pos) = atoi(str_pos.c_str());
	
	(*name) = "";
	while(line[index] != '\t')
	{
		(*name) += line[index];
		index++;
	}
	while(line[index] == '\t')
	{
		index++;
	}
	
	string ref = "";	
	while(line[index] != '\t')
	{
		ref += line[index];
		index++;
	}
	while(line[index] == '\t')
	{
		index++;
	}
	
	string alt = "";
	while(line[index] != '\t')
	{
		alt += line[index];
		index++;
	}
	
	while(line[index] != '\t')
	{
		alt += line[index];
		index;
	}
	allele[0] = ref;
	allele[1] = alt;
	
	return true;
}

void filter_snps(string file_name , string last_gene_name , vector<snps> *p_origin_typed_snps ,
					vector<typed_snp> *p_typed_snps , 
						vector<snps> *p_ignore_snps,
							 long long int window[2],
						map<string,long long int>* p_VcfIndex,
  						 vector<string>* p_VcfFile)
{
	ifstream fin(file_name.c_str());
	int origin_snps_size = (*p_origin_typed_snps).size();
	int counter = 0;
	int counter1 = 0;
	long long int left = (window[0] < window[1]) ? window[0] : window[1];
	long long int right = (window[0] > window[1]) ? window[0] : window[1];
	for(int i = 0;i < origin_snps_size;i++)
	{
		char num[20]; 
		long long int the_pos = (*p_origin_typed_snps)[i].snp_pos;
		if(the_pos < left || the_pos > right)
		{
			counter++;
			(*p_ignore_snps).push_back((*p_origin_typed_snps)[i]);
			continue;
		}  

		sprintf(num , "%lld" ,the_pos);	
		if((*p_VcfIndex).count(string(num)) == 0)//map中不存在信息 
		{
			(*p_ignore_snps).push_back((*p_origin_typed_snps)[i]);
			//record it into ignore snps 
		}
		else // map中存在信息 
		{
				char num[20];
				sprintf(num , "%lld" ,(*p_origin_typed_snps)[i].snp_pos );
				long long int v_index = (*p_VcfIndex)[num];
				string line = (*p_VcfFile)[v_index];
				
				string allele[2];
				string name;
				string *p_name = &name; 
				long long int pos;
				long long int *p_pos = &pos;
				bool ans = read_pos_name_alleles(line , p_pos , p_name ,allele);				
				if(ans && ((allele[0] == (*p_origin_typed_snps)[i].ref_allele
				&& allele[1] == (*p_origin_typed_snps)[i].alt_allele)
			||(allele[0] == (*p_origin_typed_snps)[i].alt_allele
				&& allele[1] == (*p_origin_typed_snps)[i].ref_allele )) )
				
				{
					if(isnan((*p_origin_typed_snps)[i].zscore) == true)
					{
						(*p_ignore_snps).push_back((*p_origin_typed_snps)[i]);
					}
					else
					{
						typed_snp tem;
						tem.snp_pos = (*p_origin_typed_snps)[i].snp_pos;
						tem.snp_name = name;
						tem.ref_allele = allele[0];
						tem.alt_allele = allele[1];
						tem.zscore = (*p_origin_typed_snps)[i].zscore;
						
						(*p_typed_snps).push_back(tem);
					//record it into typed snps
						
					}
		
				}
				else
				{
					(*p_ignore_snps).push_back((*p_origin_typed_snps)[i]);
					//record it into ignore snps
				}
			
		}

	}//for
	
	fin.close();
	return;					
}

void calculate_window(long long int window[2],string last_gene_name ,
				map<string,long long int*> *gene_pos_map)
{

	long long int l_pos = (*gene_pos_map)[last_gene_name][0];
	long long int r_pos = (*gene_pos_map)[last_gene_name][1];
	long long int pos = (l_pos + r_pos) / 2;
	long long int pos1 = pos - (window_size) / 2;
	if(pos1 < 0)
	{
		pos1 = 0;
	}
	window[0] = pos1;
	window[1] = pos + (window_size) / 2; 
					
}

void read_hap(string line,vector<string>* p_hap)
{		
	int index = 0;
	int len = line.length();
	for(int i = 0;i < 9;i++)
	{
		while(line[index] != '\t')
		{
			index++;
		}
		while(line[index] == '\t')
		{
			index++;
		}
	}
	
	string hap_format = "";
	while(index < len)
	{
		if(line[index] == '0')
		{
			hap_format += '1'; 
		}
		if(line[index] == '1')
		{
			hap_format += '2';
		}
		
		index += 2;
		if(index >= len)
		{
			break;
		}
		if(line[index] == '0')
		{
			hap_format += '1';
		}
		if(line[index] == '1')
		{
			hap_format += '2';
		}
		index += 2;		
	}
	(*p_hap).push_back(hap_format);
}

void gen_map_hap(string ref_file , string last_gene_name ,  vector<ref_snp> *p_snp_map , 
			vector<string>*p_hap , long long int window[2],
			map<string,long long int>* p_VcfIndex,
  						 vector<string>* p_VcfFile)
{

	ifstream fin(ref_file.c_str());
	size_t index = 0;
	
	long long int left = (window[0] < window[1]) ? window[0] : window[1];
	long long int right = (window[0] > window[1]) ? window[0] : window[1];
	
	char num[50];
	long long int pos = left;
	sprintf(num , "%lld" , pos);
	while(((*p_VcfIndex).count(string(num)) == 0 ) && (pos <= right))
	{
		pos++;
	    sprintf(num , "%lld" , pos);	
	}
	
	if(pos <= right)
	{
		long long int v_index = (*p_VcfIndex)[string(num)];
		string line = (*p_VcfFile)[v_index];
		long long int pos1;
		long long int *p_pos = &pos1;
		string name1;
		string *name = &name1;
		string allele[2];	
		bool right_snp = read_pos_name_alleles(line , p_pos , name, allele);
		while((pos1 >= left) && (pos1 <= right))
		{
			if(right_snp)
			{
				ref_snp temp;
				temp.idx = index;
				temp.snp_name = name1;
				temp.snp_pos = pos1;
				temp.ref_allele = allele[0];
				temp.alt_allele = allele[1];
				(*p_snp_map).push_back(temp);
				index++;
				read_hap(line , p_hap);
						
			}
			v_index++;
			if(v_index < (*p_VcfIndex).size())
			{
				line = (*p_VcfFile)[v_index];
			}
			else
			{

				break;
			}
			if(line == "")
			{

				break;
			}
			right_snp = read_pos_name_alleles( line , p_pos , name , allele);
		}
		
	}
	fin.close();
}

void clean_all_vector(      vector<snps> *p_origin_typed_snps , 
							vector<typed_snp> *p_typed_snps ,
							vector<snps>*p_ignore_snps , 
							vector<ref_snp>*p_snp_map,
							vector<string>* p_hap,
							vector<long long int>* p_useful_typed_snps)
{
	(*p_origin_typed_snps).clear();
	(*p_typed_snps).clear();
	(*p_ignore_snps).clear();
	(*p_snp_map).clear();
	(*p_hap).clear();
	(*p_useful_typed_snps).clear();
								
}

void travel_eqtl(long long int start , long long int end ,
			string eqtl_path ,vector<long long int>* p_batch_bonder)
{
	cout << "Dividing chrom into batches....\n";
	ifstream fin(eqtl_path.c_str());
	string line;
	fin.seekg(start);
	int counter = 1;
	getline(fin,line);
	string line_list[9];
	split_line(line_list,line);
	string gene_name = line_list[0];
	string last_gene_name;
	while(fin.tellg() != end)
	{
		last_gene_name = gene_name;
		
		getline(fin,line);
		if(line == "")
		{
			break;
		}
		init_line_list(line_list);
		split_line(line_list,line);
		gene_name = line_list[0];
		
		if(gene_name == last_gene_name)
		{
			
		}
		else
		{
			counter++;
		}
	}

//	counter = 2397;
	fin.close();
	cout << "Tere are " << counter  << " gene on the chrom\n";
	///////////////////////////////////////////////////
	
	ifstream fin1(eqtl_path.c_str());
	int batch_size = counter / batch;
	if(counter % batch == 0)
	{
		real_batch = batch;
	}
	else
	{
		real_batch = batch + 1;
	}
	
	(*p_batch_bonder).push_back(start);
		
	fin1.seekg(start);
	int gene_num = 1;
	getline(fin1,line);
	string line_list1[9];
	split_line(line_list1,line);
	gene_name = line_list1[0];
	long long int pos;
	int b = 0;
	long long int position = 0;
	while(true)
	{
		last_gene_name = gene_name;
		if(fin1.tellg() == end)
		{
			pos = fin1.tellg();
			(*p_batch_bonder).push_back(pos);
			//cout << "last " << pos << endl;
			break;
		}
		
		position = fin1.tellg();
		getline(fin1 , line);
		init_line_list(line_list1);
		split_line(line_list1,line);
		gene_name = line_list1[0];
		
		if(gene_name == last_gene_name)
		{
			
		}
		else
		{
			gene_num++;
			if((gene_num - 1) % batch_size == 0)
			{
			    b++;
				(*p_batch_bonder).push_back(position);
				(*p_batch_bonder).push_back(position);
				//cout << "batch" << b << " end at "<< position<<"(file pos)" << endl;
				
				
			}
		}
	}
	
	fin1.close();
	
}

int get_chrom(string line)
{
	int index = 0;
	while(line[index] != '\t')
	{
		index++;
	}
	index++;
	string chrom = "";
	while(line[index] != '_')
	{
		chrom += line[index];
		index++;
	}
	int ans = atoi(chrom.c_str());
	return ans;
	
}

void split_chrom(string eqtl_path ,long long int chrom[])
{
	for(int i = 0;i <= 22;i++)
	{
		chrom[i] = 0;
	}
	cout << "Dividing eqtl file into several chroms.....\n" ;
	ifstream fin(eqtl_path.c_str());
	string line;
	getline(fin , line);
	//erase head
	chrom[0] = fin.tellg();
/*	chrom[1] = 1359280869;
	chrom[2] = 2302908974;
*/	int current_chrom = 1;
	long long int pos = 0;
	while(current_chrom <= 22)
	{
		pos = fin.tellg();
		getline(fin , line);
		if(line == "")
		{
			chrom[current_chrom] = pos;
			//cout << "chrom" << current_chrom << " end at " << last_pos <<"(file pos)" << endl;
			break;
		}
		int now_chrom = get_chrom(line);
		if(current_chrom == now_chrom)
		{
			//do nothing
		}
		else
		{
			chrom[current_chrom] = pos;
	//		cout << "chrom" << current_chrom << " end at " << pos <<"(file pos)" << endl;
			current_chrom++;	
		} 
	}
	
}

void *main_process(void *threadarg)
{
		struct thread_data *my_data;
		my_data = (struct thread_data *) threadarg;
   		map<string,long long int>* p_VcfIndex = my_data -> p_VcfIndex;
   		vector<string>* p_VcfFile = my_data -> p_VcfFile;
   		map<string,long long int*> *gene_pos_map = my_data -> gene_pos_map;
   		
   		string ref_file = my_data -> ref_file; 
   		string eqtl_path = my_data -> eqtl_path; 
   		string out_dir = my_data -> out_dir;
   		long long int start = my_data -> start; 
   		long long int end = my_data -> end;
   		int chrom = my_data -> chrom;
		
		ifstream fin(eqtl_path.c_str());
		fin.seekg(start);
		string line;
		vector<snps> origin_typed_snps;
		vector<snps> *p_origin_typed_snps = &origin_typed_snps;
		vector<typed_snp> typed_snps;
		vector<typed_snp> *p_typed_snps = &typed_snps;
		vector<snps> ignore_snps;
		vector<snps> *p_ignore_snps = &ignore_snps;
		vector<ref_snp> snp_map;
		vector<ref_snp> *p_snp_map = &snp_map;
	
		vector<string> hap;
		vector<string> *p_hap = &hap; 
	 
		//start extract
		getline(fin,line);
		string line_list[9];
		split_line(line_list,line);
		string gene_name = line_list[0];
		string last_gene_name;
		//record the gene
	
		map<long long int , int> m_all_snps ;
		map<long long int , int> m_typed_snps ;
		map<long long int , int> *p_m_all_snps = &m_all_snps;
		map<long long int , int> *p_m_typed_snps = &m_typed_snps;
		double** last_sigma_it = NULL;
		ans values;
		values.last_sigma_it = NULL;
		values.weight = NULL;
		int flag = 0;
		while(fin.tellg() <= end)
		{
			last_gene_name = gene_name;
			getline(fin,line);
			if(line == "")
			{
				flag = 1;
			}
			init_line_list(line_list);
			split_line(line_list,line);
			gene_name = line_list[0];

			if(last_gene_name != gene_name || 1 == flag)
			{
			//calculate window
				long long int window[2]; 
				calculate_window(window , last_gene_name , gene_pos_map);
			
			// 1.split origin_typed_snps into typed_snps
			//and unuseful snps(wrong and special type)
				filter_snps(ref_file , last_gene_name,p_origin_typed_snps ,
						 p_typed_snps , p_ignore_snps,window , p_VcfIndex,p_VcfFile); 
			//2.search gene annotation map and pos file dict to generate 
			//.map and .hap vector with .vcf
				gen_map_hap(ref_file , last_gene_name , p_snp_map , p_hap , window,p_VcfIndex,p_VcfFile);
			//3.impute process		
				size_t num_typed_snps = typed_snps.size();
				size_t num_total_snps = snp_map.size();
				vector<char> convert_flags;
				convert_flags.resize(num_typed_snps, 0);
				vector<char> impute_flags;
				impute_flags.resize(num_total_snps, 1);
				mark_snps(typed_snps, snp_map, convert_flags, impute_flags);
    			// convert z-score
				for(size_t i = 0; i < num_typed_snps; i++) 
				{
					if(convert_flags[i] == 1) 
					{
						typed_snps[i].zscore *= -1.0;
					}
				}	
				vector<long long int> useful_typed_snps;
				vector<long long int>* p_useful_typed_snps = &useful_typed_snps;
				vector<int> snps_flag;
				vector<int>* p_snps_flag = &snps_flag;	
				values = zgenbt(p_snps_flag ,typed_snps, snp_map , convert_flags , impute_flags , 
											hap , p_useful_typed_snps,
				p_m_all_snps,p_m_typed_snps,values.last_sigma_it);

				impz(chrom , out_dir , last_gene_name , p_snps_flag , values.weight ,typed_snps, snp_map , impute_flags ,p_useful_typed_snps );	
				//final.start new recorder and record gene_name snp
				clean_all_vector(p_origin_typed_snps , p_typed_snps ,
								p_ignore_snps , p_snp_map , p_hap ,p_useful_typed_snps);
				
//				cout << last_gene_name << " recorded" << endl;
				if(1 == flag)
				{
					break;
				}
				record_all(p_origin_typed_snps , line_list[1] , line_list[7] , line_list[8]);				
			} 
			else
			{
				record_all(p_origin_typed_snps , line_list[1] , line_list[7] , line_list[8]);
		
				//record it into origin_typed_snps
			}
		} 
		fin.close();
		sem_post(my_data -> bin_sem); //信号量+1

}

void make_output_dir(char *Out)
{
	string path = string(Out);
	char tem[10];
	for(int i = 1;i <= 22;i++)
	{
		string new_path = path;
		sprintf(tem , "%d" , i);
		string new_path1 = (new_path += string(tem));
		mkdir(new_path1.c_str() , 0777);	
	}
	
}

int main(int argc , char *argv[]) 
{
	char opt;
	char* Gene_annotation;
	char* Eqtl_path;
	char* Vcf_prefix;
	char* Out;
	char* BATCH;
	while((opt = getopt(argc,argv,"g:e:v:o:t:")) != -1) 
	{
		switch(opt) 
		{
			case 'g': // prefix used in gen_beta
				Gene_annotation = (char *) strdup(optarg);
				break;
			case 'e': // snp mapping file
				Eqtl_path = (char *) strdup(optarg);
				break;
			case 'v': // typed snp file
				Vcf_prefix = (char *) strdup(optarg);
				break;
			case 'o': // output file name
				Out = (char *) strdup(optarg);
				break;
			case 't':
				BATCH = (char*)strdup(optarg);
				break;		
		}
	}
	
	string gene_annotation = string(Gene_annotation);
//	string gene_annotation = "/media/userdisk1/jjpeng/yinquanwei/gencode_v19_gene_annotation.txt";
	string eqtl_path = string(Eqtl_path);
//	string eqtl_path = "/media/userdisk1/jjpeng/yinquanwei/Brain_Amygdala.allpairs.txt";
	string vcf_prefix = string(Vcf_prefix);
//	string vcf_prefix = "/media/userdisk1/jjpeng/yinquanwei/";
	string out = string(Out);
//	string out = "/media/userdisk1/jjpeng/yinquanwei/output/";
	batch = atoi(BATCH);
	make_output_dir(Out);

	//load gene_pos_map
	map<string,long long int*> pos1;
	map<string,long long int*> *gene_pos_map = &pos1;
	printf("Reading gene annotation file......\n"); 
	load_gene_pos_map(gene_annotation , gene_pos_map);
	cout << "Gene pos information loaded" << endl;		
		
	//seperate chrom
	long long int chrom[24];
	split_chrom(eqtl_path , chrom);
		

	for(int i = 1;i <= 22;i++)
	{
		long long int start = chrom[i - 1];
		long long int end = chrom[i];
		
		if(start == 0 || end == 0)
		{
			continue;
		}
			
		//load pos_file_map
		printf("reading vcf file of NO %d chrom......\n",i);
		map<string,long long int> VcfIndex;
		map<string,long long int>* p_VcfIndex = &VcfIndex;
		vector<string> VcfFile;
		vector<string>* p_VcfFile = &VcfFile;
		
		char tem[3];
		sprintf(tem , "%d" , i);
		string ref_file = vcf_prefix + "EUR.MAF01.Biallele.chr" + string(tem) + ".PASS_only.vcf.gz";
		gzLoadVcfFile(ref_file.c_str() ,p_VcfIndex ,p_VcfFile );
		cout << "VcfFile loaded" << endl;
				
		
		/////////////////////////////////////////
			
		vector<long long int> batch_bonder;
		vector<long long int>* p_batch_bonder = &batch_bonder;
		travel_eqtl(start ,end ,eqtl_path ,p_batch_bonder);
				
		///////////////////////////////////////// 
		
		sem_t  bin_sem;    //set Semaphore
		int res = sem_init(&bin_sem, 0, 0);//init Semaphore
		if (res != 0)
    	{
        	perror("Semaphore initialization failed");
    	}
		
		pthread_t tids[real_batch];
		struct thread_data td[real_batch];

		int chrom = i;
		printf("Imputing the %dst chrom......\n",chrom);
		for(int ii = 0;ii < real_batch;ii++)
		{
			td[ii].bin_sem = &bin_sem; 
			td[ii].gene_pos_map = gene_pos_map;
			td[ii].p_VcfIndex =  p_VcfIndex;
			td[ii].p_VcfFile =  p_VcfFile;
			td[ii].eqtl_path = eqtl_path;
			td[ii].ref_file = ref_file;
			td[ii].start = batch_bonder[2 * ii];
			td[ii].end = batch_bonder[2 * ii + 1];
			td[ii].out_dir = out;
			td[ii].chrom = chrom;

			int ret = pthread_create(&tids[ii], NULL, main_process, (void *)&td[ii]);
			if (ret != 0)
        	{
           		cout << "pthread_create error: error_code=" << ret << endl;
       		}	 
		}
	//防止循环结束时子程序还没有结束 
		for(int j = 0;j < real_batch;j++)
    	{
        	sem_wait(&bin_sem);    //循环等待次，相当于个线程都执行完了
    	}
    	sem_destroy(&bin_sem);        //释放信号量
		printf("the %dst chrom finished\n",chrom);
	}
	cout << "finish imputation\n";
	
	return 0;
}
