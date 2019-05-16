#include "comp.h"

string read_pos(string line)
{
	string pos = "";
	int i = 0;
	while(line[i] != '\t' && line[i] != ' ')
	{
		i++;
	}
	while(line[i] == '\t' || line[i] == ' ' )
	{
		i++;
	}
	while(line[i] != '\t' && line[i] != ' ')
	{
		pos += line[i];
		i++;
	}
	return pos;
	
} 

void split_line(string line_list[],string line)
{
	int flag = 0;
	int index = 0;
	
	for(int i = 0;i < line.length();i++)
	{
		flag = 0;
		if(line[i] == '\t' || line[i] == ' ' )
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
void extract(long long int* pos,string* name,string line)
{
	//read start pos
	int index = 0;
	while(line[index] != '\t' && line[index] != ' ' )
	{
		(*name) += line[index];
		index++;
	}

	while(line[index] == '\t' || line[index] == ' ')
	{
		index++;
	}
	string pos1 = "";
	while(line[index] != '\t'&& line[index] != ' ')
	{
		pos1 += line[index];
		index++;
	}
	while(line[index] == '\t'|| line[index] == ' ')
	{
		index++;
	} 

	pos[0] = atoi(pos1.c_str());
	//read end pos
	string pos2 = "";
	while( (index < line.length()) && (line[index] != '\t' && line[index] != ' ') )
	{
		pos2 += line[index];
		index++; 
	}
	pos[1] = atoi(pos2.c_str());
	

}
bool read_pos_name_alleles(string line , long long int* pos , string* name,string allele[2])
{
	int index = 0;
	while(line[index] != '\t' && line[index] != ' ')
	{
		index++;
	}
	while(line[index] == '\t' || line[index] == ' ')
	{
		index++;
	}
	
	
	string str_pos = "";
	while(line[index] != '\t' && line[index] != ' ' )
	{
		str_pos += line[index];
		index++;
	}
	while(line[index] == '\t'  || line[index] == ' ')
	{
		index++;
	}
	(*pos) = atoi(str_pos.c_str());
	
	(*name) = "";
	while(line[index] != '\t' && line[index] != ' ')
	{
		(*name) += line[index];
		index++;
	}
	while(line[index] == '\t' || line[index] == ' ')
	{
		index++;
	}
	
	string ref = "";	
	while(line[index] != '\t' && line[index] != ' ')
	{
		ref += line[index];
		index++;
	}
	while(line[index] == '\t' || line[index] == ' ')
	{
		index++;
	}
	
	string alt = "";
	while(line[index] != '\t' && line[index] != ' ')
	{
		alt += line[index];
		index++;
	}
	
	while(line[index] != '\t' || line[index] == ' ')
	{
		alt += line[index];
		index;
	}
	allele[0] = ref;
	allele[1] = alt;
	
	return true;
}

void calculate_window(int window_size , long long int window[2],string last_name ,
				map<string,long long int*> *pos_map)
{

	long long int l_pos = (*pos_map)[last_name][0];
	long long int r_pos = (*pos_map)[last_name][1];
	long long int pos = (l_pos + r_pos) / 2;
	long long int pos1 = pos - (window_size) / 2;
	if(pos1 < 0)
	{
		pos1 = 0;
	}
	window[0] = pos1;
	window[1] = pos + (window_size) / 2; 
					
}

void clean_all_vector(      		vector<typed_snp> *p_maf_snps,	
							vector<snps> *p_origin_typed_snps , 
							vector<typed_snp> *p_typed_snps ,
							vector<snps>*p_ignore_snps , 
							vector<ref_snp>*p_snp_map,
							vector<string>* p_hap,
							vector<long long int>* p_useful_typed_snps,
								vector<int>* p_snps_flag)
{
	(*p_maf_snps).clear();
	(*p_maf_snps).shrink_to_fit();
	(*p_origin_typed_snps).clear();
	(*p_origin_typed_snps).shrink_to_fit();
	(*p_typed_snps).clear();
	(*p_typed_snps).shrink_to_fit();
	(*p_ignore_snps).clear();
	(*p_ignore_snps).shrink_to_fit();
	(*p_snp_map).clear();
	(*p_snp_map).shrink_to_fit();
	(*p_hap).clear();
	(*p_hap).shrink_to_fit();
	(*p_useful_typed_snps).clear();
	(*p_useful_typed_snps).shrink_to_fit();
	(*p_useful_typed_snps).clear();
	(*p_useful_typed_snps).shrink_to_fit();
								
}
int get_chrom(string line)
{
	int index = 0;
	string chrom = "";
	while(line[index] != '\t' && line[index] != ' ')
	{
		chrom += line[index];
		index++;
	}

	int ans = atoi(chrom.c_str());
	return ans;
	
}

int split_chrom(string Xqtl_path ,long long int chrom[])
{
	for(int i = 0;i <= 1260;i++)
	{
		chrom[i] = 0;
	}
	cout << "Scaning xQTL file...\n" ;
	ifstream fin(Xqtl_path.c_str());
	string line;
	getline(fin , line);
	//erase head
	chrom[0] = fin.tellg();
	int current_chrom = 1;
	long long int pos = 0;
	while(true)
	{
		pos = fin.tellg();
		getline(fin , line);
		if(line == "")
		{
			chrom[current_chrom] = pos;
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
			current_chrom++;	
		} 
	}
		

	return current_chrom;
	
}

void make_output_dir(int chrom_num , char *Out)
{
	string path = string(Out);
	char tem[10];
	for(int i = 1;i <= chrom_num;i++)
	{
		string new_path = path;
		sprintf(tem , "%d" , i);
		string new_path1 = (new_path + string(tem));
		mkdir(new_path1.c_str() , 0777);	
	}
	
}

void organize_files(int chrom_num , string out , map<string,long long int*> pos_map)
{
	
	for(int i = 1;i <= chrom_num;i++)
	{
		string new_path = out;
		char tem[10];
		sprintf(tem , "%d" , i);
		string new_path1 = out + string(tem) + "/";
		string out_file = out + "chr"+ string(tem) + "_zscores.txt";
	//	cout << out_file << endl;
		FILE* fp = fopen(out_file.c_str() , "w");
		fprintf(fp, "Chr Molecular_ID Molecular_Start Molecular_End Variant_ID Variant_pos Variant_Ref Variant_Alt Z-Statistic R2pred Imputation_flag\n");
		struct dirent* ent = NULL;
   	 	DIR *pDir;
   	 	pDir=opendir(new_path1.c_str());
    		while (NULL != (ent=readdir(pDir)))
    		{
  
 	           	if (ent->d_type==8)
 			 {
				string file_path = new_path1 + string(ent -> d_name);
				ifstream fin(file_path.c_str());
				string line;
				getline(fin,line);
				getline(fin,line);
				long long int start = pos_map[string(ent -> d_name)][0];
				long long int end = pos_map[string(ent -> d_name)][1];
				while(line != "")
				{
				    fprintf(fp, "%s %lld %lld %s %s\n"  ,tem , ent -> d_name , start , end , line.c_str());
				    getline(fin , line);
				}
				fin.close();
				remove(file_path.c_str());
				
			}
            		else
            		{
            		}
      
  		}

		fclose(fp);
		rmdir(new_path1.c_str());

	}

}

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: xQTLImp options\n ");
	fprintf (stream,
	"-h --help\t\t Display this usage information. \n"
	"-x --xQTL\t\t The file path of xQTL summary statistics.\n"
	"-m --molecule\t\t The file path of molecule annotation file.\n"
	"-v --VCF files\t\t The folder path of genome reference panel, such as 1000G VCF files.\n"
	"-o --output\t\t The folder path of output results.\n"
	"-t --num_threads\t Number of threads, 1 in default\n"
	"-f --MAF_cutoff\t\t Minimum MAF threshold for variants in genome reference panel, 0.01 in default.\n"
	"-l --lambda_value\t A constant value used to added with var-covariance matrix to gurantee the matrix is invertible, 0.1 in default \n"
	"-w --window_size \t Window size N, +-N/2 apart from molecular start pos, 500Kb in default.\n");

	exit (exit_code);
}




