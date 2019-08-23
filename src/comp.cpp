#include "comp.h"
#include <algorithm>
 bool my_cmp(file_record a,file_record b)
{
	if(a.chrom == b.chrom)
	{
		if(a.l_win == b.l_win)
		{
			if(a.r_win == b.r_win)
			{
				return strcmp(a.name.c_str(),b.name.c_str());
			}
			return a.r_win < b.r_win;	
		}

       		 return a.l_win < b.l_win;
		
	}
	return a.chrom < b.chrom;


}

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
	int len = line.length();
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

	if((*pos_map).count(last_name) == 0)
	{
		cout << "Can't find record about "+ last_name+ " in annotation file!\n";
		exit(0);
	}
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
	for(int i = 0;i < chrom.length();i++)
	{
		if(chrom[i] >= '0' && chrom[i] <= '9')
		{

		}
		else
		{
			cout << "Unrecognized chromosome number " + chrom  << " in xQTL file"<< endl;
			cout << "Chromosomes should be named by number." << endl;
			exit(0);
		}
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
	ifstream fin(Xqtl_path.c_str());
	if(fin.fail())
	{
		cout << "Can't find xQTL file\n";
		exit(0);
	}
	string line;
	if(!getline(fin , line))
	{
		cout << "Empty xQTL file!\n";
		exit(0);
	}
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
			for(int i = current_chrom + 1;i < now_chrom;i++)
			{
				chrom[i] = chrom[i - 1];
			}
			current_chrom = now_chrom;	
		} 
	}
		

	return current_chrom;
	
}

void make_output_dir(long long int chrom[] , int chrom_num , char *Out , int chr)
{
	string path = string(Out);
	char tem[10];
	if(chr != -1)
	{
		string new_path = path;
		sprintf(tem , "%d" , chr);
		string new_path1 = (new_path + '.' + string(tem));
		mkdir(new_path1.c_str() , 0777);	
			
	}
	else
	{
		for(int i = 1;i <= chrom_num;i++)
		{
			if(chrom[i] == chrom[i - 1])
			{
				continue;
			}
			string new_path = path;
			sprintf(tem , "%d" , i);
			string new_path1 = (new_path + '.' + string(tem));
			if(NULL == opendir(new_path1.c_str()))
			{
				mkdir(new_path1.c_str() , 0777);	
			}
		}


	}

	
}

void organize_files(vector<string>* p_files , int chrom , string out , map<string,long long int*> pos_map)
{
		int cursor = 0;
		string new_path = out;
		char tem[10];
		sprintf(tem , "%d" , chrom);
		string new_path1 = out + '.' + string(tem) + "/";
		string out_file = out + "chr"+ string(tem) + "_zscores.txt";
	//	cout << out_file << endl;
		FILE* fp = fopen(out_file.c_str() , "w");
		fprintf(fp, "Chr Molecular_ID Molecular_Start Molecular_End Variant_ID Variant_pos Variant_Ref Variant_Alt Z-Statistic R2pred Imputation_flag\n");


		for(int i = 0 ; i < (*p_files).size() ; i++)
		{
			string file_path = new_path1 + (*p_files)[i];
			ifstream fin(file_path.c_str());
			string line;
			getline(fin,line);
			getline(fin,line);
			long long int start = pos_map[(*p_files)[i]][0];
			long long int end = pos_map[(*p_files)[i]][1];
			while(line != "")
			{
				fprintf(fp, "%s %s %lld %lld %s\n"  ,tem , (*p_files)[i].c_str() , start , end , line.c_str());
				getline(fin , line);
			}
			fin.close();
			remove(file_path.c_str());			
		}
		
		

		fclose(fp);
		rmdir(new_path1.c_str());



}

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: xQTLImp options\n ");
	fprintf (stream,
	"-h --help\t\t null, to display this usage. \n"
	"-x --xQTL\t\t string, the file path of xQTL summary statistics.\n"
	"-m --molecule\t\t string, the file path of molecule annotation file.\n"
	"-v --VCF files\t\t string, the folder path of genome reference panel, such as 1000G VCF files.\n"
	"-o --output\t\t string, the folder path of output results. \n"
	"-t --num_threads\t int, number of threads, 1 in default.\n"
	"-f --MAF_cutoff\t\t double, minimum MAF threshold for variants in genome reference panel, 0.01 in default.\n"
	"-l --lambda_value\t double, a constant value used to added with var-covariance matrix to gurantee the matrix is invertible, 0.1 in default.  \n"
	"-c --chr\t\t int, specify which chromosome will be imputed.\n"
	"-w --window_size \t int, Window size N, +-N/2 apart from molecular center pos, in base pair, N=500000bp in default.\n"
	"-e --exclude\t\t int:int-int, specify a genome region in which variants will be ignored during imputation process.\n"
	"-b --exclude_file\t string, multiple genome regions user want to mask during imputation process.\n"
	"-s --sort\t\t sort the xQTL summary statistics by chromosome, molecular_ID, and variant_pos in increasing order prior imputation (required), FALSE in default");

	exit (exit_code);
}




