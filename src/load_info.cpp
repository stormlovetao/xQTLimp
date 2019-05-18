#include "load_info.h"
#include "comp.h"


bool gzLoadVcfFile(char* tem , const char* gzfn , map<string,long long int>* p_VcfIndex ,vector<string>* p_VcfFile )
{

	struct dirent* ent = NULL;
   	 DIR *pDir;
   	 pDir=opendir(gzfn);
	gzFile gzfp = NULL;
    	while (NULL != (ent=readdir(pDir)))
    	{
  
 	           if (ent->d_type==8)
 		 {
			int i = 0;
			for(i = 0;i < strlen(tem);i++)
			{
				if(tem[i] == ent -> d_name[4 + i])
				{
				}
				else
				{
					break;
				}
			}
			if(i == strlen(tem) && ent -> d_name[4 + strlen(tem)] == '.')
			{
				string tem = string(gzfn) +  string(ent -> d_name);
				gzfp = gzopen(tem.c_str(),"rb");
				break;
			}
 	           }
            	else
            	{
            	}
      
  	}
	//open .gz file
	if(!gzfp)
	{
		cout << "\ncan't find reference panel!\n";
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

//record postions of all molecular
void load_pos_map(string annotation , map<string,long long int*> *pos_map)
{
	ifstream fin(annotation.c_str());  
	if(fin.fail())
	{
		cout << "\nCan't find " + annotation << endl;
		exit(0);
	} 
	string line;
	if(!getline(fin,line))
	{
		cout << "\nEmpty annotation file!\n";
		exit(0);
	}

	if(!getline(fin,line))
	{
		cout << "\nEmpty annotation file!\n";
		exit(0);
	}
	
	while(line != "")
	{
		//pos0:start pos  pos1:end pos
		long long int* pos = (long long int*)malloc(2 * sizeof(long long int));
		string name = "";
		string* cur_name = &name;
		//extract info from the line
		extract(pos , cur_name ,line);
		//record it into map
		(*pos_map)[name] = pos;
		getline(fin,line);
	}
	fin.close();
}


void record_all(vector<snps> *p_origin_typed_snps , string  pos , string ref, string alt  ,   string z)
{
	snps temp;
	temp.snp_pos = atoi(pos.c_str());
	temp.ref_allele = ref;
	temp.alt_allele = alt;
	temp.zscore = double(atof(z.c_str()));
	(*p_origin_typed_snps).push_back(temp);	
	
	
}

void filter_snps(string file_name , string last_name , vector<snps> *p_origin_typed_snps ,
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
void read_hap(string line,vector<string>* p_hap)
{		
	int index = 0;
	int len = line.length();
	for(int i = 0;i < 9;i++)
	{
		while(line[index] != '\t' && line[index] != ' ')
		{
			index++;
		}
		while(line[index] == '\t' || line[index] == ' ')
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

void gen_map_hap(string ref_file , string last_name ,  vector<ref_snp> *p_snp_map , 
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
int travel_Xqtl(long long int start , long long int end ,
			string Xqtl_path ,vector<long long int>* p_batch_bonder , int batch)
{
	int real_batch = 0;
	ifstream fin(Xqtl_path.c_str());
	string line;
	fin.seekg(start);
	int counter = 1;
	getline(fin,line);
	string line_list[100];
	split_line(line_list,line);
	string name = line_list[1];
	string last_name;
	while(fin.tellg() != end)
	{
		last_name = name;
		
		getline(fin,line);
		if(line == "")
		{
			break;
		}
		init_line_list(line_list);
		split_line(line_list,line);
		name = line_list[1];
		
		if(name == last_name)
		{
			
		}
		else
		{
			counter++;
		}
	}

	fin.close();
	///////////////////////////////////////////////////
	
	ifstream fin1(Xqtl_path.c_str());
	int batch_size = counter / batch;
	if(batch_size == 0)
	{
		real_batch = counter;
	}
	else
	{
		if(counter % batch == 0)
		{
			real_batch = batch;
		}
		else
		{
			real_batch = batch + 1;
		}	
	}
	
	(*p_batch_bonder).push_back(start);
		
	fin1.seekg(start);
	int num = 1;
	getline(fin1,line);
	string line_list1[100];
	split_line(line_list1,line);
	name = line_list1[1];
	long long int pos;
	int b = 0;
	long long int position = 0;
	while(true)
	{
		last_name = name;
		if(fin1.tellg() == end)
		{
			pos = fin1.tellg();
			(*p_batch_bonder).push_back(pos);
			break;
		}
		
		position = fin1.tellg();
		getline(fin1 , line);
		init_line_list(line_list1);
		split_line(line_list1,line);
		name = line_list1[1];
		
		if(name == last_name)
		{
			
		}
		else
		{
			num++;
			if((num - 1) % batch_size == 0)
			{
			    b++;
				(*p_batch_bonder).push_back(position);
				(*p_batch_bonder).push_back(position);			
			}
		}
	}
	
	fin1.close();
	return real_batch;
	
}



