#include <pthread.h> 
#include <semaphore.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "impute.h"
#include "load_info.h"
#include "comp.h"
#include "utils.h" 


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
  double maf;
  double lam;
};

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
		double maf = my_data -> maf;
		double lam  = my_data -> lam;
		
		ifstream fin(eqtl_path.c_str());
		fin.seekg(start);
		string line;
		vector<snps> origin_typed_snps;
		vector<snps> *p_origin_typed_snps = &origin_typed_snps;
		vector<typed_snp> typed_snps;
		vector<typed_snp> *p_typed_snps = &typed_snps;
		vector<typed_snp> maf_snps;
		vector<typed_snp> *p_maf_snps = &maf_snps;
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
		string gene_name = line_list[1];
		string last_gene_name;
		//record the gene
	
		map<long long int , int> m_all_snps ;
		map<long long int , int> m_typed_snps ;
		map<long long int , int> *p_m_all_snps = &m_all_snps;
		map<long long int , int> *p_m_typed_snps = &m_typed_snps;
		
		ans values;
		values.last_sigma_it = NULL;
		values.weight = NULL;
		values.yaoyao = 0;
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
			gene_name = line_list[1];

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
				values = zgenbt(p_maf_snps , maf , lam , p_snps_flag ,typed_snps, snp_map , convert_flags , impute_flags , 
											hap , p_useful_typed_snps,
				p_m_all_snps,p_m_typed_snps,values.last_sigma_it,values.yaoyao);

				impz(p_maf_snps , p_ignore_snps , chrom , out_dir , last_gene_name , p_snps_flag , values.weight ,typed_snps, snp_map , impute_flags ,p_useful_typed_snps );	
				//final.start new recorder and record gene_name snp
				clean_all_vector(p_maf_snps , p_origin_typed_snps , p_typed_snps ,
								p_ignore_snps , p_snp_map , p_hap ,p_useful_typed_snps , p_snps_flag);
				
				if(1 == flag)
				{
					break;
				}
				record_all(p_origin_typed_snps , line_list[2] , line_list[3] , line_list[4] , line_list[5]);				
			} 
			else
			{
				record_all(p_origin_typed_snps , line_list[2] , line_list[3] , line_list[4] , line_list[5]);
		
				//record it into origin_typed_snps
			}
		} 
		fin.close();
		
		for(int i = 0;i < (values.yaoyao);i++)
		{
			free(values.last_sigma_it[i]);
		}
		free(values.last_sigma_it);
		m_all_snps.clear();
		m_typed_snps.clear();
		
		sem_post(my_data -> bin_sem); //信号量+1

}

int main(int argc , char *argv[]) 
{
	char opt;
	char* Gene_annotation;
	char* Eqtl_path;
	char* Vcf_prefix;
	char* Out;
	char* BATCH = NULL;
	char* MAF = NULL;
	char* LAMBDA = NULL;

	while((opt = getopt(argc,argv,"m:x:v:o:t:f:l:")) != -1) 
	{
		switch(opt) 
		{
			case 'm': // prefix used in gen_beta
				Gene_annotation = (char *) strdup(optarg);
				break;
			case 'x': // snp mapping file
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
			case 'f':
				MAF = (char*)strdup(optarg);
				break;
			case 'l':
				LAMBDA = (char*)strdup(optarg);
				break;
					
		}
	}
	
	string gene_annotation = string(Gene_annotation);
	string eqtl_path = string(Eqtl_path);
	string vcf_prefix = string(Vcf_prefix);
	string out = string(Out);
	int batch = 1;
	double maf = 0.01;
	double lam = 0.1;

	if(BATCH == NULL)
	{
		batch = 1;
	}
	else
	{
		batch = atoi(BATCH);
	}

	if(MAF == NULL)
	{
		maf = 0.01;
	}
	else
	{
		maf = double(atof(MAF));
	}

	if(LAMBDA == NULL)
	{
		lam = 0.1;
	}
	else
	{
		lam = double(atof(LAMBDA));
	}



	make_output_dir(Out);

	//load gene_pos_map
	map<string,long long int*> pos1;
	map<string,long long int*> *gene_pos_map = &pos1; 
	cout << "Reading molecular annotation file..." << endl;	
	load_gene_pos_map(gene_annotation , gene_pos_map);	
	cout << "Done!" << endl;
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
		string ref_file = vcf_prefix;
		gzLoadVcfFile(tem , ref_file.c_str() ,p_VcfIndex ,p_VcfFile );
		cout << "VcfFile loaded" << endl;
				
		
		/////////////////////////////////////////
			
		vector<long long int> batch_bonder;
		vector<long long int>* p_batch_bonder = &batch_bonder;
		int real_batch = travel_eqtl(start ,end ,eqtl_path ,p_batch_bonder , batch);
				
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
			td[ii].maf = maf;
			td[ii].lam = lam;

			int ret = pthread_create(&tids[ii], NULL, main_process, (void *)&td[ii]);
			if (ret != 0)
        	{
           		cout << "pthread_create error: error_code=" << ret << endl;
       		}	 
		}
	//防止循环结束时子程序还没有结束 
		for(int j = 0;j < real_batch;j++)
    	{
        	sem_wait(&bin_sem);    //循环等待n次，相当于线程都执行完了
    	}
    	sem_destroy(&bin_sem);        //释放信号量
    	
		printf("the %dst chrom finished\n",chrom);
	}
	
	printf("organizing files......\n");
	
	organize_files(out);
	
	cout << "finish imputation\n";
	
	return 0;
}
