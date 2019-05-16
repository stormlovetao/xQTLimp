#include <pthread.h> 
#include <semaphore.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include "impute.h"
#include "load_info.h"
#include "comp.h"
#include "utils.h" 


struct thread_data{
   map<string,long long int>* p_VcfIndex;
   vector<string>* p_VcfFile;
   map<string,long long int*> *pos_map;
   int chrom; 
   string Xqtl_path; 
   string ref_file;
   string out_dir;
   long long int start; 
   long long int end;
   sem_t * bin_sem;
  double maf;
  double lam;
  int window_size;
};

void *main_process(void *threadarg)
{
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
   	map<string,long long int>* p_VcfIndex = my_data -> p_VcfIndex;
   	vector<string>* p_VcfFile = my_data -> p_VcfFile;
   	map<string,long long int*> *pos_map = my_data -> pos_map;
   		
   	string ref_file = my_data -> ref_file; 
   	string Xqtl_path = my_data -> Xqtl_path; 
   	string out_dir = my_data -> out_dir;
   	long long int start = my_data -> start; 
   	long long int end = my_data -> end;
	int window_size = my_data -> window_size;
   	int chrom = my_data -> chrom;
	double maf = my_data -> maf;
	double lam  = my_data -> lam;
		
	ifstream fin(Xqtl_path.c_str());
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
	string line_list[100];
	split_line(line_list,line);
	string name = line_list[1];
	string last_name;
	record_all(p_origin_typed_snps , line_list[2] , line_list[3] , line_list[4] , line_list[5]);	
	//record the molecular
	
	map<long long int , int> m_all_snps ;
	map<long long int , int> m_typed_snps ;
	map<long long int , int> *p_m_all_snps = &m_all_snps;
	map<long long int , int> *p_m_typed_snps = &m_typed_snps;
		
		ans values;
		values.last_sigma_it = NULL;
		values.weight = NULL;
		values.row = 0;
		int flag = 0;
		while(fin.tellg() <= end)
		{
			last_name = name;
			getline(fin,line);
			if(line == "")
			{
				flag = 1;
			}
			init_line_list(line_list);
			split_line(line_list,line);
			name = line_list[1];

			if(last_name != name || 1 == flag)
			{
			//calculate window
				long long int window[2]; 
				calculate_window(window_size , window , last_name , pos_map);
			
			// 1.split origin_typed_snps into typed_snps
			//and unuseful snps(wrong and special type)
				filter_snps(ref_file , last_name,p_origin_typed_snps ,
						 p_typed_snps , p_ignore_snps,window , p_VcfIndex,p_VcfFile); 
			//2.search  annotation map and pos file dict to generate 
			//.map and .hap vector with .vcf
				gen_map_hap(ref_file , last_name , p_snp_map , p_hap , window,p_VcfIndex,p_VcfFile);
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
				p_m_all_snps,p_m_typed_snps,values.last_sigma_it,values.row);

				impz(p_maf_snps , p_ignore_snps , chrom , out_dir , last_name , p_snps_flag , values.weight ,typed_snps, snp_map , impute_flags ,p_useful_typed_snps );	
				//final.start new recorder
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
		
		for(int i = 0;i < (values.row);i++)
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
	char* Annotation;
	char* XQTL_path;
	char* Vcf_prefix;
	char* Out;
	char* BATCH = NULL;
	char* MAF = NULL;
	char* LAMBDA = NULL;
	char* WIN = NULL;

	const char *const short_options = "hx:m:v:o:t:f:l:w:";
	const struct option long_options[] = {
	{"help", 0, NULL, 'h'},
	{"xQTL", 1, NULL, 'x'},
	{"molecule", 1, NULL, 'm'},
	{"VCF", 1, NULL, 'v'},
	{"output", 1, NULL, 'o'},
	{"threads", 1, NULL, 't'},
	{"MAF_cutoff", 1, NULL, 'f'}, 
	{"lambda", 1, NULL, 'l'},
	{"window_size", 1, NULL, 'w'},
	{NULL, 0, NULL, 0 }
	};

	char *program_name;
	char input_filename[256], output_directory[256];

	int v = 0;
	strcpy(output_directory, "result");
	int next_option = 1;

	int counter = 0;
	program_name = argv[0];
	do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);

		switch (next_option)
		 {
			case 'h':
				print_usage (stdout, 0);
				break;

			case 'm': // Molecular annotation file
				Annotation = (char *) strdup(optarg);
				break;
			case 'x': // xQTL statistic summary
				XQTL_path = (char *) strdup(optarg);
				break;
			case 'v': // genome reference panel, such as 1000G VCF
				Vcf_prefix = (char *) strdup(optarg);
				break;
			case 'o': // output folder path
				Out = (char *) strdup(optarg);
				break;
			case 't': // threads number， 1 in default
				BATCH = (char*)strdup(optarg);
				break;	
			case 'f': // MAF threshold for variants in reference panel, 0.01 in default
				MAF = (char*)strdup(optarg);
				break;
			case 'l': // lambda value, 0.1 in default
				LAMBDA = (char*)strdup(optarg);
				break;
			case 'w': // window size in total(upstream and downstream of TSS), 500kb in default
				WIN = (char*)strdup(optarg);
				break;

			case '?':
				print_usage (stdout, 1);
				break;

			case -1:/* Done with options. */
				if(counter ==0)
				{
				        print_usage (stdout, 0);
				}
				break;

			default:/* Something else: unexpected. */
				print_usage (stderr, -1);
		 }
		counter++;
	} while (next_option != -1);
	
	string annotation = string(Annotation);
	string Xqtl_path = string(XQTL_path);
	string vcf_prefix = string(Vcf_prefix);
	string out = string(Out);
	int batch = 1;
	double maf = 0.01;
	double lam = 0.1;
	int window_size = 500000;

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

	if(WIN == NULL)
	{
		window_size = 500000;
	}
	else
	{
		window_size = int(atoi(WIN));
	}

	//load pos_map
	cout << "Loading molecular annotation file ... ";
	map<string,long long int*> pos1;
	map<string,long long int*> *pos_map = &pos1; 
	load_pos_map(annotation , pos_map);	
	 cout << "Done!" << endl;	
	//seperate chrom
	long long int chrom[1261];
	int chrom_num = split_chrom(Xqtl_path , chrom);
	cout << chrom_num << " chromosomes detected!" << endl;

	make_output_dir(chrom_num , Out);

	for(int i = 1;i <= chrom_num;i++)
	{
		long long int start = chrom[i - 1];
		long long int end = chrom[i];
		
		if(start == 0 || end == 0)
		{
			continue;
		}
			
		//load pos_file_map
		printf("Loading reference VCF file of chromosome No.%d ... ",i);
		map<string,long long int> VcfIndex;
		map<string,long long int>* p_VcfIndex = &VcfIndex;
		vector<string> VcfFile;
		vector<string>* p_VcfFile = &VcfFile;
		
		char tem[3];
		sprintf(tem , "%d" , i);
		string ref_file = vcf_prefix;
		gzLoadVcfFile(tem , ref_file.c_str() ,p_VcfIndex ,p_VcfFile );
		cout << "Done!\n";	
		
		/////////////////////////////////////////
			
		vector<long long int> batch_bonder;
		vector<long long int>* p_batch_bonder = &batch_bonder;
		int real_batch = travel_Xqtl(start ,end ,Xqtl_path ,p_batch_bonder , batch);
				
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
		printf("Performing xQTL imputation on chromosome No.%d ...\n",chrom);
		printf("Partitioning into %d threads ... \n" , batch);
		for(int ii = 0;ii < real_batch;ii++)
		{
			td[ii].bin_sem = &bin_sem; 
			td[ii].pos_map = pos_map;
			td[ii].p_VcfIndex =  p_VcfIndex;
			td[ii].p_VcfFile =  p_VcfFile;
			td[ii].Xqtl_path = Xqtl_path;
			td[ii].ref_file = ref_file;
			td[ii].start = batch_bonder[2 * ii];
			td[ii].end = batch_bonder[2 * ii + 1];
			td[ii].out_dir = out;
			td[ii].chrom = chrom;
			td[ii].maf = maf;
			td[ii].lam = lam;
			td[ii].window_size = window_size;

			int ret = pthread_create(&tids[ii], NULL, main_process, (void *)&td[ii]);
			if (ret != 0)
        	{
           		cout << "pthread_create error: error_code=" << ret << endl;
       		}	 
		}
	// To make sure all sub-threads are finished in each chromosome, because they share the same memory of the reference panel
		for(int j = 0;j < real_batch;j++)
    	{
        	sem_wait(&bin_sem);    
    	}
    	sem_destroy(&bin_sem);        //release sem
    	
		printf("Imputation on chromosome No.%d finished!\n",chrom);
	}
	
	printf("Finalizing output files ... ");
	
	organize_files(chrom_num , out , pos1);
	cout << "Done!" << endl;
	
	cout << "Imputation processes have been Successfully Finished!\n";
	
	return 0;
}
