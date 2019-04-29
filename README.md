# xQTLImp
## introduction
xQTLImp is an open source software that implements xQTL(such as eQTL, mQTL, haQTL et. al) statistics (i.e Z statistics) imputation across the genome. xQTLImp accepts xQTL summary statistics without the need of individual-level genotypes and molecular features, and could accurately impute novel xQTL associations based on genetic reference panel. The imputation process is performed by modeling variants, within same LD and associated with a certain molecular trait, by using multivariate normal distribution. Novel QTL statistics of variants associated with a molecular trait, say gene *G's* expression level, will be calculated as a weighted linear combination of known statistics of variants associated with G, with the weights reflecting LD relationships (i.e. LD r^2) among those variants (See Methods of Paper). Genetic reference panel such as 1000G, HapMap or user-provided refrence panel in gzipped VCF format is required for LD calculations.

In general, xQTLImp can handle any kinds of molecular traits (not limited to gene expression, DNA methylation or histone acetylation) that can be physically mapped onto genomic regions, without limitation in species. Genomic variants such as SNV and small Indel(Insertion/deletion) are both supported. And xQTLImp can be applied in multiple scenarios, such as in performing meta-analysis of multiple eQTL studies which have reported different, but correlated sets of variants. 

xQTLImp is implemented in C++ and released under GNU GPL license. The source code and sample data are freely available for download at current webpage, and can be run on Linux/Unix/windows with C++ environment. To get better performance, xQTLImp can be executed in parallel mode, during which each chromosome will be broken into *N* chunks (*N* = number of threads) with each chunk has similar number of molecules. And on average, xQTLImp will need 8Gb memory to execute, which is easy to be distributed on PC and server.
</br>

##  Building xQTLImp
 To read gzipped VCF files, **zlib** package is required, and to install:
```bash
sudo apt-get install zlib1g-dev
```
Then compiling:
```bash
cd ./src/
make #Any C++11 compiler should work.
```
## Usage
### Requirements for input files.
The following files and format are required as input:
#### 1. Genome reference panel in gzipped VCF format. 
For example, the 1000G human genome reference panel [(Available here)](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/), or HapMap3 reference panel [(Available here)](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html).</br>

The genotype VCF files should be separated by chromosomes, and named in the format of 'chr.*N*.XXX.vcf.gz' within a folder. *N* represents for chromosome number, i.e. 1~26 (X->23, Y->24, XY (Pseudo-autosomal region of X) ->25, MT (Mitochondrial) ->26); XXX represents for any user defined string. Other domain should be freezed in required format.  </br>
</br>
#### 2.Molecular trait file
Molecular trait file must start with a line that contains column labels – molecular_ID, start_pos, end_pos ...(option) </br>followed by lines of data entries. Each field of data entries must be separated by white spaces.</br>
##### Example:
`molecular_ID` `start_pos` `end_pos`</br>
ENSG00000223972.4	11869	14412</br>
ENSG00000227232.4	14363	29806</br>
ENSG00000243485.2	29554	31109</br>
ENSG00000237613.2	34554	36081</br>
......</br>
</br>
#### 3.xQTL file
xQTL file must start with a line that contains column labels – chromosome , molecular_ID, variant_pos , Ref_allele , Alt_allele , z_statistics ...(option)</br> followed by lines of data entries. Each field of data entries must be separated by white spaces.</br>
##### Example:
`chromosome` `molecular_ID` `variant_pos` `Ref_allele`  `Alt_allele` `z_statistics`</br>
1 ENSG00000223972.4 13417 G C 1.5</br>
1 ENSG00000227232.4 17559 A G 2.6</br>
1 ENSG00000227232.4 54421 G A -1.0</br>
......</br>
</br>
### Parameter Description：
-m : the path of Molecular trait file</br>
-x : the path of xQTL file</br>
-v : the 1000G files folder</br>
-o : the output folder</br>
-t : the num of threads</br>
</br>
### demon
There is a demon in sample folder.
#### step1: make under src folder
#### step2: create a new folder for output
#### step3: Execute the command line under src folder
```bash
./xQTLImp -m /sample/gencode_v19_gene_annotation.txt -x /sample/Brain_Amygdala.allpairs.txt -v /sample/ -o (your output folder) -t 6
```
### the output file format
22 subfolders will appear in the output folder,and each folder corresponds to the imputing result on the chromosome.</br>
A file corresponds to the impute result of a gene</br>
#### Example:
`SNP_name` `SNP_pos` `Ref_Allele` `Alt_Allele` `Z-Score` `r2pred`</br>
rs142006308 31757791 G A 0.328344 1.000000</br>
rs71563368 31758240 G A 0.335234 1.000000</br>
rs6899983 31758931 A C 0.279963 1.000000</br>
......</br>

## Reference
During the implementation of xQTLImp, we referenced the work of:
* Pasaniuc, B. et al., 2014. Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), pp.2906–2914.
* Kwan, J.S.H. et al., 2016. FAPI: Fast and Accurate P-value Imputation for genome-wide association study. European Journal of Human Genetics, 24(5), p.761.
* Han, B., Kang, H.M. & Eskin, E., 2009. Rapid and accurate multiple testing correction and power estimation for millions of correlated markers. PLoS genetics, 5(4), p.e1000456.




