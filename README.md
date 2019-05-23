# xQTLImp
## Introduction
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
For example, the 1000G human genome reference panel [(Available here)](http://www.internationalgenome.org/data/), or HapMap3 reference panel [(Available here)](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html).</br>

```bash
# Codes for downloading 1000G Chr1 to Chr22 VCF files.
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;
for chr in {1..22} ; do
    wget $prefix$chr$suffix  $prefix$chr$suffix.tbi ;
done
```
The genotype reference VCF files should be separated by chromosomes, and named in the format of 'chr.*N*.XXX.vcf.gz' within a folder. *N* represents for chromosome number, for example 1~26 (X->23, Y->24, XY (Pseudo-autosomal region of X) ->25, MT (Mitochondrial) ->26); XXX represents for any user defined string. Other domain should be freezed in required format. See *Tricks* if you do not want to change your file name. </br>
</br>

*Notes:* To get better performance and accuracy, we recommend the users to preprocess the genotype reference VCF files. For example, keeping only samples of EUR population if the input xQTL statistics are from subjects of European ancestry, or/and filtering rare and non-biallelic variants since most of the current xQTL studies still focus on common variants due to limited sample size. The example source codes with *[VCFtools](https://github.com/vcftools/vcftools)* are as follows:

```bash
for chr in {1..22} ; do
vcftools --gzvcf $prefix$chr$suffix \
--maf 0.01                          \ # keep common variants with MAF > 0.01
--min-alleles 2 --max-alleles 2     \ # keep biallelic variants
--keep $EUR_file                    \ # keep EUR population samples, not provided.
--remove-filtered-all --recode --stdout | gzip -c > $out_prefix$chr$out_suffix;
done
```

#### 2. Molecule annotation file
Molecule annotation file gives the physical position of each molecule on reference genome. This annotation file should start with a column name line that contains at least three columns – molecular_ID, start_pos, end_pos ...(optional columns), followed by lines of data entries. Each field of data entries must be separated by tab or white spaces. Header names can be changed.</br>
##### Example:

`molecular_ID`| `start_pos`| `end_pos`
------------------|-------|-------
ENSG00000223972.4	| 11869 |	14412
ENSG00000227232.4	| 14363	| 29806
ENSG00000243485.2	| 29554	| 31109
ENSG00000237613.2	| 34554	| 36081
......</br>


#### 3. xQTL summary statistics
This file gives summary statistics associated with pairs of variants and molecular traits. The xQTL file should start with a line that contains at least 6 columns – chromosome , molecular_ID, variant_pos , ref_allele , alt_allele , z_statistics ...(optional columns),followed by lines of data entries. Each field of data entries must be separated by tab or white spaces. Header names can be changed. </br>
##### Example:
`chromosome` | `molecular_ID` | `variant_pos` | `ref_allele` | `alt_allele` | `z_statistics`
--|--|--|--|--|--
1 | ENSG00000223972.4 | 13417 | G | C    | 1.5
1 | ENSG00000227232.4 | 17559 | A | AGCC | 2.6
1 | ENSG00000227232.4 | 54421 | G | A    | -1.0
......</br>

*Notes:* </br>
* xQTLImp only handle cis-xQTL associations because of the nature of LD, so variant and molecule should be on the some chromosome. </br>
* The coordinates of variants and molecules should be in same genome build (e.g. hg19).
* The variant positions associated with each molecule should be in order.
* The Z statistic should be computed as the effect of the same type of allele (Ref or Alt) for all data entries (i.e. The effect allele is Ref allele or Alt allele).  And the effect allele of imputed variants will be the same type of effect allele as the inputs.
* The sign of Z statistics will be converted by xQTLImp if the user provided <ref, alt> alleles are in opposite to that in reference panel.
* Data entries (input lines) would not be used for imputation if xQTLImp cannot locate the variant on reference panel, but they will still be written into outputs without loss of information. 


### xQTLImp parameters：
The parameters can be specified by short options (e.g. -h) or long options (e.g. --help), which have same effect. Please see the candidate parameters and explanations as follows:
```bash
xQTLImp
-h, --help           null           # null, to display this usage.
-x, --xQTL           file_path      # string, the file path of xQTL summary statistics.
-m, --molecule       file_path      # string, the file path of molecule annotation file.
-v, --VCF            folder_path    # string, the folder path of genome reference panel, such as 1000G VCF files.
-o, --output         folder_path    # string, the folder path of output results. 
-t, --num_threads    num_threads    # int, number of threads, 1 in default.
-f, --MAF_cutoff     MAF_cutoff     # double, minimum MAF threshold for variants in genome reference panel, 0.01 in default.
-l, --lambda_value   lambda_value   # double, a constant value used to added with var-covariance matrix to gurantee the matrix is invertible, 0.1 in default. 
-w, --window_size    window_size    # int, Window size N, +-N/2 apart from molecular start pos, 500Kb in default.
```

### Running sample data
The sample data is generated from [GTEx](https://gtexportal.org/home/index.html) cis-eQTL results from Brain Amygdala downloaded in January 2019. We selected 49 genes with all associated variants (Pvalue < 1) from chromosome 1.

```bash
cd ./src ; make; cd .. # compliling under src folder
mkdir sample_output    # creat a new folder for output
# Running xQTLImp
./src/xQTLImp                                  \
-m ./sample/gencode_v19_gene_annotation.txt    \
-x ./sample/Brain_Amygdala.allpairs.sample.txt \
-v ./sample/                                   \
-o ./sample_output/                            \
-t 2                                           \
-w 500000
```
### xQTLImp output file format
*N* files corresponding with imputation results of *N* chromosomes will be created in the output folder.</br> 
#### Example:
`Chr` | `Molecular_ID` | `Molecular_Start` | `Molecular_End` | `Variant_ID` | `Variant_pos` | `Variant_Ref` | `Variant_Alt` | `Z-Statistic` | `R2pred` | `Imputation_flag`
--|--|--|--|--|--|--|--|--|--|--
1 | ENSG00000231709.1 | 521369 | 523833 | 1_13417_C_CGAGA | 13417 | C | CGAGA | -0.338210 | 0.98 | 1
1 | ENSG00000231709.1 | 521369 | 523833 | 1_17559_G_C | 17559 | G | C | 1.605512 | 1.000000 | 0
1 | ENSG00000231709.1 | 521369 | 523833 | 1_54421_A_G | 54421 | A | G | 1.069012 | 1.000000 | 0
......

## Reference
During the implementation of xQTLImp, we referenced the following works:
* Pasaniuc, B. et al., 2014. Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), pp.2906–2914.
* Kwan, J.S.H. et al., 2016. FAPI: Fast and Accurate P-value Imputation for genome-wide association study. European Journal of Human Genetics, 24(5), p.761.
* Han, B., Kang, H.M. & Eskin, E., 2009. Rapid and accurate multiple testing correction and power estimation for millions of correlated markers. PLoS genetics, 5(4), p.e1000456.


## Bug shooting
```bash
# Mac user might occur the following error during compiling:
# xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun
# which means you need to install proper C++ compiling environment. Installing XCode can fix this problem, open a Terminal and run this command:
xcode-select --install
```

## Tricks
If you do not want to rename the VCF file names of the reference panel, you can simply creat softlinks for them in our required format. For example:
```bash
# Source VCF files are in the path of:
${source_vcf_path}/EUR.MAF01.Biallele.chr[1..22].PASS_only.vcf.gz
# And you can create softlinks for those VCF files in a new path, and use the new folder path for input:
for chr in {1..22} ; do
	ln -s ${source_vcf_path}/EUR.MAF01.Biallele.chr${chr}.PASS_only.vcf.gz ${soft_links_path}/chr.${chr}.EUR_MAF01_Biallele_PASS_only.vcf.gz
done
```
