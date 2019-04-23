# xQTLImp
## introduction
 xQTLImp is an open source software that implements imputing xQTL(eQTL,mQTL,haQTL) information across the genome. xQTLImp borrows the functions of ImpG-master in its implementation. ImpG-master implements impute function from a window given by user, and needs to provide files containing all snps information and typed snps zscore in the window . Based on this, the software reduces the complexity of user input and automatically supplements the information needed to impute the location using 1000G files. In terms of performance, the method of turning on multithreading and preserving LD information is used to speed up the impution. xQTLImp is suitable for handling more general cases. The tool extends the type of processing to Snps type and Indel type and supports x molecular features, such as gene/transcript/exon expression level, DNA methylation, histone acylation…
 </br>
##  Building xQTLImp
 To read gz format, you need zlib.
```bash
sudo apt-get install zlib1g-dev
```
Any C++11 compiler should work.
```bash
make #under src folder
```
## Usage
### Requirements for input files
There input files are required:
#### 1.1000G files in vcf format 
[Available download address](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/)</br>
The 1000G files should be compressed into gz format</br>
The 1000G files should be named like chrom[1-22].vcf.gz in a folder</br>
</br>
#### 2.Molecular trait file
Molecular trait file must start with a line that contains column labels – molecular_ID, start_pos, end_pos ...(option) </br>followed by lines of data entries. Each field of data entries must be separated by white spaces.</br>
##### Example:
molecular_ID start_pos end_pos</br>
ENSG00000223972.4	11869	14412</br>
ENSG00000227232.4	14363	29806</br>
ENSG00000243485.2	29554	31109</br>
ENSG00000237613.2	34554	36081</br>
</br>
#### 3.xQTL file
xQTL file must start with a line that contains column labels – chromosome , molecular_ID, variant_pos , Ref_allele , Alt_allele , z_statistics ...(option)</br> followed by lines of data entries. Each field of data entries must be separated by white spaces.</br>
##### Example:
chromosome molecular_ID variant_pos Ref_allele  Alt_allele z_statistics</br>
1 ENSG00000223972.4 13417 G C 1.5</br>
1 ENSG00000227232.4 17559 A G 2.6</br>
1 ENSG00000227232.4 54421 G A -1.0</br>
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
22 subfolders will appear in the output folder,and each folder corresponds to the imputing result on the chromosome.
