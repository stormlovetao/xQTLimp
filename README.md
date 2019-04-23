# xQTLImp
## introduction
 xQTLImp is an open source software that implements imputing eQTL information across the genome. The software borrows the method of using the known eQTL statistics to predict unknown eQTL statistics in the ImpG-master software by using the linkage disequilibrium information between SNPs, and completes the imputing of eQTL information. ImpG-master implements the function of completing eQTL imputing in a given window of the user.Based on this, the software reduces the complexity of user input and automatically supplements the information needed to impute the location using 1000G files. In terms of performance, the method of turning on multithreading and preserving LD information is used.
Speed up the impution.
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
eg:./xQTLImp -g /media/userdisk1/jjpeng/yinquanwei/gencode_v19_gene_annotation.txt -e /media/userdisk1/jjpeng/yinquanwei/Brain_Amygdala.allpairs.txt -v /media/userdisk1/jjpeng/yinquanwei/ -o /media/userdisk1/jjpeng/yinquanwei/output3/ -t 55

###demon
There is a demon in sample folder.

