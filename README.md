# impute_tool
## introduction
Impute_tool is an open source software that implements imputing eQTL information across the genome. The software borrows the method of using the known eQTL statistics to predict unknown eQTL statistics in the ImpG-master software by using the linkage disequilibrium information between SNPs, and completes the imputing of eQTL information. ImpG-master implements the function of completing eQTL imputing in a given window of the user.Based on this, the software reduces the complexity of user input and automatically supplements the information needed to impute the location using 1000G files. In terms of performance, the method of turning on multithreading and preserving LD information is used.
Speed up the impution.
