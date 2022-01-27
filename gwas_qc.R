# this script performs basic QC of the base data together with the awk_scripts
#load the packages
install.packages('R.utils')
library(R.utils)
require(data.table)
library(data.table)

# read in the file Height.txt (your trait of interest)
filename= "/imaging/projects/cbu/calm-genomics/diandra_wd/tutorial/Height.gwas.txt.gz" #put the correct file to your gwas summary stats here
dat = fread(input=filename)

# filter out SNPs
result = dat[INFO > 0.8 & MAF >0.01]

# output the gz file
fwrite(result, "Height.gz", sep="\t")
