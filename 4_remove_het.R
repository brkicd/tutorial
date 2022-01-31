# remove heterozygosity from the target sample

library(data.table)
#read in file 
filename= "/imaging/projects/cbu/calm-genomics/diandra_wd/tutorial/EUR.QC.het" #put the correct path to your file here
dat = fread(input=filename)
# get samples wih F coefficient within 3SD of the poulation mean
valid = dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)]
# print FID and IID for valid samples
fwrite(valid[,c("FID", "IID")], "EUR.valid.sample", sep="\t")
