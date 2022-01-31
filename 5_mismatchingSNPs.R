# 5 mismatching SNPs
# 1. exclude mismatching SNPs (strand flipping)
.libPaths(c('~/R_libs', .libPaths()))
library(magrittr)
library(data.table)
# read the bimfile
filename= "/imaging/projects/cbu/calm-genomics/diandra_wd/tutorial/EUR.bim"

bim = fread (input= filename) %>%
  # replace the original column names by the new names
  setnames (.,colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1","B.A2")) %>%
  #change alleles to upper cases 
  .[,c("B.A1", "B.A2"):=list(toupper(B.A1), toupper(B.A2))]

# read in the summary statistics data 
height = fread("Height.QC.gz") %>%
  # change the alleles to upper cases
  .[,c("A1", "A2"):=list(toupper(A1), toupper(A2))]

# read in the QCed SNPS
filename2= "/imaging/projects/cbu/calm-genomics/diandra_wd/tutorial/EUR.QC.snplist"
qc = fread(input=filename2, header=F)

# 2. identify SNPs that require strand flipping
# merge summary stats with target
info = merge(bim, height, by=c("SNP", "CHR", "BP")) %>%
# filter QCed SNPs
.[SNP %in% qc[,V1]]
# function to calculate the complementary allele
complement <- function(x){
  switch(x, 
        "A"="T",
        "C"="G",
        "T"="A",
        "G"="C", 
        return(NA)
         )
}

# get SNPs that have the same alleles across base and target
info.match=info[A1==B.A1 & A2 == B.A2, SNP]
#identify snps that are complementary btw base & target
com.snps = info [sapply(B.A1, complement) == A1 &
                   sapply(B.A2, complement) == A2, SNP]

#info$C.A1 = sapply(info$B.A1, complement)
#info$C.A2 = sapply (info$B.A2, complement)
#info.complement = subset (info, A1==C.A1 & A2 == C.A2)

#update the bim file
bim[SNP %in% com.snps, c("B.A1", "B.A2"):=
      list(sapply(B.A1, complement),
           sapply(B.A2, complement))]

# 3. Identify SNPs that need recoding
recode.snps = info [B.A1==A2 & B.A2==A1, SNP]
# update the bim file
bim[SNP %in% recode.snps, c("B.A1", "B.A2"):=
      list(B.A2, B.A1)]
#identify SNPs that need recoding & complement
com.recode = info[sapply(B.A1, complement)== A2 &
                    sapply(B.A2, complement)== A1, SNP]
# now update the bim file
bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
      list(sapply(B.A2, complement), 
           sapply(B.A1, complement))]
# write the updated bim file
fwrite(bim[,c("SNP", "B.A1")], "EUR.a1", col.names=F, sep="\t")

# 4. identify SNPs that have different alleles in base and target
mismatch = bim[!(SNP %in% info.match |
                 SNP %in% com.snps  |
                 SNP %in% recode.snps |
                   SNP %in% com.recode), SNP]
write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)

# now you have to do the sex chromosome check

q()
