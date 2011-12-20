#from Jonas Klasen

args<-commandArgs(trailingOnly=TRUE)

infile<-args[1]
outfile<-args[2]

### R qtl
library("qtl")

### Read data
data.normal <- read.cross("csv", ".", infile, sep="\t")
                          #, genotypes = c("L", "H", "C", "not LL", "not CC"), alleles=c("L", "C")

### Composid intervall mapping
data.normal <- calc.genoprob(data.normal)
#perm.cim.normal <- cim(data.normal, n.marcovar=4, n.perm=3, window=300, method="hk")
fit.cim.normal <- cim(data.normal, n.marcovar=4, window=300, method="hk")
#all.marker.cim.normal <- fit.cim.normal[fit.cim.normal[,3] >= mean(perm.cim.normal), ]
#plot(fit.cim.normal)
#abline(h=mean(perm.cim.normal))
write.table(fit.cim.normal,outfile,sep="\t",row.names=FALSE,col.names=FALSE)
