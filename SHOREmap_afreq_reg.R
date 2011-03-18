# Read in command line options
args <- commandArgs()
print(args[5])
chrsize<-read.table(args[5])

pdf(file=args[6], width=20, height=5)
layoutmat=matrix(data=c(1), ncol=1, nrow=1)
layout(layoutmat)

options(scipen=999999999)

data<-read.table(args[7])
winstep=args[8]
winsize=args[9]
regchr=args[10]
regbegin=as.numeric(args[11])
regend=as.numeric(args[12])
regfreq_min=as.numeric(args[13])
regfreq_max=as.numeric(args[14])
#path=args[15]


#source(paste(path,"SHOREmap_confInt.R",sep="/"))

print (regfreq_min)
print (regfreq_max)

ls=round((regend-regbegin+1)/10)
max=1


# Plot

for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]

	if (chrname == regchr) {

		plot(data$V2[data$V1[]==chrname], data$V5[data$V1[]==chrname], ylim=c(regfreq_min, regfreq_max), xlim=c(regbegin, regend), type="p", pch=20, axes=F, xlab=paste("Chromosome ", chrname, sep=""), ylab="", main=paste("winstep:", winstep, " winsize:", winsize, sep=""))

		#ciData<- data[data[,1]==chrname&data[,2]>=regbegin&data[,2]<=regend,]
		#chromosome<-ciData[,1]
		#positions<-ciData[,2]
		#background_count<-ciData[,3??]
		#forground_count<-ciData[,4??]
		#error_count<-ciData[,6??]
		#ci<-ShoreMap.confint(chromosome,positions,background_count,foreground_count, error_count,foreground_frequency=1,level=0.99,minWindow=10)
		#plot??
		
		
		for (bgl in seq(0.1, 1, 0.1)) {
			abline(h=bgl, col="lightgrey")
		}

		labels=c(regbegin, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), regend)
		axis(1, label=labels, at=labels)
		axis(2, las=1, labels=c("0", "1"), at=c(0, max))

	}

}


# Thank you for riding:
dev.off()

