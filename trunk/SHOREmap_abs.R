# Read in command line options
args <- commandArgs()
chrsize<-read.table(args[5])

pdf(file=args[6], width=10, height=20)
layoutmat=matrix(data=c(1,2,3,4,5), ncol=1, nrow=5)
layout(layoutmat)

options(scipen=999999999)
ls=5000000

data<-read.table(args[7])
max=max(max(data$V5), max(data$V6))

source(SHOREmap_confInt.R)

winstep=args[8]
winsize=args[9]


# Plot

for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]

	plot(data$V2[data$V1[]==chrname], data$V5[data$V1[]==chrname], ylim=c((-1)*max, max), xlim=c(0, max(chrsize$V2)), col="red", type="h", axes=F, xlab=paste("Chromosome ", chrname, sep=""), ylab="", main=paste("winstep:", winstep, " winsize:", winsize, sep=""))

	#ciData<- data[data[,1]==chrname,]
	#chromosome<-ciData[,1]
	#positions<-ciData[,2]
	#background_count<-ciData[,3??]
	#forground_count<-ciData[,4??]
	#error_count<-ciData[,6??]
	#ci<-ShoreMap.confint(chromosome,positions,background_count,foreground_count, error_count,foreground_frequency=1,level=0.99,minWindow=10)
	#plot??

	x=data$V2[data$V1[]==chrname]
	y=data$V6[data$V1[]==chrname]
	for (i in (1:length(x))) {
		lines(c(x[i],x[i]), c(0,(-1)*y[i]), col="blue")
	}

	labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
	axis(1, label=labels, at=labels)
	axis(2, las=1, labels=c(-max, "0", max), at=c(-max, 0, max))

	#abline(v=16240000)

}


# Thank you for riding:
dev.off()

