# Read in command line options
args <- commandArgs()
chrsize<-read.table(args[5])

pdf(file=args[6], width=10, height=12)
layoutmat=matrix(data=c(1,2,3,4,5), ncol=1, nrow=5)
layout(layoutmat)

options(scipen=999999999)
lsby=1000000

data<-read.table(args[7])

max=max(data$V5)
min=min(data$V5)
extr=max(max, (-1)*min)

winstep=args[8]
winsize=args[9]


# Plot

for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]

	plot(data$V2[data$V1[]==chrname], data$V5[data$V1[]==chrname], ylim=c((-1)*extr, extr), xlim=c(0, max(chrsize$V2)), type="l", axes=F, xlab=paste("Chromosome ", chrname, " [Mb]", sep=""), ylab=paste(args[10], args[11], sep="    "), main=paste("winstep:", winstep, " winsize:", winsize, sep=""))

	print(lsby)
	print(chrsize$V2[chrsize$V1[]==chrname])
	labelsat=c(0, seq(lsby, chrsize$V2[chrsize$V1[]==chrname], by=lsby), chrsize$V2[chrsize$V1[]==chrname])
	labels=labelsat/lsby
	axis(1, label=labels, at=labelsat)
	axis(2, las=1, labels=c("-1", "0", "1"), at=c((-1)*extr, 0, extr))

	#abline(h=0)

}


# Thank you for riding:
dev.off()

