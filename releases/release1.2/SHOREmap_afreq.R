# Read in command line options
args <- commandArgs()
chrsize<-read.table(args[5])

pdf(file=args[6], width=10, height=20)
layoutmat=matrix(data=c(1,2,3,4,5), ncol=1, nrow=5)
layout(layoutmat)

options(scipen=999999999)
ls=5000000

data<-read.table(args[7])
max=1




winstep=args[8]
winsize=args[9]


# Plot

for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]

	plot(data$V2[data$V1[]==chrname], data$V5[data$V1[]==chrname], ylim=c(0, 1), xlim=c(0, max(chrsize$V2)), type="p", pch=20, axes=F, xlab=paste("Chromosome ", chrname, sep=""), ylab="", main=paste("winstep:", winstep, " winsize:", winsize, sep=""))

	for (bgl in seq(0.1, 1, 0.1)) {
		abline(h=bgl, col="lightgrey")
	}

	labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
	axis(1, label=labels, at=labels)
	axis(2, las=1, labels=c("0", "1"), at=c(0, max))

	#abline(v=16240000)

}


# Thank you for riding:
dev.off()

