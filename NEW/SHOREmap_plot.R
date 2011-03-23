##########################################
# Read in command line options
args <- commandArgs()

chrsize<-read.table(args[5])
fpdf<-read.table(args[6])
zoom<-read.table(args[7])
windowsizes<-read.table(args[8])
winstep<-as.numeric(args[9])
path<-args[10]

##########################################
# Set up graphical device
# and graphics related parameter

pdf(file=fpdf, width=15, height=8*length(sizes))
layoutdat=seq(1, 2*length(windowsizes$V1))
layoutmat=matrix(data=layoutdat, ncol=1, nrow=2*length(windowsizes$V1))
layout(layoutmat)

options(scipen=999999999)
ls=5000000
max=1


##########################################
# Source confidence interval stats

source(paste(path,"SHOREmap_confInt.R", sep="/"))


##########################################
# Plotting

for (winsize_i in 1:(length(windowsizes$V1))) {

	winsize=windowsizes$V1[winsize_i]
	data<-read.table(paste(path, "/SHOREmap.winsize", winsize, ".txt", sep=""))

	for (chr in 1:(length(chrsize$V1))) {

		chrname=chrsize$V1[chr]

		plot(data$V2[data$V1[]==chrname], data$V5[data$V1[]==chrname], ylim=c(0, 1), xlim=c(0, max(chrsize$V2)), type="p", pch=20, axes=F, xlab=paste("Chromosome ", chrname, sep=""), ylab="", main=paste("winstep:", winstep, " winsize:", winsize, sep=""))

		ciData<- data[data[,1]==chrname,]
		chromosome<-ciData[,1]
		positions<-ciData[,2]
		background_count<-ciData[,3??]
		forground_count<-ciData[,4??]
		error_count<-ciData[,6??]
		ci<-ShoreMap.confint(chromosome,positions,background_count,foreground_count, error_count,foreground_frequency=1,level=0.99,minWindow=10)
		#plot??
	
		for (bgl in seq(0.1, 1, 0.1)) {
			abline(h=bgl, col="lightgrey")
		}

		labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
		axis(1, label=labels, at=labels)
		axis(2, las=1, labels=c("0", "1"), at=c(0, max))

		#abline(v=16240000)

}


# close graphics device
dev.off()

