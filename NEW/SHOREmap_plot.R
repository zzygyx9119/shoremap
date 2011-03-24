library(bbmle)

##########################################
# Read in command line options
args <- commandArgs()

target<-args[5] # todo set type
chrsize<-read.table(args[6])
fpdf<-args[7]
#zoom<-read.table(args[8]) # todo check for existance
windowsizes<-read.table(args[9])
winstep<-as.numeric(args[10])
path<-args[11]
outputpath<-args[12]

##########################################
# Set up graphical device
# and graphics related parameter

pdf(file=fpdf, width=15, height=8*length(windowsizes$V1))
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

for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]

	for (winsize_i in 1:(length(windowsizes$V1))) {

		winsize=windowsizes$V1[winsize_i]
		data<-read.table(paste(outputpath, "/SHOREmap.winsize", winsize, ".txt", sep=""))

		# Set up plot data
		freq = data$V3[data$V1[]==chrname] / ( data$V3[data$V1[]==chrname] + data$V4[data$V1[]==chrname] )
		temp_winstep = winstep
		if (winsize == 1) {
			temp_winstep = 1
		}

		# Calc confidence interval
                ciData<- data[data[,1]==chrname,]
                ci_chromosome<-ciData[,1]
                ci_positions<-ciData[,2]
                ci_background_count<-ciData[,4]
                ci_forground_count<-ciData[,3]
                ci_error_count<-ciData[,5]
                ci<-ShoreMap.confint(ci_chromosome, ci_positions, ci_background_count, ci_forground_count, ci_error_count, foreground_frequency=target, level=0.999, recurse=T)

		# Plot
		plot(data$V2[data$V1[]==chrname], freq, ylim=c(0, 1.2), xlim=c(0, max(chrsize$V2)), type="p", pch=20, axes=F, ylab="", main=paste("Winstep:", temp_winstep, " Winsize:", winsize, sep=""))
		#plot(data$V2[data$V1[]==chrname], freq, ylim=c(0, 1), xlim=c(0, max(chrsize$V2)), type="p", pch=20, axes=F, xlab=paste("Chromosome ", chrname, sep=""), ylab="", main=paste("winstep:", winstep, " winsize:", winsize, sep=""))
		#plot confidence interval
		if (ci[3, 1] <= 1) {
                	for (ci_i in 1:(length(ci[1,]))) {
				rect(ci[1, ci_i], 1.01, ci[2, ci_i], 1.1)
				# add positions
                	}
		}
		else {
			# Indicate that there is no interval
		}



		# Making it pretty
		# TODO: put to background
		for (bgl in seq(0.1, 1, 0.1)) {
			abline(h=bgl, col="lightgrey")
		}

		labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
		axis(1, label=labels, at=labels)
		axis(2, las=1, labels=c("0", "1"), at=c(0, max))

		#abline(v=16240000)

	}
}


# close graphics device
dev.off()

