##########################################
# Load libraries

library(bbmle)

##########################################
# Read in command line options

args <- commandArgs()

target<-args[5] # todo set type
chrsize<-read.table(args[6])
fpdf<-args[7]
zoomf<-args[8] # todo check for existance
windowsizes<-read.table(args[9])
winstep<-as.numeric(args[10])
path<-args[11]
outputpath<-args[12]
filterOutliers<-args[13] # size of window around marker for outlier assessment
filterPValue<-args[14] # p-value for outlier removal

##########################################
# Set up zoom interval if it exists

z_chr = -1

if(file.exists(zoomf)) {
	zoom=read.table(zoomf)
	z_chr = zoom$V1
	z_beg = zoom$V2
	z_end = zoom$V3
	z_min = zoom$V4
	z_max = zoom$V5
}

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

	if (z_chr == -1 || chrname == z_chr) {

		for (winsize_i in 1:(length(windowsizes$V1))) {
	
			winsize=windowsizes$V1[winsize_i]
			data<-read.table(paste(outputpath, "/SHOREmap.winsize", winsize, ".txt", sep=""))

#			print(paste("Analysing: ",outputpath, "/SHOREmap.winsize", winsize, ".txt", sep=""))

			# Set up plot data
			freq = data$V3[data$V1[]==chrname] / ( data$V3[data$V1[]==chrname] + data$V4[data$V1[]==chrname] )

			# Calc confidence interval
        	        ciData<- data[data[,1]==chrname,]
                	ci_chromosome<-ciData[,1]
	                ci_positions<-ciData[,2]
        	        ci_background_count<-ciData[,4]
                	ci_forground_count<-ciData[,3]
	                ci_error_count<-ciData[,5]
        	        ci_result<-ShoreMap.confint(ci_chromosome, ci_positions, ci_background_count, ci_forground_count, ci_error_count, foreground_frequency=target, level=0.999, recurse=F,forceInclude=T,allowAdjustment=0.05,filterOutliers=filterOutliers, filterPValue=filterPValue)

        	        ci<-ci_result$confidenceInterval
                        ci_filtered <- ci_positions %in% ci_result$excluded

        	        #extract removed positions
                        

			########################################################

			# Plot frequencies
			y_min = 0
			y_max = 1.0
			x_min = 1
			x_max = max(chrsize$V2)

			if (z_chr != -1) {
				y_min = z_min
				y_max = z_max
				x_min = z_beg
				x_max = z_end	
			}

			plot(data$V2[data$V1[]==chrname], freq, ylim=c(y_min, y_max+0.2), xlim=c(x_min, x_max), type="n", axes=F, xlab="", ylab="Allele Frequency", main=paste("Chromosome:", chrname, " (Using window size of ", winsize, " reporting every ", winstep, " bp.)", sep=""))

			for (bgl in seq(0.1, 1, 0.1)) {
                                abline(h=bgl, col="lightgrey")
                        }

			points(data$V2[data$V1[]==chrname], freq, ylim=c(y_min, y_max+0.2), col=ifelse(ci_filtered,"red","black"), xlim=c(x_min, x_max), pch=20)
		
			# Plot confidence interval
			if (ci[3, 1] <= 1) {
                		for (ci_i in 1:(length(ci[1,]))) {
			
					rect(ci[1, ci_i], y_max+0.02, ci[2, ci_i], y_max+0.06, col="orange", border="orange")
					#rect(ci[1, ci_i], y_max+0.02, ci[2, ci_i], y_max+0.06, col="orange", border="orange")
					#rect(ci[1, ci_i], y_max+0.02, ci[2, ci_i], y_max+0.06, col="orange", border="orange")
					
					# add size
					#size = ci[2, ci_i] - ci[1, ci_i] + 1
					#text(c(ci[1, ci_i] + size/2), c(y_max+0.035), labels=c(paste(size, sep="")))
					# add positions
					text(c(ci[1, ci_i]), c(y_max+0.09), labels=c(paste(ci[1, ci_i], sep="")), pos=2)
					text(c(ci[2, ci_i]), c(y_max+0.09), labels=c(paste(ci[2, ci_i], sep="")), pos=4)
        	        	}
			}
			else {
				# Indicate that there is no interval
			}

			# debug:
			# Stig:
			abline(v=16702262, col="limegreen")

			# Vini:
			#abline(v=18816001, col="limegreen")

			labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
			axis(1, label=labels, at=labels)
			axis(2, las=1, labels=c(paste(y_min,sep=""), paste(y_max, sep="")), at=c(y_min, y_max))


			#######################################################

			# Plot support
	
			max_read = max(max(data$V3[data$V1[]==chrname]), max(data$V4[data$V1[]==chrname]))

			plot(data$V2[data$V1[]==chrname], data$V3[data$V1[]==chrname], type="h", col="darkblue", xlim=c(x_min, x_max), ylim=c((-1)*max_read, max_read), axes=F, ylab="Read count", xlab="")
			points(data$V2[data$V1[]==chrname], (-1)*data$V4[data$V1[]==chrname], type="h", col="darkred");

			axis(1, label=labels, at=labels)
			axis(2, las=1, labels=c(paste((-1)*max_read), "0", paste(max_read,sep="")), at=c((-1)*max_read, 0, max_read))

		}
	}
}


# close graphics device
dev.off()

