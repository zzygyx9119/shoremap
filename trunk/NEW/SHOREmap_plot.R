##########################################
# Read in command line options

args <- commandArgs()

target<-args[5] # todo set type
chrsize<-read.table(args[6])
fpdf<-args[7]
zoomf<-args[8] # todo check for existance
#windowsizes<-read.table(args[9], as.is=c(1))
winsize<-as.numeric(args[9])
winstep<-as.numeric(args[10])
path<-args[11]
outputpath<-args[12]
filterOutliers<-as.numeric(args[13]) # size of window around marker for outlier assessment
filterPValue<-as.numeric(args[14]) # p-value for outlier removal
conf_level<-as.numeric(args[15])
misscored<-as.numeric(args[16])
minMarker<-as.numeric(args[17])
minCov<-as.numeric(args[18])
maxCov<-as.numeric(args[19])
rMax<-as.numeric(args[20])
plotR<-as.numeric(args[21])
boostMax<-as.numeric(args[22])
plotBoost<-as.numeric(args[23])
peakwinsize<-as.numeric(args[24])
peakwinstep<-as.numeric(args[25])
runid<-as.numeric(args[26])

##########################################
# Load libraries

if (conf_level != 2) {
        library(bbmle)
}

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

pdf(file=fpdf, width=17, height=10) #*length(windowsizes$V1))
#png(filename=fpdf, width=17, height=10, units="in", res=72) 
layoutmat=matrix(data=c(1,2), ncol=1, nrow=2) #*length(windowsizes$V1))
layout(layoutmat)

options(scipen=999999999)
ls_default=5000000
max=1

##########################################
# Source confidence interval stats

source(paste(path,"SHOREmap_confInt.R", sep="/"))

##########################################
# Plotting

data<-read.table(paste(outputpath, "/SHOREmap.winsize1.txt", sep=""), as.is=c(1))

realRmax<-0
realBoostmax<-0

chrres = list()

for (chr in 1:(length(chrsize$V1))) {
	chrname=chrsize$V1[chr]
        ciData<- data[data[,1]==chrname,]
	ci_chromosome<-ciData[,1]
       	ci_positions<-ciData[,2]
        ci_background_count<-ciData[,4]
    	ci_forground_count<-ciData[,3]
        ci_error_count<-ciData[,5]
      	#ci_result<-ShoreMap.confint(ci_chromosome, ci_positions, ci_background_count, ci_forground_count, ci_error_count, foreground_frequency=target, level=2, recurse=F, forceInclude=T, allowAdjustment=misscored, filterOutliers=0, filterPValue=filterPValue,winSize=winsize,winStep=winstep,minMarker=minMarker,minCoverage=minCov,peakFinding=3)
	ci_result<-ShoreMap.confint(ci_chromosome, ci_positions, ci_background_count, ci_forground_count, ci_error_count, foreground_frequency=target, level=conf_level, recurse=F, forceInclude=T, allowAdjustment=misscored, filterOutliers=filterOutliers, filterPValue=filterPValue, winSize=winsize,winStep=winstep, minMarker=minMarker, minCoverage=minCov,peakFinding=3, peakWinSize=peakwinsize, peakWinStep=peakwinstep)
	realBoostmax<-max(ci_result$averaged[,3], realBoostmax)
	#realRmax<-max(ci_result$averaged[,4], realRmax, peakWinSize=peakwinsize, peakWinStep=peakwinstep)
	realRmax<-max(ci_result$averaged[,4], realRmax)

	chrres[[chr]] = ci_result
}


for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]
	ls = min(ls_default, chrsize$V2[chr])

	if (z_chr == -1 || chrname == z_chr) {

		if (length(data$V2[data$V1[]==chrname])!= 0) {

			# Set up plot data
			freq = data$V3[data$V1[]==chrname] / ( data$V3[data$V1[]==chrname] + data$V4[data$V1[]==chrname] )

			## Calc confidence interval
       		        #ciData<- data[data[,1]==chrname,]
	       		#ci_chromosome<-ciData[,1]
       	        	#ci_positions<-ciData[,2]
       		        #ci_background_count<-ciData[,4]
       	        	#ci_forground_count<-ciData[,3]
        	        #ci_error_count<-ciData[,5]
	        	#ci_result<-ShoreMap.confint(ci_chromosome, ci_positions, ci_background_count, ci_forground_count, ci_error_count, foreground_frequency=target, level=conf_level, recurse=F, forceInclude=T, allowAdjustment=misscored, filterOutliers=filterOutliers, filterPValue=filterPValue, winSize=winsize,winStep=winstep, minMarker=minMarker, minCoverage=minCov,peakFinding=3, peakWinSize=peakwinsize, peakWinStep=peakwinstep)

			ci_result = chrres[[chr]]

       		        ciData<- data[data[,1]==chrname,]
       	        	ci_positions<-ciData[,2]

                	ci<-ci_result$confidenceInterval
       	               	ci_filtered <- ci_positions %in% ci_result$excluded 

       	               	#for plotting the windowed markers
			ci_avgPosFreq <- ci_result$averaged

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


			plot(data$V2[data$V1[]==chrname], freq, ylim=c(y_min, y_max+0.2), xlim=c(x_min, x_max), type="n", axes=F, xlab="", ylab="Allele Frequency", main=paste("Chromosome:", chrname, " (In black, window size of ", winsize, " bp.)", sep=""))
			
			# plot grid
			for (bgl in seq(0.1, 1, 0.1)) {
       		               	lines(c(1, chrsize$V2[chr]), c(bgl, bgl), col="lightgrey", lty=2)
               		}
	

			####################################################		
                        # boost
                        if (plotBoost == 1) {
				print("Plotting boost...")
                                lines(ci_avgPosFreq[,1], ci_avgPosFreq[,3]/realBoostmax,col="#8B8970CC")
                                #polygon(ci_avgPosFreq[,1], ci_avgPosFreq[,3]/realBoostmax, col="lemonchiffon4", border = NA)
                        }

                        #r
                        if (plotR == 1) {
				print("Plotting r...")
                                lines(ci_avgPosFreq[,1], ci_avgPosFreq[,4]/realRmax, col="#8B8970CC")
                                #polygon(ci_avgPosFreq[,1], ci_avgPosFreq[,4]/realRmax, col="#8B897033", border="#8B8970CC")
                        }

			####################################################
			# Plot AFE
			points(data$V2[data$V1[]==chrname], freq, ylim=c(y_min, y_max+0.2), col=ifelse(conf_level!=2 & ci_filtered, "purple2", "lightblue"), xlim=c(x_min, x_max), pch=ifelse(conf_level!=2 & ci_filtered, 4, 16), cex=0.75)


			# Plot confidence interval
			if (conf_level != 2 && ci[3, 1] <= 1) {
       				for (ci_i in 1:(length(ci[1,]))) {
				
					rect(ci[1, ci_i], y_max+0.02, ci[2, ci_i], y_max+0.06, col="maroon2", border="maroon2")
					
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
			#abline(v=16702262, col="limegreen")
			# Vini:
			#abline(v=18816001, col="limegreen")
			# Qtl:
			if(file.exists(paste("Sim.",runid,".qtl.txt", sep=""))) {
				qtl<-read.table(paste("Sim.",runid,".qtl.txt", sep=""), header=F)
				for (q in 1:(length(qtl$V3))) {
					if (qtl$V3[q] == chrname) {
						abline(v=qtl$V4[q])
						#text(c(qtl$V4[q]), c((y_min+(y_max-y_min)/2)), labels=c(paste(round(qtl$V5), round(qtl$V6), sep=" ")))
						text(c(qtl$V4[q]), c((y_min+(y_max-y_min)/2)), labels=c(paste(round(qtl$V5[q]), round(qtl$V6[q]), sep=" ")), col="green")
					}
				}
			}


			# Axes descriptions
			labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
			axis(1, label=labels, at=labels, col="lightgrey")
			axis(2, las=1, labels=c(paste(y_min,sep=""), paste(y_max, sep="")), at=c(y_min, y_max), col="lightgrey")


			# Sliding window
           		#points(ci_avgPosFreq[,1], ci_avgPosFreq[,2], ylim=c(y_min, y_max+0.2), col="grey2", xlim=c(x_min, x_max), pch=16, cex=0.75)
			dist=diff(ci_avgPosFreq[,1])
			bl=x_max*0.005
			indices= c(0,which(dist>bl),length(ci_avgPosFreq[,1]))
			for(i in 2:length(indices)){
				lines(ci_avgPosFreq[(indices[i-1]+1):(indices[i]),1], ci_avgPosFreq[(indices[i-1]+1):(indices[i]),2], ylim=c(y_min, y_max+0.2), col="grey2", xlim=c(x_min, x_max), cex=0.75, lwd=2)
			}


			#######################################################

			# Plot support
	
	               	#winsize=windowsizes$V1[winsize_i]
		        #data<-read.table(paste(outputpath, "/SHOREmap.winsize", winsize, ".txt", sep=""), as.is=c(1))


			max_read = max(max(data$V3[data$V1[]==chrname & data$V3[]>=minCov & data$V3[]<=maxCov]), max(data$V4[data$V1[]==chrname & data$V3[]>=minCov & data$V3[]<=maxCov]))
	
			plot(data$V2[data$V1[]==chrname & data$V3[]>=minCov & data$V3[]<=maxCov], data$V3[data$V1[]==chrname & data$V3[]>=minCov & data$V3[]<=maxCov], type="h", col="steelblue4", xlim=c(x_min, x_max), ylim=c((-1)*max_read, max_read), axes=F, ylab="Read count", xlab="")
			points(data$V2[data$V1[]==chrname & data$V4[]>=minCov & data$V4[]<=maxCov], (-1)*data$V4[data$V1[]==chrname & data$V4[]>=minCov & data$V4[]<=maxCov], type="h", col="olivedrab");
	
			axis(1, label=labels, at=labels)
			axis(2, las=1, labels=c(paste(max_read), "0", paste(max_read,sep="")), at=c((-1)*max_read, 0, max_read))
		}
	}
}


# close graphics device
dev.off()

