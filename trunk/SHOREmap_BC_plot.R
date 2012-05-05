arg <- commandArgs()
data <-read.table(arg[5])
R_out <-paste(arg[6],".pdf",sep = "")
R_out_summary <-paste(arg[6],"_summary.pdf",sep = "")
chrsize <-read.table(arg[7])
summary <-arg[8]
only_EMS <- arg[9]
other_mutant <-arg[10]
min_fre =0



color.bar <- function(mat, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
	range = (length(mat)-1)/(max-min)
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	axis(2, ticks, las=1)
	for (i in 1:(length(mat)-1)) {
		y = (i-1)/range + min
		rect(0,y,10,y+1/range, col=mat[i], border=NA)
	}
}

## setting the varibles for plotting ##
## plotting symboles are, EMS is pch=20, Non-EMS is pch=4(x) and if the study is other_mutant then pch="*"

col_martix<-c(rainbow(15, start=.17, end=.3),rainbow(9, start=.31, end=.4),rainbow(11, start=.7, end=.8),rainbow(5, start=.9, end=1))
ylabels=c(seq(min_fre, 1, by=0.2), 1)
xmax = max(chrsize$V2)



## plotting summary file ##
if(summary==1){
	pdf(file = R_out_summary, width=17, height=10)
	par(mfrow=c(3,2))
	for (chr in 1:(length(chrsize$V1))) {
        	chrname=chrsize$V1[chr]
		ls = min(5000000, chrsize$V2[chr])
		if(!(other_mutant==1)){
			EMS<-data[((data[,4]=="G" & data[,5]=="A") | (data[,4]=="C" & data[,5]=="T")) ,]
			EMS_no<-data[((data[,4] =="G" & data[,5]!="A") | (data[,4]=="C" & data[,5]!="T") | (data[,4]=="A") | data[,4]=="T" | data[,4]=="-") ,]
			plot(EMS$V3[EMS$V2==chrname],EMS$V8[EMS$V2==chrname], col=col_martix[EMS$V6[EMS$V2==chrname]],pch=20, main=paste("Chromosome",chrname) ,ylab="Frequency" ,  xlim=c(1,xmax), ylim=c(min_fre,1),axes=F, xlab="")
			if(!(only_EMS ==1)) {
				points(EMS_no$V3[EMS_no$V2==chrname],EMS_no$V8[EMS_no$V2==chrname], col=col_martix[EMS$V6[EMS$V2==chrname]],pch=4)
			}
			
		}
		else{
			plot(data$V3[data$V2==chrname],data$V8[data$V2==chrname], col=col_martix[data$V6[data$V2==chrname]],pch="*", main=paste("Chromosome",chrname) ,ylab="Frequency" ,  xlim=c(1,xmax), ylim=c(min_fre,1),axes=F, xlab="")
		}
		labels=c(1, seq(ls, chrsize$V2[chrname], by=ls), chrsize$V2[chrname])
		axis(1, label=labels, at=labels, col="lightgrey")
		axis(2, label=ylabels, at=ylabels, col="lightgrey")
	}
	color.bar(col_martix, 0, 40, title='Shore score', nticks=5)
	dev.off()
}


## plotting each chromosome ##
pdf(file = R_out,width=17, height=10)
for (chr in 1:(length(chrsize$V1))) {
        chrname=chrsize$V1[chr]
	ls = min(5000000, chrsize$V2[chr])
	if(!(other_mutant==1)){
		EMS<-data[((data[,4]=="G" & data[,5]=="A") | (data[,4]=="C" & data[,5]=="T")) ,]
		EMS_no<-data[((data[,4] =="G" & data[,5]!="A") | (data[,4]=="C" & data[,5]!="T") | (data[,4]=="A") | data[,4]=="T" | data[,4]=="-") ,]
		plot(EMS$V3[EMS$V2==chrname],EMS$V8[EMS$V2==chrname], col=col_martix[EMS$V6[EMS$V2==chrname]],pch=20, main=paste("Chromosome",chrname) ,ylab="Frequency" , xlim=c(1,xmax), ylim=c(min_fre,1),axes=F, xlab="")
		if(!(only_EMS ==1)) {
			points(EMS_no$V3[EMS_no$V2==chrname],EMS_no$V8[EMS_no$V2==chrname], col=col_martix[EMS$V6[EMS$V2==chrname]],pch=4)
		}
	}
	else{
		plot(data$V3[data$V2==chrname],data$V8[data$V2==chrname], col=col_martix[data$V6[data$V2==chrname]],pch="*", main=paste("Chromosome",chrname) ,ylab="Frequency" , xlim=c(1,xmax), ylim=c(min_fre,1),axes=F, xlab="")
	}
	labels=c(1, seq(ls, chrsize$V2[chrname], by=ls), chrsize$V2[chrname])
	axis(1, label=labels, at=labels, col="lightgrey")
	axis(2, label=ylabels, at=ylabels, col="lightgrey")
}
color.bar(col_martix, 0, 40, title='Shore score', nticks=5)
dev.off()

