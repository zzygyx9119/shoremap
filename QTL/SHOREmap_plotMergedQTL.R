#grep -v '#' */*v10.qtlEstimates.csv | sed s/^.*Sim/Sim/|sed s/_v...qtlEstimates.csv:/"\t"/|awk 'BEGIN{OFS="\t"} {split($1,l,"."); print "Sim."l[2]"\t"l[3]"\t"$0}' |cat ../../../header.txt - > merged_v10_all.csv

#/usr/bin/R --vanilla --args merged_v10_all.csv merged_v10_all < ~/shoreMap/NEW/SHOREmap_plotMergedQTL.R

library(reshape2)
args<- commandArgs(trailingOnly=TRUE)
data<-read.table(args[1],sep="\t",fill=TRUE,header=TRUE)

outprefix<-args[2]

coverages<-sort(unique(data[,2]))

covCount<-table(unique(data[,1:2])[,1])
data<-data[!(data[,1] %in% names(covCount)[covCount<length(coverages)]),]



#remove uncomplete simulations




nrOfsimulations<-length(unique(data[,1]))
if(nrOfsimulations>2500){
 data<-data[data[,1] %in% matrix(unique(data[,1])[1:2500]),]
}
nrOfsimulations<-length(unique(data[,1]))


nrOfQtls<-sum(!is.na(unique(data$id)))
all_colors<-c("#cf5f3b", "#e78166", "#ffa391", "#f0c17f", "#e0de6d", "#bbee9e", "#95fed0","#51f0e8", "#0ce3ff", "#06acf1", "#0075e2")
colors<-all_colors[1:nrOfQtls]

#sensitivity
tp<-data[data$judgement=="TP",]
tp_counts<-tapply(tp[,2],tp[,8],table)

pdf(paste(outprefix,".sensitivity.pdf",sep=""))
plot(c(),xlim=c(1,length(coverages)),ylim=c(0,100),ylab="%",main="Sensitivity",xlab="coverage",xaxt="n")
axis(1,at=1:length(coverages),labels=coverages)
for(i in 1:length(tp_counts)){
 lines(tp_counts[[i]]/nrOfsimulations*100,col=colors[i],lwd=2,type="b")
}
legend(1, 20, legend=1:nrOfQtls, col=colors, ncol=5, lty=1, lwd=2, title="Rank")#, horiz=TRUE)

dev.off()

#specificity
tp_dist<-tapply(1:nrow(tp),tp[,2],function(cov) tapply(abs(tp[cov,17]-tp[cov,6]),tp[cov,8],function(x) x))

pdf(paste(outprefix,".specificity.full.pdf",sep=""))
top<-max(c(tp_dist,recursive=TRUE))
bp<-boxplot(do.call(c,tp_dist), col=rep(colors, length(coverages)), main="Specificity", xaxt="n", xlab="coverage groups", ylab="distance to true position",ylim=c(-top/5,1.05*top))
axis(1,at=seq(nrOfQtls/2,nrOfQtls*length(coverages),nrOfQtls),labels=coverages,las=2)
legend(nrOfQtls*length(coverages)*0.25,-top/70,legend=1:nrOfQtls,col=colors,lty=1,lwd=10,title="Rank",ncol=5)
dev.off()


pdf(paste(outprefix,".specificity.noOutliers.pdf",sep=""))
top<-max(bp$stats)
boxplot(do.call(c,tp_dist), col=rep(colors, length(coverages)), main="Specificity", xaxt="n", xlab="coverage groups", ylab="distance to true position",outline=FALSE,ylim=c(-top/5,1.05*top))
axis(1,at=seq(nrOfQtls/2,nrOfQtls*length(coverages),nrOfQtls),labels=coverages,las=2)
legend(nrOfQtls*length(coverages)*0.25,-top/70,legend=1:nrOfQtls,col=colors,lty=1,lwd=10,title="Rank",ncol=5)
dev.off()

#correl
freqDiff_dist<-tapply(1:nrow(tp),tp[,2],function(cov) tapply(tp[cov,16],tp[cov,8],function(x) x))

pdf(paste(outprefix,".correlation.pdf",sep=""))
boxplot(do.call(c,freqDiff_dist),col=rep(colors,length(coverages)),ylim=c(0,1),main="Correlation",xaxt="n",xlab="coverage groups",ylab="frequency difference")
axis(1,at=seq(nrOfQtls/2,nrOfQtls*length(coverages),nrOfQtls),labels=coverages,las=2)
legend(nrOfQtls*length(coverages)*0.25,1,legend=1:nrOfQtls,col=colors,lty=1,lwd=10,title="Rank",ncol=5)
dev.off()

freqDiff_median<-tapply(1:nrow(tp),tp[,2],function(cov) tapply(tp[cov,16],tp[cov,8],median))


names=c(sapply(sort(unique(tp[,2])),function(i) c("","",i,"","")))

#false positives
fp<-data[data$judgement=="FP"& data$spec=="close",]
fp_cat_count<-sapply(1:length(coverages),function(i) table(sapply(fp[fp[,2]==coverages[i],16],function(x) which.min(abs(freqDiff_median[[i]]-x)) )))

pdf(paste(outprefix,".false.positives.close.pdf",sep=""))
points<-barplot(fp_cat_count/nrOfsimulations,col=colors,legend=1:nrOfQtls,main="False positives (interval border within 2Mb)",xlab="coverages",ylab="average number of False Positives")
axis(1,at=points,labels=coverages,las=2)
dev.off()

fp<-data[data$judgement=="FP"& data$spec!="close",]
fp_cat_count<-sapply(1:length(coverages),function(i) table(sapply(fp[fp[,2]==coverages[i],16],function(x) which.min(abs(freqDiff_median[[i]]-x)) )))

pdf(paste(outprefix,".false.positives.far.pdf",sep=""))
points<-barplot(fp_cat_count/nrOfsimulations,col=colors,legend=1:nrOfQtls,main="False positives (no interval border within 2 Mb)",xlab="coverages",ylab="average number of False Positives")
axis(1,at=points,labels=coverages,las=2)
dev.off()

fp<-data[data$judgement=="FP",]
fp_cat_count<-sapply(1:length(coverages),function(i) table(sapply(fp[fp[,2]==coverages[i],16],function(x) which.min(abs(freqDiff_median[[i]]-x)) )))

pdf(paste(outprefix,".false.positives.all.pdf",sep=""))
points<-barplot(fp_cat_count/nrOfsimulations,col=colors,legend=1:nrOfQtls,main="False positives",xlab="coverages",ylab="average number of False Positives")
axis(1,at=points,labels=coverages,las=2)
dev.off()


tableData<-data.frame(judgement=data$judgement,cov=data$cov,rank=data$rank,linked=ifelse(data$lg>0 & !is.na(data$lg),"yes","no"))

countTable<-dcast(melt(tableData,id=c("cov","judgement","rank","linked")),judgement+cov~rank+linked,length)
write.table(countTable,file=paste(outprefix,".counts.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE)





















