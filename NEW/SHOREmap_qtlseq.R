

args<-commandArgs(trailingOnly=TRUE)

#command line arguments
high_file<-args[1]
low_file<-args[2]

if(length(args)==4){
 qtl_file<-args[3]
 
 outprefix<-args[4]
}else if(length(args)==3){
 outprefix<-args[3]
 qtl_file=""
}else{
 print("Usage:\nR --args seqfile1 seqfile2 qtlfile outprefix < SHOREmap_qtlseq.R\nor\nR --args seqfile1 seqfile2 qtlfile outprefix < SHOREmap_qtlseq.R")
 q("no",616,FALSE)
}


#Parameters
winSize<-500000
winStep<-50000
minMarkersInWindow<-10
p.valueCutoff<-0.2
minCoverage <- 4
maxCoverage <- 300

#prep data
hs<-read.table(high_file)
ls<-read.table(low_file)

#extract instersection of markers
chrmod<-10**ceiling(log10(max(hs[,2],ls[,2])))
hs.mod<-hs[,1]*chrmod+hs[,2]
ls.mod<-ls[,1]*chrmod+ls[,2]

intersect.mod<-intersect(hs.mod,ls.mod)
hs<-hs[hs.mod %in% intersect.mod,]
ls<-ls[ls.mod %in% intersect.mod,]


hs.cov<-rowSums(hs[,3:5])
ls.cov<-rowSums(ls[,3:5])

goodCov<-(hs.cov>minCoverage & hs.cov<maxCoverage)&(ls.cov>minCoverage &ls.cov<maxCoverage)

hs<-hs[goodCov,]
ls<-ls[goodCov,]
qtls<-data.frame(pos=-1,chr=-1)
if(qtl_file==""){
 qtls<-data.frame(pos=-1,chr=-1)
}else{
 qtls<-read.table(qtl_file)
 qtls<-data.frame(pos=qtls[,4],chr=qtls[,3])
}


##merge
hs.freq<-hs[,3]/rowSums(hs[,3:5])
ls.freq<-ls[,3]/rowSums(ls[,3:5])

data<- cbind(hs[,1:2],hs.freq,hs[,3:5],ls.freq,ls[,3:5])

#fisher test
p<-apply(data,1,function(x) fisher.test(matrix(as.numeric(x[c(4,5,8,9)]),nrow=2))$p.value)

pdata<-cbind(data[,1:2],p,data[,3]-data[,7])


estimates<-c()

for (chr in unique(data[,1])){
 x<-which(data[,1]==chr)
 shifts<-seq(0,winSize-1,winStep)

 #for each shift value, calculate the p score for each window
 unsorted<-sapply(shifts,function(shift){
  #total count in each window
  markerCount<-table(floor((pdata[x,2]+shift)/winSize))
  markerCount<-cbind(as.numeric(rownames(markerCount)),markerCount)
  #filter on p-value. Only use markers with p-values below the cutoff
  passedData<-pdata[x,][pdata[x,3]<p.valueCutoff,]
  #calculate windows
  windows<-floor((passedData[,2]+shift) /winSize)
  #only use windows with more than the minimum number of markers
  passedCount<- table(windows)
  windowsToUse<-as.numeric(rownames(passedCount)[passedCount>minMarkersInWindow])
  markersToUse<- windows %in% windowsToUse
  passedData<-passedData[markersToUse,]
  windows<-windows[markersToUse]
  #calculate scores
#  aa<-tapply(passedData[,4],windows,function(window){sum(window)})/markerCount[markerCount[,1] %in% unique(windows),2] #needed for calculation of allele
  tp<-tapply(passedData[,3],windows,function(window){(sum(window)^2)/(sum(window))})/markerCount[markerCount[,1] %in% unique(windows),2]
#  pos<-(unique(windows)*winSize)+(winSize+shift)/2
  pos<-tapply(passedData[,2],windows,mean)
#  cbind(pos,tp,aa) #needed for calculation of allele
  cbind(pos,tp)
 },simplify=F)
 #extract data and sort
 pos<-c(sapply(unsorted,function(i) i[,1]),recursive=T)
 tp<-c(sapply(unsorted,function(i) i[,2]),recursive=T)
# aa<-c(sapply(unsorted,function(i) i[,3]),recursive=T) #needed for calculation of allele
 o<-sort.int(pos,index.return=T)
# sorted<-t(sapply(o$ix,function(i) c(pos[i],tp[i],aa[i]))) #needed for calculation of allele
 sorted<-t(sapply(o$ix,function(i) c(pos[i],tp[i])))
 if(dim(sorted)[2]>0){
  #second smoothing
  pscore.loess<-loess(pscore~x,span=0.6,data.frame(x=sorted[,1],pscore=sorted[,2]))
  pscore.predict<- predict(pscore.loess,data.frame(x=sorted[,1]))

  #plot and print
  png(paste(outprefix,"_chr",chr,"_winsize",winSize,"_winstep",winStep,".png",sep=""))

  plot(sorted[,1],sorted[,2],type="l",main=paste("pointy peak, chr",chr,"(winsize: ",winSize," bp, winstep: ",winStep," bp)"),xlab="pos",ylab="p score")


  lines(sorted[,1],pscore.predict,col="green",type="l")

  #identify local min and max values
  maxDiff<-sign(diff(c(-Inf,pscore.predict,-Inf)))
  maxID<-which(maxDiff!=0)
  maxIndex<-maxID[which(diff(maxDiff[maxID])==-2)]
  minDiff<-sign(diff(c(-Inf,-pscore.predict,-Inf)))
  minID<-which(minDiff!=0)
  minIndex<-minID[which(diff(minDiff[minID])==-2)]

  if(length(maxIndex)>0&&length(minIndex)>0){
   #filter peaks, a good peak must be twice as large as the flanking minimas and at least 1% of the total maxima on the chromosome
   minShift<-0
   if(maxIndex[1]<minIndex[1]){
    #starts with a local maxima
    minShift<--1
   }else{
    #starts with a local minima
    minShift<-0
   }
   for(i in 1:length(maxIndex)){
    abline(v=sorted[maxIndex[i],1],col="blue")
    if(pscore.predict[maxIndex[i]]*100>max(pscore.predict[maxIndex])){
     checkLow<-FALSE
     checkHigh<-FALSE
     minI<-i+minShift
     if(minI>0&minI<=length(minIndex)){
      checkLow<-pscore.predict[maxIndex[i]]>2*pscore.predict[minIndex[minI]]
     }else{
      checkLow<-TRUE
     }
     minI<-minI+1
     if(minI>0&minI<=length(minIndex)){
      checkHigh<-pscore.predict[maxIndex[i]]>2*pscore.predict[minIndex[minI]]
     }else{
      checkHigh<-TRUE
     }
     if(checkLow && checkHigh){
      estimates<-rbind(estimates,c(chr,sorted[maxIndex[i],1]))
      points(sorted[maxIndex[i],1],pscore.predict[maxIndex[i]],col="blue",pch=8,lwd=3,cex=2)
     }
    }
   }
  }
  
  #mark the true QTLs
  for(pos in qtls[qtls$chr==chr,1]){
   abline(v=pos,col="red")
  }
 
  dev.off()
 }
}

#print estimates
if(length(estimates[1,])>0){
 colnames(estimates)<-c("#chr","position")
 write.table(estimates,file=paste(outprefix,"_qtlEstimates.csv",sep=""),row.names=F,sep="\t",quote=F)
}

q("no",0,F)
