library(bbmle)

estimateFreq_one<-function(interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 f<-colSums(shoremap_qtlData[interval,cols])
 f[1]/sum(f)
}


ll_one_fun<-function(f,interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 nrOfPlants<-500
 fc<-min(max(floor(f*nrOfPlants),1),499)/nrOfPlants #fulhack
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==min(interval) & shoremap_qtlmem[,3]==max(interval) & shoremap_qtlmem[,4]==memory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==min(interval) & shoremap_qtlmem[,3]==max(interval)& shoremap_qtlmem[,4]==memory,5]
 }else{
  l<--sum(apply(shoremap_qtlData[interval,cols],1,function(x){dbinom(x[1],size=sum(x),prob=fc,log=TRUE)}))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,min(interval),max(interval),memory,l)),".GlobalEnv")
  l
 }
}

ll_rest_two<-function(f1,f2,minDiff){
 if(f1<0||f1>1){
  110000+abs(f1-0.5)
 }else if(f2<0 || f2>1){
  120000+abs(f2-0.5)
 }else if(abs(f1-f2)<minDiff){
  130000+abs(f1-f2)
 }else{
  ll_one_fun(f1,shoremap_interval,1)+ll_one_fun(f2,shoremap_interval,2)
 }
}


maxConf<-function(x,minDiff,level,validMin=1,validMax=Inf,maxStart=Inf,minStop=-Inf){
# print(paste(minDiff,level,sep=" ## "))
 start=floor(x[1])
 stop=floor(x[2])
 validMin<-max(1,validMin)
 validMax<-min(length(shoremap_qtlData[,1]),validMax)
# print(paste(start,stop,sep=" ## "))
 if(start>stop){
  1000000000+start-stop
 }else if(start<validMin){
  1100000000+validMin-start
 }else if(stop>validMax){
  1200000000+stop-validMax
 }else if(start>maxStart){
  1300000000+start-maxStart
 }else if(stop<minStop){
  1400000000+minStop-stop
 }else if(stop-start<3){
  1500000000+start-stop
 }else{
  interval<-start:stop
  f1<-estimateFreq_one(interval,1)
  f2<-estimateFreq_one(interval,2)
  diff<-abs(f1-f2)
  if(diff>=minDiff){
   -(shoremap_qtlData[stop,2]-shoremap_qtlData[start,2])
  }else{
   assign("shoremap_interval",interval,".GlobalEnv")
   f1e<-f1
   f2e<-f2
   if(f1e<f2e){
    f1e<-max(0,f1e-(minDiff-diff)/2)
    f2e<-f1e+minDiff
   }else{
    f1e<-min(f1e+(minDiff-diff)/2,1)
    f2e<-f1e-minDiff
   }
   restricted<-mle2(ll_rest_two,start=list(f1=f1e,f2=f2e),fixed=list(minDiff=minDiff))
   full<-ll_one_fun(f1,interval,1)+ll_one_fun(f2,interval,2)
   p<-pchisq(-2*(full-restricted@min),1)
#   print(paste(p,start,stop,sep=" ## "))
   if(p<=level){
    ll<--(shoremap_qtlData[stop,2]-shoremap_qtlData[start,2]+log(1-p))
#    print(paste(p,start,stop,shoremap_qtlData[start,2],shoremap_qtlData[stop,2],ll,sep=" ## "))
    ll
   }else{
    shoremap_qtlData[stop,2]-shoremap_qtlData[start,2]
   }
  }
 }
}

windowedScoreOld<-function(data2,shifts,minMarkersInWindow,p.valueCutoff,winSize){
 unsorted<-sapply(shifts,function(shift){
 #total count in each window
 windowsAll<-floor((data2[,2]+shift)/winSize)
 markerCount<-table(windowsAll)
 markerCount<-cbind(as.numeric(rownames(markerCount)),markerCount)
 #filter on p-value. Only use markers with p-values below the cutoff
 passedData<-data2[data2[,12]<p.valueCutoff,]
 #calculate windows
 windows<-floor((passedData[,2]+shift) /winSize)
 #only use windows with more than the minimum number of markers
 passedCount<- table(windows)
 passedCountMatrix<-cbind(as.numeric(rownames(passedCount)),passedCount)
 windowsToUse<-as.numeric(rownames(passedCount)[passedCount>minMarkersInWindow])
 markersToUse<- windows %in% windowsToUse
# markersAllToUse<- windowsAll %in% windowsToUse #NEW
# windowsAll<-windowsAll[markersAllToUse] #NEW
 passedData<-passedData[markersToUse,]
 windows<-windows[markersToUse]
 #calculate scores
 #adjusted allele score
# aas<-tapply(data2[,11][markersAllToUse],windowsAll,function(window){sum(window)})*sapply(unique(windowsAll), function(win) sum(passedCountMatrix[passedCountMatrix[,1]==win,2]))/markerCount[markerCount[,1] %in% unique(windowsAll),2]**2 #CHANGED
 #allele score
 aa<-tapply(passedData[,11],windows,function(window){sum(window)})/markerCount[markerCount[,1] %in% unique(windows),2] 
 #p score
# tp<-tapply(passedData[,12],windows,function(window){sum(window)})/markerCount[markerCount[,1] %in% unique(windows),2]
#  pos<-tapply(passedData[,2],windows,mean)
 pos<-sapply(unique(windows),function(i) (i+0.5)*winSize-shift)
# cbind(pos,tp,aa,aas)
 cbind(pos,aa,aa,aa)
},simplify=F)

 pos<-c(sapply(unsorted,function(i) i[,1]),recursive=T)
# tp<-c(sapply(unsorted,function(i) i[,2]),recursive=T)
 aa<-c(sapply(unsorted,function(i) i[,3]),recursive=T)
# aas<-c(sapply(unsorted,function(i) i[,4]),recursive=T)
 o<-sort.int(pos,index.return=T)
# sorted<-t(sapply(o$ix,function(i) c(pos[i],tp[i],aa[i],aas[i])))
 t(sapply(o$ix,function(i) c(pos[i],aa[i],aa[i],aa[i])))
}


windowedScore<-function(data2,winSize,winStep,minMarkersInWindow,p.valueCutoff){
 if(winSize %% winStep !=0){
  print("winsize is better an even multiple of winstep")
  windowedScoreOld(data2,seq(0,winSize-1,winStep),minMarkersInWindow,p.valueCutoff,winSize)
 }else{
  winsubrate<-winSize/winStep
  subwindows<-floor(data2[,2]/winStep)
  uniqWindows<-unique(subwindows)
  index<-(min(uniqWindows)-winsubrate+1):(max(uniqWindows))
  allIndex<-(min(uniqWindows)-winsubrate+1):(max(uniqWindows)+winsubrate-1)
  missingWindows<-allIndex[!(allIndex %in% uniqWindows)]
  extractOrder<-sort(c(uniqWindows,missingWindows),index.return=TRUE)$ix
  #generate positions
  pos<-index*winStep+winSize/2
  #generate total counts
  subCount<-table(subwindows)
  subCount<-c(subCount,rep(0,length(missingWindows)))[extractOrder]
  totCount<-sapply(1:length(index),function(i) sum(subCount[i:(i+winsubrate-1)]))
  
  passedData<-data2[data2[,12]<p.valueCutoff,]
  passedWindows<-floor(passedData[,2]/winStep)
  uniqPassedWindows<-unique(passedWindows)
  missingPassedWindows<-allIndex[!(allIndex %in% uniqPassedWindows)]
  extractOrderPassed<-sort(c(uniqPassedWindows,missingPassedWindows),index.return=TRUE)$ix
  #generate passed counts
  passedSubCount<-table(passedWindows)
  passedSubCount<-c(passedSubCount,rep(0,length(missingPassedWindows)))[extractOrderPassed]
  passedCount<-sapply(1:length(index),function(i) sum(passedSubCount[i:(i+winsubrate-1)]))
  #generate sums
  subSum<-tapply(passedData[,11],passedWindows,sum)
  subSum<-c(subSum,rep(0,length(missingPassedWindows)))[extractOrderPassed]
  passedSum<-sapply(1:length(index),function(i) sum(subSum[i:(i+winsubrate-1)]))
  aa<-passedSum/totCount

  cbind(pos,aa,aa,aa)[passedCount>minMarkersInWindow,]
 }
}

predictAllPeaksSub<-function(sorted,winSize,direction,doPlot=FALSE,span=0.17){
 direction<-sign(direction)
 if(direction==0){
  direction=1
 }
 #predicting the number of peaks
# span<-nrWin*max(sorted[,1])/winSize/nrow(sorted)
# print(span)
 y.loess<- loess(y~x, span=span,data.frame(x=sorted[,1],y=sorted[,4]))
 y.predict<-direction*predict(y.loess,data.frame(x=sorted[,1]))
 if(doPlot){
  lines(sorted[,1],direction*y.predict,col="steelblue")
 }

 maxPeakDiff<-sign(diff(c(-Inf,y.predict,-Inf)))
 maxPeakID<-which(maxPeakDiff!=0)
 maxPeakIndex<-maxPeakID[which(diff(maxPeakDiff[maxPeakID])==-2)]
 maxPeakCount<-length(maxPeakIndex)
 minPeakDiff<-sign(diff(c(-Inf,-y.predict,-Inf)))
 minPeakID<-which(minPeakDiff!=0)
 minPeakIndex<-minPeakID[which(diff(minPeakDiff[minPeakID])==-2)]
 minPeakCount<-length(minPeakIndex)

 #calculate heights
 
 lMinIndex<-minPeakIndex[1:(maxPeakCount-(minPeakIndex[1]>=maxPeakIndex[1]))]
 if(minPeakIndex[1]>=maxPeakIndex[1]){
  lMinIndex<-c(maxPeakIndex[1],lMinIndex)
 }
 lDiff <- y.predict[maxPeakIndex]-y.predict[lMinIndex]

 rMinIndex<-minPeakIndex[(minPeakCount-maxPeakCount+1+(minPeakIndex[minPeakCount]<=maxPeakIndex[maxPeakCount])):minPeakCount]
 if(minPeakIndex[minPeakCount]<=maxPeakIndex[maxPeakCount]){
  rMinIndex<-c(rMinIndex,maxPeakIndex[maxPeakCount])
 }
 rDiff<- y.predict[maxPeakIndex]-y.predict[rMinIndex]
 cbind(maxPeakIndex,direction*y.predict[maxPeakIndex],lDiff,rDiff,lMinIndex,rMinIndex,rep(direction,length(maxPeakIndex)))
}

predictAllPeaks<-function(sorted,winSize,doPlot=FALSE,span=0.17){
 minimas<-predictAllPeaksSub(sorted,winSize,-1,doPlot,span)
 maximas<-predictAllPeaksSub(sorted,winSize,1,FALSE,span)
 peaks<-rbind(maximas[maximas[,2]>0,],minimas[minimas[,2]<0,])
 matrix(peaks[sort(peaks[,1],index.return=TRUE)$ix,],ncol=7)
}


predictPeaks_p<-function(x,cutOffs){
 max(sum(cutOffs[,1]>x[2]), sum(cutOffs[,2]>max(x[3:4])))/nrow(cutOffs)
}

predictPeaks<-function(sorted,minFrequencyChange,minAbsolutePeak,cutOffs,winSize,doPlot,span){
 peaks<-predictAllPeaks(sorted,winSize,doPlot,span)
 
 mPeak<-peaks[,3]>minFrequencyChange | peaks[,4]>minFrequencyChange
 
 if(sum(mPeak)>0){
  p<-sapply(which(mPeak),function(i) predictPeaks_p(peaks[i,],cutOffs))
  cbind(matrix(c(peaks[mPeak,c(1,5:7)]),ncol=4),p)
 }else if(max(peaks[,2])>minAbsolutePeak){
  #treat it as one big peak
  p<-predictPeaks_p(peaks[which.max(abs(peaks[,2])),],cutOffs)
  matrix(c(peaks[which.max(abs(peaks[,2])),1],min(peaks[,5]),max(peaks[,6]),peaks[which.max(abs(peaks[,2])),7],p),ncol=5)
 }else{
  matrix(c(peaks[mPeak,c(1,5:7)]),ncol=5)
 }
}

optimFn<-function(x,winSize,low,high,direction){
 if(x>=low && x<=high){
  markers<-which(shoremap_qtlData[,2]>=x-winSize/2 & shoremap_qtlData[,2]<=x+winSize/2)
  -direction*(estimateFreq_one(markers,1)-estimateFreq_one(markers,2))
 }else{
  616
 }
}

optimGr<-function(x,winSize,low,high,direction){
 if(x<low){
  direction*(x-low)
 }else if(x>high){
  direction*(x-high)
 }else{
  markers<-which(shoremap_qtlData[,2]>=x-winSize/2 & shoremap_qtlData[,2]<=x+winSize/2)
  -direction*lm(y~x,weights=rowSums(shoremap_qtlData[markers,c(4:6,8:10)]),data=data.frame(x=shoremap_qtlData[markers,2],y=shoremap_qtlData[markers,11]))$coefficients[2]
 }
}






#Run.....

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
 print("Usage:\nR --args seqfile1 seqfile2 outprefix < SHOREmap_qtlseq.R\nor\nR --args seqfile1 seqfile2 qtlfile outprefix < SHOREmap_qtlseq.R")
 q("no",616,FALSE)
}


#Parameters
winSize<-1000000
winStep<-10000
p.valueCutoff<-0.2
minMarkersInWindow<-10
minCoverage <- 4
maxCoverage <- 300
minFrequencyChange<-0#0.05
minAbsolutePeak<-0#0.05
bootstrap=1000

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

data<- cbind(hs[,1:2],hs.freq,hs[,3:5],ls.freq,ls[,3:5],hs.freq-ls.freq)

#p<-apply(data,1,function(x) fisher.test(matrix(as.numeric(x[c(4,5,8,9)]),nrow=2))$p.value)
p<-rep(0,nrow(data))

data<-cbind(data,p)

shifts<-seq(0,winSize-1,winStep)

#sapply(unique(data[,1]),function(chr){
#system.time({
cutOffs<-matrix(c(sapply(unique(data[,1]),function(chr){
 data2<-data[data[,1]==chr,]
 nrOfMarkers<-nrow(data2)
 nrInWindow<-round(nrOfMarkers/max(data2[,2])*winSize)
 replicate(bootstrap,{
#  p<-data2[,2]
#  i<-sample(1:(nrOfMarkers-nrInWindow+1),size=ceiling(nrOfMarkers/nrInWindow),replace=TRUE)
#  i<-sapply(i,function(x) x:(x+nrInWindow-1))[1:nrOfMarkers]
  d2<-data2
  d2[,11]<-d2[,11]*sample(c(1,-1),nrOfMarkers,replace=TRUE)
  sorted<-windowedScore(d2,winSize,winStep,minMarkersInWindow,p.valueCutoff)
  sorted<-sorted[sorted[,1]>min(data2[,2])+winSize/2 & sorted[,1]<max(data2[,2])-winSize/2,]
  peaks<-predictAllPeaks(sorted,winSize,FALSE)
#  lines(sorted[,c(1,4)],col="red")  
  rbind(peaks[,2],pmax(peaks[,3],peaks[,4]))
 },simplify=TRUE)
}),recursive=TRUE),ncol=2,byrow=TRUE)
#})

#cutOffs2<-c(sapply(unique(data[,1]),function(chr){
# data2<-data[data[,1]==chr,]
# nrOfMarkers<-nrow(data2)
# nrInWindow<-round(nrOfMarkers/max(data2[,2])*winSize)
# replicate(bootstrap,{
#  d2<-data2
#  d2[,11]<-d2[,11]*sample(c(1,-1),nrOfMarkers,replace=TRUE)
#  sorted<-windowedScore(d2,winSize,winStep,minMarkersInWindow,p.valueCutoff)
#  sorted[sorted[,1]>min(data2[,2])+winSize/2 & sorted[,1]<max(data2[,2])-winSize/2,3]
# },simplify=TRUE)
#}),recursive=TRUE)



#cutOffs<- matrix(rep(0,20000),ncol=2)

cutOffs<-cbind(sort(abs(cutOffs[,1])),sort(cutOffs[,2]))

minAbsolutePeak<-max(cutOffs[round(nrow(cutOffs)*0.99),1],minAbsolutePeak)
minFrequencyChange<-max(cutOffs[round(nrow(cutOffs)*0.99),2],minFrequencyChange)



estimates<-c()

nrOfChrs<-length(unique(data[,1]))

pdf(paste(outprefix,".plots.pdf",sep=""))
par(mfrow=c(ceiling(nrOfChrs/2),2))



for(chr in unique(data[,1])){
 data2<-data[data[,1]==chr,]
# png()
 plot(data2[,2],data2[,11],main=paste("chr",chr),xlab="pos",ylab="frequncy difference",ylim=c(-1,1))


 sorted<-windowedScore(data2,winSize,winStep,minMarkersInWindow,p.valueCutoff)
 sorted<-sorted[sorted[,1]>min(data2[,2])+winSize/2 & sorted[,1]<max(data2[,2])-winSize/2,]
 #for each shift value, calculate the p score for each window
 lines(sorted[,c(1,4)],col="red")
 
 peaks<-predictPeaks(sorted,minFrequencyChange,minAbsolutePeak,cutOffs,winSize,TRUE,0.2)
 if(nrow(peaks)>0){
  for(peakIndex in 1:nrow(peaks)){
   #extract region
   peak<-peaks[peakIndex,1]
   direction<-peaks[peakIndex,4]
   assign("shoremap_qtlmem",matrix(c(-1,-1,-1,-1,-1),nrow=1),".GlobalEnv")
   lowerBoundary<-sorted[peaks[peakIndex,2],1]
   upperBoundary<-sorted[peaks[peakIndex,3],1]
   data3<-data2[data2[,2]>=lowerBoundary&data2[,2]<=upperBoundary,]
   assign("shoremap_qtlData",data3,".GlobalEnv")


   #pinpointing window (point estimation)
   md_long<-sapply(shifts,function(shift){
    windows<-floor((data3[,2]+shift)/winSize)
    d<-tapply(1:length(data3[,1]),windows,function(interval){
     estimateFreq_one(interval,1)-estimateFreq_one(interval,2)
    })
    s<-table(windows)
    interval<-which(windows==unique(windows)[which.max(direction*d)])
    list(md=d,best=direction*max(direction*d),interval=interval,size=s)
   })
   md<-c(md_long[1,],recursive=TRUE)
   ms<-c(md_long[4,],recursive=TRUE)
   mp<-c(sapply(shifts,function(shift){
    windows<-floor((data3[,2]+shift)/winSize)
  #  tapply(data3[,2],windows,mean)
    sapply(unique(windows),function(i) (i+0.5)*winSize-shift)
   }),recursive=TRUE)
   wins<-t(sapply(sort(mp,index.return=TRUE)$ix,function(i) c(mp[i],md[i],ms[i])))
   wins<-wins[wins[,1]>lowerBoundary+winSize/2 & wins[,1]<upperBoundary-winSize/2 & wins[,3]>minMarkersInWindow,1:2]
   lines(wins[,1],wins[,2],col="violetred")

 
   #adjust frequency
   minIndex<-which.max(direction*wins[,2])
   minDiff<-wins[minIndex,2]
   minDiff<-abs(floor(minDiff*1000)/1000)

#   interval<-md_long[3,which(c(md_long[2,],recursive=TRUE)==wins[minIndex,2])[1]]$interval
   interval<-which(data3[,2]>=wins[minIndex,1]-winSize/2 & data3[,2]<wins[minIndex,1]+winSize/2)
   roughEst<-round(mean(data3[interval,2]))

  

   #identify window
   opt<-optim(fn=maxConf,par=c(min(interval),max(interval)),minDiff=minDiff,level=0.99)
   bestValue<-Inf
   while(opt$value<bestValue){
    bestValue<-opt$value
    opt<-optim(fn=maxConf,par=c(floor(opt$par[1]),floor(opt$par[2])),minDiff=minDiff,level=0.99)
   }
   rect(data3[opt$par[1],2],0.96,data3[opt$par[2],2],1,col="limegreen",border="limegreen")

   est<-round(optim(par=roughEst,fn=optimFn,gr=optimGr,winSize=winSize,low=data3[opt$par[1],2],high=data3[opt$par[2],2],direction=direction)$par[1])
   abline(v=est,col="steelblue")

   estimates<-rbind(estimates,c(chr,data3[opt$par[1],2],data3[opt$par[2],2],minDiff,est))
  }
 }
 for(pos in qtls[qtls$chr==chr,1]){
  abline(v=pos,col="red")
 }
}

dev.off()

#print estimates
if(length(estimates[1,])>0){
 colnames(estimates)<-c("#chr","start","stop","freqDiff","roughEst")
 write.table(estimates,file=paste(outprefix,".qtlEstimates.csv",sep=""),row.names=F,sep="\t",quote=F)
}

q("no",0,F)
