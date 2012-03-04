#lib v01
#v02 - included the absolute height as cutoff for peak prediction

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
 
 mPeak<- (peaks[,3]>minFrequencyChange | peaks[,4]>minFrequencyChange) & abs(peaks[,2])>minAbsolutePeak
 
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



