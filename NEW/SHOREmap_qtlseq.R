#v09 - changes in ploting
#v10 - added the true interval of the qtl

source("~/shoreMap/NEW/SHOREmap_qtlseq_lib.R")

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
tolerance=2e6
winSize<-1000000
winStep<-10000
p.valueCutoff<-0.2
minMarkersInWindow<-10
minCoverage <- 0
maxCoverage <- 3000000
minFrequencyChange<-0.05
minAbsolutePeak<-0.05
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
 qtls<-read.table(paste(qtl_file,".qtl.txt",sep=""))
 ls.t<-read.table(paste(qtl_file,".low.AF.txt",sep=""))
 hs.t<-read.table(paste(qtl_file,".high.AF.txt",sep=""))
 a.t<-cbind(ls.t[,1],ls.t[,2],(hs.t[,3]-ls.t[,3])/100)
 marker<-qtls[,2]+1
 freqDiff<-a.t[marker,3]
 trueStart<-sapply(1:length(freqDiff),function(i){
  start<-marker[i]
  while(a.t[start,3]==freqDiff[i]&&start>1){
   start<-start-1
  }
  a.t[start,2]+1
 })
 trueEnd<-sapply(1:length(freqDiff),function(i){
  end<-marker[i]
  while(a.t[end,3]==freqDiff[i]&& end<nrow(a.t)){
   end<-end+1
  }
  a.t[end,2]-1
 })
 rank<-sort(sort(abs(qtls[,5]),index.return=TRUE,decreasing=TRUE)$ix,index.return=TRUE)$ix
 lg<-sapply(qtls[,3],function(x) ifelse(sum(qtls[,3]==x)>1,x,0))
 lgType<-sapply(lg,function(x) ifelse(x==0,0,{
  tt<-table(sign(qtls[lg==x,5]))
  ifelse(length(tt)==1,sign(sum(qtls[lg==x,5])),ifelse(tt[1]<tt[2],-tt[1]/tt[2],ifelse(tt[1]==tt[2],0,tt[2]/tt[1])))
 }))
 qtls<-data.frame(id=qtls[,1],pos=qtls[,4],chr=qtls[,3],effect=qtls[,5],rank=rank,lg=lg,lgType=lgType,trueFreqDiff=freqDiff,trueStart=trueStart,trueEnd=trueEnd)
 
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



estimates<- matrix(ncol=5)[FALSE,]
colnames(estimates)<-c("chr","start","stop","freqDiff","Est")

nrOfChrs<-length(unique(data[,1]))

pdf(paste(outprefix,".plots.pdf",sep=""))
par(mfrow=c(ceiling(nrOfChrs/2),2))

colors<-rainbow(50,start=4/6,end=1/6)

for(chr in unique(data[,1])){
 data2<-data[data[,1]==chr,]
# png()
 plot(data2[,2],data2[,11],main=paste("chr",chr),xlab="pos",ylab="frequency difference",ylim=c(-1.2,1),pch=16,cex=0.75,col="lightsteelblue3")


 sorted<-windowedScore(data2,winSize,winStep,minMarkersInWindow,p.valueCutoff)
 sorted<-sorted[sorted[,1]>min(data2[,2])+winSize/2 & sorted[,1]<max(data2[,2])-winSize/2,]
 #for each shift value, calculate the p score for each window
 lines(sorted[,c(1,4)],col="limegreen")
 
 peaks<-predictPeaks(sorted,minFrequencyChange,minAbsolutePeak,cutOffs,winSize,FALSE,0.2)
 if(nrow(peaks)>0){
  for(peakIndex in 1:nrow(peaks)){
   #extract region
   peak<-peaks[peakIndex,1]
   direction<-peaks[peakIndex,4]
   assign("shoremap_qtlmem",matrix(c(-1,-1,-1,-1,-1),nrow=1),".GlobalEnv")
   lowerBoundary<-sorted[peaks[peakIndex,2],1]
   upperBoundary<-sorted[peaks[peakIndex,3],1]
   if(upperBoundary-lowerBoundary>1.5*winSize){
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
    #lines(wins[,1],wins[,2],col="violetred")

 
    #adjust frequency
    minIndex<-which.max(direction*wins[,2])
    minDiff<-wins[minIndex,2]
    minDiff<-abs(floor(minDiff*1000)/1000)

#    interval<-md_long[3,which(c(md_long[2,],recursive=TRUE)==wins[minIndex,2])[1]]$interval
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
    abline(v=est,col="limegreen")

    estimates<-rbind(estimates,c(chr,data3[opt$par[1],2],data3[opt$par[2],2],minDiff,est))
   }
  }
 }
 if(sum(qtls$chr==chr)>0){
  apply(qtls[qtls$chr==chr,],1,function(x) {
   e<-min(50,x[4])
   e<-max(-50,e)
   e<-round((e+50)/2)
   abline(v=x[2],col=colors[e])
  })
 }
 #print scale
 chrStart<-min(data2[,2])
 chrEnd<-max(data2[,2])
 chrSize<-chrEnd-chrStart
 scaleStart<-chrStart+chrSize/5
 scaleEnd<-chrEnd-chrSize/5
 scaleSize<-scaleEnd-scaleStart
 dx<-scaleSize/50
 starts<-seq(scaleStart,scaleEnd-1,dx)  
  
 rect(starts, -1, starts+0.9*dx, -1.1, col = colors, border = "grey");
 text(scaleStart+scaleSize/2,-0.95,"effects")
 text(scaleStart+scaleSize/2,-1.15,"0")
 text(starts[1],-1.15,"-50")
 text(max(starts)+0.9*dx,-1.15,"50")

}

dev.off()


#print estimates
header<-c("chr",colnames(qtls)[c(1:2,4:ncol(qtls))],colnames(estimates)[c(2:ncol(estimates))],"judgement","spec")
header[1]<-paste("#",header[1],sep="")
qhc<-ncol(qtls)-1
ehc<-ncol(estimates)-1

toPrint<-sapply(unique(data[,1]),function(chr){
 cq<-sum(qtls$chr==chr)
 ce<-sum(estimates[,1]==chr)
 q<-matrix(c(qtls[qtls$chr==chr,c(1:2,4:ncol(qtls))],recursive=TRUE),ncol=qhc)
 e<-matrix(estimates[estimates[,1]==chr,c(2:ncol(estimates))],ncol=ehc)
 if(cq>0 && ce>0){
  #get True Positive
  pairs<-combn(cq+ce,2)
  pairs<-matrix(pairs[,pairs[1,]<=cq & pairs[2,]>cq],nrow=2)
  pairs[2,]<-pairs[2,]-cq
  tp<-apply(pairs,2,function(x){
   q[x[1],2]>=e[x[2],1] && q[x[1],2]<=e[x[2],2]
  })
  toPrint<-matrix(apply(matrix(pairs[,tp],nrow=2),2,function(x) c(chr,q[x[1],],e[x[2],],"TP","")),byrow=TRUE,ncol=length(header))
  #check if any of the qtls are close to a interval
  close<-apply(pairs,2,function(x){
   q[x[1],2]+tolerance>=e[x[2],1] && q[x[1],2]-tolerance<=e[x[2],2]
  })
  toPrint<-rbind(toPrint,matrix(apply(matrix(pairs[,close &!tp],nrow=2),2,function(x) c(chr,q[x[1],],e[x[2],],"FP","close")),byrow=T,ncol=length(header)))
  #more FP?
  for(i in 1:ce){
   if(!i %in% unique(pairs[2,close]) ){
    toPrint<-rbind(toPrint,matrix(c(chr,rep(NA,qhc),e[i,],"FP","" ),ncol=length(header)))
   }
  }
  #FN
  for(i in 1:cq){
   if(!i %in% unique(pairs[1,close]) ){
    toPrint<-rbind(toPrint,matrix(c(chr,q[i,],rep(NA,ehc),"FN","" ),ncol=length(header)))
   }
  }
  toPrint
 }else if(cq>0){
  #False Negative
  matrix(apply(q,1,function(x) matrix(c(chr,x,rep(NA,ehc),"FN","" ),ncol=length(header))),ncol=length(header),byrow=TRUE)
 }else if(ce>0){
  #False Positive
  matrix(apply(e,1,function(x) matrix(c(chr,rep(NA,qhc),x,"FP","" ),ncol=length(header))),ncol=length(header),byrow=TRUE)
 }else{
  #True Negative
  matrix(c(chr,rep(NA,qhc+ehc),"TN","" ),ncol=length(header))
 }
})

if(is.list(toPrint)){
 toPrint<-do.call(rbind,toPrint)
}else{
 toPrint<-matrix(toPrint,ncol=length(header),byrow=TRUE)
}

colnames(toPrint)<-header

write.table(toPrint,file=paste(outprefix,".qtlEstimates.csv",sep=""),row.names=F,sep="\t",quote=F,na="")

#if(length(estimates[1,])>0){
# colnames(estimates)<-c("chr","start","stop","freqDiff","roughEst")
# write.table(estimates,file=paste(outprefix,".qtlEstimates.csv",sep=""),row.names=F,sep="\t",quote=F)
#}

q("no",0,F)
