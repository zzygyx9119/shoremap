library(bbmle)

mult_two_fun <- function(locus,r,limit=0.9,memmory=3){
 if((locus%%1000)==0){
  print(locus)
 }
 dist<-abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2])
 assign("shoremap_q",sapply(dist,q_fun,r),".GlobalEnv")
# assign("shoremap_mod",ppois(dist,r,lower.tail=F),".GlobalEnv")
 assign("shoremap_toUse",shoremap_q>limit,".GlobalEnv")
# assign("shoremap_toUse",ppois(dist,r,lower.tail=F)>0,".GlobalEnv")
 if(sum(shoremap_toUse)>1){
  mle2(ll_two_fun,optimizer="optimize",start=list(f=0.5),fixed=list(locus=locus,r=r,memmory=memmory),lower=0,upper=1)@coef
 }else{
  sum(shoremap_qtlData[locus,c(4,8)])/sum(shoremap_qtlData[locus,c(4,5,8,9)])
 }

}

mult_one_fun<-function(locus,r,limit=0.9,memmory=1){
 if((locus%%1000)==0){
  print(locus)
 }
 cols=c(4,5)
 if(memmory==2){
  cols=c(8,9)
 }
 dist<-abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2])
 assign("shoremap_q",sapply(dist,q_fun,r),".GlobalEnv")
# assign("shoremap_mod",ppois(dist,r,lower.tail=F),".GlobalEnv")
 assign("shoremap_toUse",shoremap_q>limit,".GlobalEnv")
# assign("shoremap_toUse",ppois(dist,r,lower.tail=F)>0,".GlobalEnv")
 if(sum(shoremap_toUse)>1){
  mle2(ll_one_fun,optimizer="optimize",start=list(f=0.5),fixed=list(locus=locus,r=r,memmory=memmory),lower=0,upper=1)@coef
 }else{
  shoremap_qtlData[locus,cols[1]]/sum(shoremap_qtlData[locus,cols])
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
    print(paste(p,start,stop,shoremap_qtlData[start,2],shoremap_qtlData[stop,2],ll,sep=" ## "))
    ll
   }else{
    shoremap_qtlData[stop,2]-shoremap_qtlData[start,2]
   }
  }
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

estimateFreq_one<-function(interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 f<-colSums(shoremap_qtlData[interval,cols])
 f[1]/sum(f)
}

estimateFreq_two<-function(interval){
 cols<-c(4,5,8,9)
 f<-colSums(shoremap_qtlData[interval,cols])
 (f[1]+f[3])/sum(f)
}


ll_one_fun<-function(f,interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
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

ll_two_fun<-function(f,interval,memory=3){
 cols1<-c(4,5)
 cols2<-c(8,9)
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else{
  l<- ll_one_fun(fc,interval,1)+ll_one_fun(fc,interval,2)
  l
 }
}

ll_two_fun_withQ<-function(f,locus,r=50000,memory=3){
 cols1<-c(4,5)
 cols2<-c(8,9)
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory,4]
 }else{
  l<- ll_one_fun_withQ(fc,locus,r,1)+ll_one_fun_withQ(fc,locus,r,2)
#  print(paste(fc,l,sep=" ## "))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,locus,memmory,l)),".GlobalEnv")
  l
 }
}

ll_one_fun_withQ<-function(f,locus,r=50000,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
# print(summary(q))
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory,4]
 }else{
#  r<- 250000*100/nrOfPlants
#  q<- sapply(abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2]),q_fun,r)
  l<- -sum(apply(cbind(shoremap_q[shoremap_toUse],shoremap_qtlData[shoremap_toUse,cols]),1,function(x) pa_fun(fc,x[1],x[2],x[3])))
#  print(paste(fc,l,sep=" ## "))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,locus,memmory,l)),".GlobalEnv")
  l
 }
}

pa_fun<-function(f,q,p1,p2){
 #mod: supposed to be the probability that less than
 res<-0 
 if(q==0){
  #0
  res<-(1-f)^p1*f^p2
 }else{
  res<-(q*f^p1*(1-f)^p2+(1-q)*(1-f)^p1*f^p2)
 }
 ifelse(res>0,log(res),-61600)
}

q_fun<-function(d,r){
 exp(-d/r)
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
winSize<-500000
winStep<-10000
minMarkersInWindow<-10
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

data<- cbind(hs[,1:2],hs.freq,hs[,3:5],ls.freq,ls[,3:5],abs(hs.freq-ls.freq))

estimates<-c()

nrOfChrs<-length(unique(data[,1]))

pdf(paste(outprefix,".plots.pdf",sep=""))
par(mfrow=c(ceiling(nrOfChrs/2),2))

for(chr in unique(data[,1])){
 data2<-data[data[,1]==chr,]
# png()

 #predicting the number of peaks
 plot(data2[,2],data2[,11],main=paste("chr",chr),xlabel="pos",ylabel="frequncy difference")
 y.loess<- loess(y~x, span=0.4,data.frame(x=data2[,2],y=data2[,11]))
 y.predict<-predict(y.loess,data.frame(x=data2[,2]))
 lines(data2[,2],y.predict,col="skyblue")

 maxPeakDiff<-sign(diff(c(-Inf,y.predict,-Inf)))
 maxPeakID<-which(maxPeakDiff!=0)
 maxPeakIndex<-maxPeakID[which(diff(maxPeakDiff[maxPeakID])==-2)]
 minPeakDiff<-sign(diff(c(-Inf,-y.predict,-Inf)))
 minPeakID<-which(minPeakDiff!=0)
 minPeakIndex<-minPeakID[which(diff(minPeakDiff[minPeakID])==-2)]

 for(peak in maxPeakIndex){
  #extract region
  assign("shoremap_qtlmem",matrix(c(-1,-1,-1,-1,-1),nrow=1),".GlobalEnv")
  lowerBoundary<-max(c(minPeakIndex[minPeakIndex<peak],1))
  upperBoundary<-min(c(minPeakIndex[minPeakIndex>peak],length(data2[,1])))
  data3<-data2[lowerBoundary:upperBoundary,]
  shifts<-seq(0,winSize-1,winStep)
  assign("shoremap_qtlData",data3,".GlobalEnv")


  #pinpointing window (point estimation)
  md_long<-sapply(shifts,function(shift){
   windows<-floor((data3[,2]+shift)/winSize)
   d<-tapply(1:length(data3[,1]),windows,function(interval){
    abs(estimateFreq_one(interval,1)-estimateFreq_one(interval,2))
   })
   interval<-which(windows==unique(windows)[which.max(d)])
   list(md=d,best=max(d),interval=interval)
  })
  md<-c(md_long[1,],recursive=TRUE)
  mp<-c(sapply(shifts,function(shift){
   windows<-floor((data3[,2]+shift)/winSize)
   tapply(data3[,2],windows,mean)
  }),recursive=TRUE)
  wins<-t(sapply(sort(mp,index.return=TRUE)$ix,function(i) c(mp[i],md[i])))
  lines(wins[,1],wins[,2],col="violetred")

  #adjust frequency
  minDiff<-max(wins[,2])
  minDiff<-floor(minDiff*1000)/1000
  interval<-md_long[3,which.max(c(md_long[2,],recursive=TRUE))]$interval
  roughEst<-mean(c(min(interval),max(interval)))
  abline(v=roughEst,col="steelblue")  

  #identify window
  opt<-optim(fn=maxConf,par=c(min(interval),max(interval)),minDiff=minDiff,level=0.99)
  bestValue<-Inf
  while(opt$value<bestValue){
   bestValue<-opt$value
   opt<-optim(fn=maxConf,par=c(floor(opt$par[1]),floor(opt$par[2])),minDiff=minDiff,level=0.99)
  }
  rect(data3[opt$par[1],2],0,data3[opt$par[2],2],0.04,col="limegreen",border="limegreen")
  estimates<-rbind(estimates,c(chr,data3[opt$par[1],2],data3[opt$par[2],2],minDiff,roughEst))
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
